#this scrip takes all the output files from AffyAnalysis.R
#and calculate the summary for community analysuis
#process: 1) standardize the controls, 2) scale each experiment by the controls
#3)remove unresponsive probes

#get the parameters from the environment:
#Parameters:
#od = the output directory for files that are to be collected
#exptinfo= a tab delim file of experimental information
#outdir = the output directory
#organ = organ name of interest
#runID = run id value
#t1 = the quantile/percential value for threshold based on expression intensity
#t2 = the quantile/percential value for threshold based on max rank change
#k = the number of communities
#note fc= fc cutoff, so fc2 is converted using log(2,2) to become 1 in the code

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  break
  ##supply default values
  
  outdir <- '/datadrive/TG-Gates/Output/testLiv'
  chem<-'colchicine'
  organ<-'Liver'
  runID<-'testLiv'
  Repeat<-"Single"
  t1<-0.75
  t2<-0.75
  k<-50
  fc<- 2
}else{
  print(args)
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
#libraries
library(gdata)
library(doBy)
library(gplots)
library(RColorBrewer)
mypalette <- rev(brewer.pal(9,"RdBu"))
###################################
#read in the file with the meta data for the experiments, and file names
#to make it simple as these parameters are required for AffyAnalysis to create the file name
fname<-paste("RunInfo", organ, paste(runID, 'txt', sep='.'), sep='-')
#read in the file info
info<-read.table(file.path(outdir, fname), sep='\t', header=FALSE, na.strings='', quote="\"")
colnames(info)<-c("runID", "fileName","SPECIES","ORGAN_ID", "CHEMICAL", "DOSE_LEVEL", "SACRIFICE_PERIOD","TEST-TYPE","SINGLE_REPEAT_TYPE", "folder", "treatment", "experiment")
#adding in factorization of some levels
info$DOSE_LEVEL<-factor(info$DOSE_LEVEL, levels=c("Control","Low","Middle","High"), ordered=TRUE)
info$SACRIFICE_PERIOD<-factor(info$SACRIFICE_PERIOD, levels=c("3 hr", "6 hr", "9 hr", "24 hr", "4 day","8 day", "15 day", "29 day"), ordered=TRUE)

#output for verification:
info[1,]
#now get in all the output files from ArryAnalysis and assemble into one dataframe
files<-as.character(unique(info$fileName))
for(i in 1:length(files)){
  if(i == 1){
    dataTemp<-read.table(file.path(outdir,"ProcArray",files[i]), sep='\t', header=TRUE, na.strings='', quote="\"")
  }
  else{
    inFile<-read.table(file.path(outdir,"ProcArray",files[i]), sep='\t', header=TRUE, na.strings='', quote="\"")
    dataTemp<-cbind(dataTemp, inFile)
  }  
}
#changing class to facilitate the downstream steps
dataTemp<-as.matrix(dataTemp)
#extract the controls and so that we use only data with controls
#NOTE that any experiments that do NOT have a matched control will be excluded!!!!!
ctrlInfo<-subset(info, DOSE_LEVEL == 'Control')
#ctrls<-ctrlInfo$treatment
treatInfo<-subset(info, DOSE_LEVEL != 'Control')
#treats<-treatInfo$treatment
validExpts<-intersect(treatInfo$experiment, ctrlInfo$experiment)
infoNew<-subset(info, experiment %in% validExpts)
#reorder data according to experimental conditions etc
dataTemp<-dataTemp[,colnames(dataTemp) %in% infoNew$treatment]
ctrlData<-dataTemp[, colnames(dataTemp) %in% ctrlInfo$treatment]
###########
#write out the excluded experimental units
noCtrl<-setdiff(treatInfo$experiment, ctrlInfo$experiment)
if (length(noCtrl) >0){
  noCtrlData<-subset(info, experiment %in% noCtrl)
  #save output
  fname<-paste("NoControls", runID, organ, paste(gsub('-','',Sys.Date()), 'txt', sep='.'), sep='-')
  write.table(noCtrlData, file=file.path(outdir, fname), sep='\t', col.names=TRUE, quote=FALSE)
}
rm(ctrlInfo, treatInfo, noCtrlData,validExpts, inFile )
gc()
######################
#now to scale the data. Assumption is that the probe intensity distribtuion from
#control arrays is the same. using the mean value to scale
#figure out the replication factor for scaling
exptReps<-summaryBy(treatment~experiment, data=infoNew, FUN=length, keep.names=TRUE)
#summaryBy I believe does a sort, but just to be sure
exptReps<-exptReps[order(exptReps$experiment),]
#get order by EXPERIMENT so that the timing of the experiments are better ordered
#because the treatments are the column names
ctrlInfo<-subset(infoNew, DOSE_LEVEL == 'Control')
ctrlData<-dataTemp[,as.character(ctrlInfo[do.call(order, subset(ctrlInfo, select=(experiment))),]$treatment)]
orderAssays<-dataTemp[,as.character(infoNew[do.call(order, subset(infoNew, select=c(experiment))),]$treatment)]

#now to get the median of the control and scaling factor
ctrlMed<-apply(ctrlData,2,median)
ctrlGMed<-median(apply(ctrlData,2,as.numeric))
sf<-ctrlGMed/ctrlMed
fullsf<-rep(sf, exptReps$treatment)
#and scale the data
scaleArray<-sweep(orderAssays,2,FUN='*',fullsf)
dim(scaleArray)
dim(orderAssays)
summary(c(scaleArray))
summary(c(orderAssays))
scaleCtrl<-scaleArray[,as.character(ctrlInfo$treatment)]
summary(c(scaleCtrl))
#scaleCtrl<-scaleCtrl[,as.character(ctrlInfo[do.call(order, subset(ctrlInfo, select=c(CHEMICAL, DOSE_LEVEL, SACRIFICE_PERIOD))),]$treatment)]
#scaleArray<-scaleArray[,as.character(infoNew[do.call(order, subset(infoNew, select=c(CHEMICAL, DOSE_LEVEL, SACRIFICE_PERIOD))),]$treatment)]

###################
#plot
#write plot to file
pnm<-paste("DensityPlotv2",runID, paste(gsub('-','',Sys.Date()), 'png', sep='.'), sep='-')
png(file=file.path(outdir, pnm), width=7, height=7, units="in",pointsize=16, bg='white', res=150)
plot(density(orderAssays), main='density of log2 intensities')
lines(density(scaleArray), col=2)
lines(density(scaleCtrl), col=3)
legend("topright",legend=c('All data', 'Scaled data', 'Scaled controls'), col=c(1,2,3), lty=1)
dev.off()
##################################
#Data thresholding:
#remove the probes that are consistantly low expressors
#implying that they are non-responsive to any treatment
#then remove probes that are not changing over all experiments, incd control

#first threshold
#Does a probe ever have an intensity > t1 quantile intensity of all data
#asked a different way, if the probe always in the bottom t1 quantile for expression
quantT1<-scaleArray - quantile(scaleArray,t1)
quantT1[quantT1 <=0]<-NA
temp<-quantT1[!apply(is.na(quantT1),1,all),]
scaleArrayT1<-scaleArray[row.names(temp),]
dim(scaleArrayT1)
summary(c(scaleArrayT1))
#####
#remove probes that are nonresponsive to treatment (differentially expressed)
#get order by EXPERIMENT so that the timing of the experiments are better ordered
#because the treatments are the column names
#note the order should have been preserved, this is just insurance :)
ctrlData<-scaleArrayT1[,as.character(ctrlInfo[do.call(order, subset(ctrlInfo, select=(experiment))),]$treatment)]
orderAssays<-scaleArrayT1[,as.character(infoNew[do.call(order, subset(infoNew, select=(experiment))),]$treatment)]
#And getting full control array:
#repeat each scaling factor according to the number of treatments in an experiment
#doing so buy getting column names
tmp<-rep(1:ncol(ctrlData), times=exptReps$treatment)
#now getting the corresponding columns (columns are repeated)
fullCtrl<-ctrlData[,tmp]
######################
#and calculate the DE#
######################
#in log2 space so subtract
deArray<-orderAssays - fullCtrl
#remove those that are not DE based on (log2 fc <1.5)
det2<-deArray
det2[abs(deArray) < log(fc,2)] <-NA
temp<-det2[!apply(is.na(det2),1,all),]
scaleArrayT1FC<-orderAssays[row.names(temp),]
scaleArrayT1FC<-as.matrix(scaleArrayT1FC[,as.character(infoNew[do.call(order, subset(infoNew, select=c(CHEMICAL, DOSE_LEVEL, SACRIFICE_PERIOD))),]$treatment)])
#
deArray<-deArray[,as.character(infoNew[do.call(order, subset(infoNew, select=c(CHEMICAL, DOSE_LEVEL, SACRIFICE_PERIOD))),]$treatment)]
deData<-as.matrix(deArray[row.names(temp),colnames(deArray) %in% subset(infoNew, DOSE_LEVEL != 'Control')$treatment])
#save output
fname<-paste("ProcessedDatav2", paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),organ, paste(runID, 'txt', sep='.'), sep='-')
write.table(scaleArrayT1FC, file=file.path(outdir, fname), sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE) # Writes expression values to text file in working directory.
#
fname<-paste("DiffExprDatav2", paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),organ, paste(runID, 'txt', sep='.'), sep='-')
write.table(deData, file=file.path(outdir, fname), sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE) # Writes expression values to text file in working directory.
#generates the files needed for downstream analusis (updated 01/2015)

#descritize the data to decrease the search space
deDesc<-deData
#using a cutoff of 2 fold change
deDesc[deDesc <= -1]<--1
deDesc[deDesc > -1]<-0
deDesc[deData >= 1]<-1

fname<-paste("Desc2FC",paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),organ, paste(runID, 'txt', sep='.'), sep='-')
write.table(deDesc, file=file.path(outdir, fname), sep='\t', col.names=TRUE, row.names=TRUE,quote=FALSE )

#
###########
#plotting
pnm<-paste("DataThresholdPlotv2",paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),runID, paste(gsub('-','',Sys.Date()), 'png', sep='.'), sep='-')
png(file=file.path(outdir, pnm), width=7, height=7, units="in",pointsize=16, bg='white', res=150)
plot(density(scaleArray), main="Data Processing", col="black", xlab='log2 intensity')
lines(density(scaleArrayT1), col="red")
lines(density(scaleArrayT1FC),col='blue')
legend("topright", legend=c('Inital data', paste("Post T1=",t1), paste("Post FC=",fc)), col=c('black','red','blue'), lty=1)
dev.off()

########
mapkey<-as.data.frame(cbind(key=c(1:length(colnames(deData))), label=colnames(deData)))
pnmk<-paste("heatmapKey",t1, fc,runID, paste(gsub('-','',Sys.Date()), 'txt', sep='.'), sep='-')
write.table(mapkey, file=file.path(outdir, pnmk), sep='\t', col.names=TRUE, quote=FALSE)

dcor<-cor(t(deData), method='pearson')
pctree<-hclust(as.dist(1-dcor), method='complete')
pnm<-paste("DiffExprplot",runID, paste(gsub('-','',Sys.Date()), 'png', sep='.'), sep='-')
png(file=file.path(outdir, pnm), width=8, height=8, units="in",pointsize=16, bg='white', res=150)
heatmap.2(deData, dendrogram='row', Rowv=as.dendrogram(pctree), trace='none', col=rev(brewer.pal(9,"RdBu")), 
          symbreaks=TRUE,labRow='', labCol=mapkey$key, main=bquote("DiffExpr,"~T1==.(t1)~FC==.(fc)))
dev.off()
AssayTree<-hclust(dist(t(scaleArrayT1FC), method="euclidean"), method="ward")
##################
#plotting intensity heatmap
Acor<-cor(t(scaleArrayT1FC), method='pearson')
ProbeTree<-hclust(as.dist(1-Acor), method='complete')
pnm<-paste("IntensityHeatmap",paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),runID, paste(gsub('-','',Sys.Date()), 'png', sep='.'), sep='-')
png(file=file.path(outdir, pnm), width=12, height=8, pointsize=14,units='in', bg='white', res=150)
heatmap.2(as.matrix(scaleArrayT1FC), Rowv=as.dendrogram(ProbeTree),Colv='', col=mypalette[5:9], dendrogram='row', trace='none', labRow='',labCol='', main=bquote("Scaled Intensities,"~T1==.(t1)~FC==.(fc)))
dev.off()




