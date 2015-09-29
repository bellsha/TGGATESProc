#formatting the attributes file
data<-read.table("/wikihomedir/sbell/TG-Gates/Data/Open-tggates_AllAttribute.txt",sep="\t", header=TRUE, na.strings='NA', quote="\"")
#fixes to chemicals (removing spaces and hypens)
chems<-gsub("[[:space:]:]", '_', data$COMPOUND_NAME)
chems<-gsub("['()]", '', chems)

#chems<-gsub(' ','_',data$COMPOUND_NAME)
#chems<-gsub('-','',chems)
#bc TNFa is weird
chems<-gsub("TNF.*", 'TNFa', chems)
#to make consistent with the label in treatment
#chems<-gsub("[[:space:],%-/]","", chems)
#
treatment<-gsub("[[:space:],%-/]","",as.factor(paste(data$ORGAN_ID, chems,data$DOSE_LEVEL, data$SACRI_PERIOD, data$SIN_REP_TYPE, sep="_")))
#treatment<-paste('Liver', chems, dataRL2$DOSE_LEVEL, gsub(' ','', dataRL2$SACRI_PERIOD), dataRL2$SIN_REP_TYPE, sep='_')
experiment<-gsub("[[:space:],%-/]","",as.factor(paste(data$ORGAN_ID, chems, data$SACRI_PERIOD, data$SIN_REP_TYPE, sep="_")))
labData<-cbind(treatment, experiment,chems, data)
dataRL<-subset(labData, ORGAN_ID =='Liver' & SPECIES == 'Rat' & TEST_TYPE =='in vivo')
labLiver<-dataRL[,c(1,2,21,25,26,30:67)]
write.table(labLiver, file="/wikihomedir/sbell/TG-Gates/Data/Lab/LiverLabData201411.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
#
#now to start processing the data- need to account for scores that are different than controls
library(doBy)
#first order the data b/c it will be reordered by doby functions
labLiverO<-labLiver[do.call(order, subset(labLiver, select=c(experiment, treatment))),]
#calculate the treatment-wise mean
labM<-summaryBy(.~treatment, data=subset(labLiverO, select=-c(experiment, SACRI_PERIOD)), FUN=mean,na.rm=TRUE, keep.names=TRUE)
#calculate the upper and lower limits for the whole dataset, by sacrifice period
#this is bc there are bigger changes that happen wrt age of animal
#allowing for up 10% of the data to be FP
ctrlQ1<-summaryBy(.~SACRI_PERIOD, data=labLiverO[,3:ncol(labLiverO)], FUN=quantile,probs=0.025, na.rm=TRUE, keep.names=TRUE)
ctrlQ3<-summaryBy(.~SACRI_PERIOD, data=labLiverO[,3:ncol(labLiverO)], FUN=quantile,probs=0.975, na.rm=TRUE, keep.names=TRUE)

upperLM<-merge(unique(labLiverO[,c(1,3)]), ctrlQ3, all.x=TRUE)
upperLM2<-upperLM[do.call(order, subset(upperLM, select=treatment)),]    
lowerLM<-merge(unique(labLiverO[,c(1,3)]), ctrlQ1, all.x=TRUE)
lowerLM2<-lowerLM[do.call(order, subset(lowerLM, select=treatment)),]    
#
labD<-labM
labD[,2:ncol(labD)][labM[,2:ncol(labM)] > upperLM2[,3:ncol(upperLM2)]] <-1
labD[,2:ncol(labD)][lowerLM2[,3:ncol(lowerLM2)] <= labM[,2:ncol(labM)] & labM[,2:ncol(labM)] <= upperLM2[,3:ncol(upperLM2)]] <-0
labD[,2:ncol(labD)][labM[,2:ncol(labM)] < lowerLM2[,3:ncol(lowerLM2)]] <--1

#as a test see if any controls have a value
ctrl<-subset(dataRL, DOSE_LEVEL =='Control')$treatment
x<-subset(labD, treatment %in% ctrl)
summary(as.factor(c(as.matrix(labD[,2:ncol(labD)]))))
#-1      0      1 
#2088 183553   2563 
summary(as.factor(c(as.matrix(x[,2:ncol(x)]))))
#-1     0     1 
#148 47121   203 
#remove Rows/Columns which are empty
labD<-labD[!apply(labD[,2:ncol(labD)],1,function(x) sum(abs(x), na.rm=TRUE)==0),]
labD<-labD[,!apply(labD[,2:ncol(labD)],2,function(x) sum(abs(x), na.rm=TRUE)==0)]
write.table(labD, file="/wikihomedir/sbell/TG-Gates/Data/Lab/LiverLabDescData201411.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)

################
#and chemcial->labscore edges
liv<-t(subset(labD, select =-treatment))
liv[is.na(liv)]<-0
colnames(liv)<-labD$treatment
n<-matrix(nrow=nrow(liv), ncol=ncol(liv))
for(i in 1:nrow(liv)){
  #note for negative relationships it will just maintain a negative sign
  n[i,]<-gsub('1',rownames(liv)[i], liv[i,])
}

n[n==0]<-NA
n2<-as.data.frame(n)
lcdat<-NULL
lcnam<-colnames(liv)
#for(i in 1:2){
for(i in 1:ncol(n2)){
  tmp<-na.omit(n2[,i])
  tmp2<-cbind(Treat=rep(lcnam[i], length(tmp)), Out=as.character(tmp))
  lcdat<-rbind(lcdat, tmp2)
}
lcdat<-as.data.frame(lcdat)
c2l<-merge(unique(dataRL[,c(1,3,16,21,24)]), unique(lcdat), by.x='treatment',by.y='Treat', all.y=TRUE)
test<-gsub('-.*','Decrease', c2l$Out)
test2<-gsub('^(?!Decrease).*','Increase', test, perl=TRUE)
c2l$Change<-test2
c2l$Out<-gsub('-', '', c2l$Out)
write.table(c2l, file="/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverLabChemEdged201411.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)


