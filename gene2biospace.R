#this is ment to get the gene expression (changes) based on chemcial
#as the dataset requires descretization it will be based on array/treatment (treatment is easier)
#and dose_level='Control' will not be included in final results as all the values will be zero
#libraries
library(doBy)
library(igraph);library(arules)
#get the differential experss data
deDesc<-read.table("../Output/NetRun1501/Desc2FC-T1-0.75-FC-2-Liver-NetRun1501.txt", sep='\t', header=TRUE) #DE data from https://github.com/bellsha/TGGATESProc/blob/master/CombineArrayDE.R
deDA<-abs(deDesc)
#############################33
#get the experimental mappings
probe<-read.table("../Files/rat2302.probe.entrez.go_20150515.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=FALSE) #from: https://github.com/bellsha/TGGATESProc/blob/master/ProbeAnnotation.R
colnames(probe)<-c("ProbeID","ENTREZID","GOID", "Evi","GOprocess", "UniprotID", "GOTerm", "GODef")
tmp2<-subset(probe, ProbeID %in% row.names(deDesc))
path<-read.table("../Files/Reactome/ReactomePathways2UniProtRAT.txt", sep='\t',comment.char='', na.strings='', quote="\"", header=TRUE) #from: https://github.com/bellsha/Reactome2Network/blob/master/ReactomeClassv2.R
#remove the NA
path<-subset(path, UniProtID != "NA")
#fix the whitespace issue
########################
# Abstracting from probe to "biological space"
#there are a few different levels to explore.
#Merging everything to keep that into the options availabl
################
pathMap<-merge(tmp2[,c(1:3,6)], path[,c(1,2,4,7,9)], by.x='UniprotID', by.y='UniProtID')  

############################################
#lets build the pathway enrichment, with each row corresponding to a "community" 
m<-matrix(nrow=nrow(deDA), ncol=ncol(deDA))
for(i in 1:nrow(deDA)){
	m[i,]<-gsub('1',rownames(deDA)[i], deDA[i,])
}
m[m==0]<-NA
m2<-as.data.frame(m)
cdat<-NULL
cnam<-colnames(deDA)
for(i in 1:ncol(m2)){
	tmp<-na.omit(m2[,i])
	tmp2<-cbind(Treat=rep(cnam[i], length(tmp)), Probe=as.character(tmp))
	cdat<-rbind(cdat, tmp2)
}
cdat<-as.data.frame(cdat)
#Calculate the number of de probes for a given treatment
#Note that we only care about the probes for which we have pathway information
cdat2<-subset(cdat, Probe %in% unique(pathMap$ProbeID))
sizes<-summaryBy(Probe~Treat,data=cdat2, FUN=length, keep.names=TRUE) 
cdata<-merge(cdat2, sizes, by.x="Treat", by.y="Treat")
colnames(cdata)<-c('Treat','Probe','Size')
#merging using uniprotID as the identifier to go back to pathMap
cdataM<-merge(cdata, unique(pathMap[,1:2]), by.x="Probe", by.y="ProbeID")
sdataM<-na.omit(cdataM)
###########
#for reactome
source("../Code/HyperEnrich.R")

#by reactome pathway...order is 1st column=clusterID and second column=featureID and thrid column= size
#need to go by id
anno<-ClusReport(sdataM[,c(2,4,3)], background=sdataM$UniprotID, RefData=pathMap[,c(5,1)])
#anno<-ClusReport(sdataM[,c(2,4,3)], background=sdataM$UniprotID, RefData=pathMap[,c(6,1)])
#now make a dataframe for the items
pwdf<-matrix(nrow=length(unique(anno$CLID)), ncol=length(unique(anno$Label)))
ID<-unique(anno$CLID)
labID<-unique(anno$Label)
colnames(pwdf)<-labID
rownames(pwdf)<-ID

for(i in 1:length(ID)){
	tmp<-subset(anno, CLID ==ID[i])
	pwdf[i,]<-labID %in% tmp$Label
}
pwdf2<-as.data.frame(apply(pwdf,2,as.numeric))
pwdf2$treatment<-as.factor(row.names(pwdf))
#adding in the name to the anno file (bc enrichment function doesnt have it)
anno2<-merge(anno, unique(path[,1:2]), by.x="Label", by.y="ReactomeID", all.x=TRUE)
#all the enriched pathways
write.table(anno2,file="../Output/NetRun1501/Update/PathEnrichReactomeAnnoTable.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )
write.table(pwdf2,file="../Output/NetRun1501/Update/PathEnrichReactomeDataFrame.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )




