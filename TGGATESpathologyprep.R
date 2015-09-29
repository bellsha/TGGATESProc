#Processing of pathology data for TG Gates
library(doBy)
#pathology data is kinda a mess- I am assuming 5 slides per treatment
slides<-read.table("/wikihomedir/sbell/TG-Gates/Data/Pathology/open_tggates_pathology.csv", sep=',', header=TRUE, na.strings='', quote="\"")
#minor changes to make consistent with other data
#path lacks single/repeat type...
repInfo<-cbind(SACRIFICE_PERIOD=c('3 hr', '6 hr', '9 hr','24 hr', '4 day', '8 day', '15 day','29 day'), SIN_REP_TYPE = c('Single','Single','Single','Single','Repeat','Repeat','Repeat','Repeat'))
slides$chems<-gsub("[[:space:]:]", '_', slides$COMPOUND_NAME)
slides$chems<-gsub("['()]", '', slides$chems)
slidesL<-subset(slides, ORGAN == 'Liver')
slidesL<-merge(repInfo, slidesL, by.x='SACRIFICE_PERIOD', by.y='SACRIFICE_PERIOD')
experiment<-gsub("[[:space:],%-/]","",as.factor(paste(slidesL$ORGAN, slidesL$chems, slidesL$SACRIFICE_PERIOD, slidesL$SIN_REP_TYPE, sep="_")))
treatment<-gsub("[[:space:],%-/]","",as.factor(paste(slidesL$ORGAN, slidesL$chems,slidesL$DOSE_LEVEL, slidesL$SACRIFICE_PERIOD, slidesL$SIN_REP_TYPE, sep="_")))
slidesL<-cbind(experiment, treatment, slidesL)
slidesLO<-slidesL[do.call(order, subset(slidesL, select=c(experiment, treatment))),]
score<-as.character(slidesLO$GRADE_TYPE)
score[score=='P']<-1
score[score=='minimal']<-2
score[score=='slight']<-3
score[score=='moderate']<-4
score[score=='severe']<-5
slidesLO$score<-as.numeric(score)
slidesLO$indID<-paste(slidesLO$EXP_ID,slidesLO$GROUP_ID, slidesLO$INDIVIDUAL_ID, sep='.')
#testing max number of slides per treatment
#is should be 5 based on the data in TG gates (a maximum of 5 individuals/expt)
tmp<-unique(slidesLO[,c(1:11,16,18)])
x<-summaryBy(INDIVIDUAL_ID ~ EXP_ID + GROUP_ID, data=tmp, FUN=length)
summary(x[,3])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   2.000   2.754   5.000   5.000 
####################################################
#transform the data frame, there will be warnings
s2<-reshape(na.omit(slidesLO[,c('indID','FINDING_TYPE', 'score')]), idvar="indID",timevar="FINDING_TYPE", direction="wide")
colnames(s2)<-gsub('score.','', colnames(s2))
colnames(s2)<-gsub('[^[:alnum:]]','.', colnames(s2))
s3<-merge(unique(slidesLO[,c('indID','treatment','experiment')]), s2, by.x='indID', by.y='indID', all.y=TRUE)
s3O<-subset(s3[do.call(order, subset(s3, select=c(experiment, treatment))),],select=-indID)
#calculate the treatment-wise mean based on 5 reps
spS<-summaryBy(.~experiment+treatment, data=s3O, FUN=sum,na.rm=TRUE, keep.names=TRUE)
spM<-cbind(spS[,1:2],spS[,3:ncol(spS)]/5)
#descretization, based on a score of 1
pathDesc<-spM
pathDesc[,3:67][pathDesc[,3:67] <= 1]<-NA
pathDesc[,3:67][pathDesc[,3:67] > 1]<-1
pathDesc<-pathDesc[!apply(is.na(pathDesc[,3:67]),1,all),]
pathDesc<-pathDesc[,!apply(is.na(pathDesc),2,all)]
write.table(pathDesc, file="/wikihomedir/sbell/TG-Gates/Data/Pathology/open_tggates_pathologyDESC.txt", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(spM, file="/wikihomedir/sbell/TG-Gates/Data/Pathology/open_tggates_pathologySUM.txt", sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

##############3
#and chemcial->patholgy edges
liv<-t(subset(pathDesc, select =-c(treatment, experiment)))
colnames(liv)<-pathDesc$treatment
n<-matrix(nrow=nrow(liv), ncol=ncol(liv))
for(i in 1:nrow(liv)){
  n[i,]<-gsub('1',rownames(liv)[i], liv[i,])
}
n[n==0]<-NA
n2<-as.data.frame(n)
lcdat<-NULL
lcnam<-colnames(liv)
for(i in 1:ncol(n2)){
  tmp<-na.omit(n2[,i])
  tmp2<-cbind(Treat=rep(lcnam[i], length(tmp)), Out=as.character(tmp))
  lcdat<-rbind(lcdat, tmp2)
}
lcdat<-as.data.frame(lcdat)
treat2pheno<-merge(unique(slidesLO[,c('chems','SACRIFICE_PERIOD','DOSE_LEVEL','treatment','SIN_REP_TYPE')]),unique(lcdat), by.x='treatment', by.y='Treat')
write.table(treat2pheno, file="/wikihomedir/sbell/TG-Gates/Output/NetRun1410/LiverPathologyChemicalEdges201411.txt", sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE )


