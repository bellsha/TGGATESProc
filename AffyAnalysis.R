#Affy array analysis
#performs RMA normalization for all arrays in a folder then
#Calcualtes the average probe intensity for each treatment
#Parameters:
#d= the working directory
#od = the output directory for files that are to be collected
#exptinfo= a tab delim file of experimental information
#outdir = the output directory
#chem = the chemical name of interest
#organ = organ name of interest
#runID = run id value
#NOTE This code is modified for use with TG-Gates data. 
#Changes are wrt the input file for experimental parameters

################################
#get the parameters from the environment:
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  break
  ##supply default values
  d <-"/TG-Gates/Data/Rat/in_vivo/Liver/Single/colchicine.Rat.in_vivo.Liver.Single/celfiles"
  exptinfo<-'/TG-Gates/Data/Microarray_Data_RatLabels_20141030.txt'
  outdir <- '/TG-Gates/Output/test'
  chem<-'colchicine'
  organ<-'Liver'
  runID<-'test'
  Repeat<-"Single"
  
}else{
  print(args)
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
#Import libraries
library(affy)
#read in experimental paramters for all the microarray experiments
allinfo<-read.table(exptinfo, sep='\t', header=TRUE, na.strings='', quote="\"", stringsAsFactors=FALSE)
#
#Move to the folder containing the data of interest
setwd(d)
print(getwd())

#analysis for affy data with rma normalization
data<-ReadAffy()
eset <- rma(data,verbose=FALSE, background=FALSE )#rma with no bkgd correction
#NOTE: eset from rma returns obj with expression values in log2 space
#create a file name and save the eset object
fname<-paste("eset", organ, chem, paste(runID, 'txt', sep='.'), sep='-')
write.exprs(eset, file=fname) # Writes expression values to working directory.
###############################################
#get the average for each treatment level including controls
#convert to dataframe
edata <- as.data.frame(exprs(eset))
#note the columns are the cel file names
#use this to remove the information for the arrays from the information table
info<-subset(allinfo, BARCODE %in% as.character(gsub('.CEL','',colnames(edata)) ))
info$Array_ID<-paste(info$BARCODE, '.CEL', sep='')
#note that to get the control and all treatments corresponding to the chemical treatment group,
#I need to use the chem variable and not the CHEMICAL column. 
#A value of 'NA" for dose indicates controls
treatment<-as.factor(paste(info$ORGAN_ID, chem, info$DOSE_LEVEL, info$SACRIFICE_PERIOD, info$SINGLE_REPEAT_TYPE, sep="_"))
#remove special characters from treatment as it will become column name
treatment<-gsub("[[:space:],%-/]","",treatment)
info$treatment<-treatment
#adding in experiment information as it is needed downstream and easier to put in here
experiment<-gsub("[[:space:],%-/]","",as.factor(paste(info$ORGAN_ID, chem, info$SACRIFICE_PERIOD, info$SINGLE_REPEAT_TYPE, sep="_")))
info$experiment<-experiment
treats<-unique(treatment)
ArrayAve<-NULL
for(i in 1:length(treats)){
  grp_arrays<-subset(info, treatment == treats[i])$Array_ID
  #updating to catch occurances where samples are missing
  if(length(grp_arrays) <2){
  	grp_ave<-rep(NA, nrow(edata))
  }
  else{
	grp_ave<-apply(edata[,colnames(edata) %in% grp_arrays],1,FUN=mean, na.rm=TRUE)
  }
  ArrayAve<-cbind(ArrayAve, grp_ave)
}
#addind in the identifiers, which correspond to the treatment parameters
colnames(ArrayAve)<-treats
#removal of columns lacking data
ArrayAve<-ArrayAve[,!apply(is.na(ArrayAve),2,all)]
#now write the file for use in integrating datasets
Repeat<-unique(info$SINGLE_REPEAT_TYPE)
fname<-paste("ExpAVE",organ,chem,Repeat, paste(runID, 'txt', sep='.'), sep='-')
#write object to table to facilitate input and mergeing
#the pattern of filenames should facilitate the retreival and then the checking
#write.table(ArrayAve, file=file.path(outdir, fname), sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE) # Writes expression values to text file in working directory.
write.table(ArrayAve, file=file.path(outdir, "ProcArray", fname), sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE) 
#####
#need to make a new outputs file for the next steps- standardizing across the experiments
#and other things that may pop up where sorting based on column name might be a PITA
#inclusion of runID for check when wanting to collect several days worth post-run
info2<-cbind(runID,fileName=fname,unique(info[,c("SPECIES","ORGAN_ID", "COMPOUND_NAME", "DOSE_LEVEL", "SACRIFICE_PERIOD","TEST_TYPE","SINGLE_REPEAT_TYPE", "folder", "treatment", "experiment")]))
info2<-subset(info2, treatment %in% colnames(ArrayAve))
fname<-paste("RunInfo", organ, paste(runID, 'txt', sep='.'), sep='-')
#append the existing table (or make new if yet to be created)
#need to remove the column headings so that it appends clean
colnames(info2)
#[1] "runID"              "fileName"           "SPECIES"            "ORGAN_ID"           "COMPOUND_NAME"     
#[6] "DOSE_LEVEL"         "SACRIFICE_PERIOD"   "TEST_TYPE"          "SINGLE_REPEAT_TYPE" "folder"            
#[11] "treatment"          "experiment" 
write.table(info2, file=file.path(outdir, fname), sep='\t', append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)

