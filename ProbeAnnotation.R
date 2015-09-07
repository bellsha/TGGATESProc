#getting the annotation, creating the table for enrichment
#must have the annotation from a cel file, any cell will do
#For finding out the annotation
#this might be something that is spit out in a file?
#library(affy)
#if you are unsure of the annotation you need, go into a direcroty with your affy 
#temp<-ReadAffy #be in a data containing directory
#annotation(temp)
#[1] "rat2302"
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("rat2302.db", "GO.db"))
#cant figure out how to automate this!!
library(rat2302.db)
#options for data to retreive
columns(rat2302.db)
#[1] "PROBEID"      "ENTREZID"     "PFAM"         "IPI"          "PROSITE"     
#[6] "ACCNUM"       "ALIAS"        "CHR"          "CHRLOC"       "CHRLOCEND"   
#[11] "ENZYME"       "PATH"         "PMID"         "REFSEQ"       "SYMBOL"      
#[16] "UNIGENE"      "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "GENENAME"    
#[21] "UNIPROT"      "GO"           "EVIDENCE"     "ONTOLOGY"     "GOALL"       
#[26] "EVIDENCEALL"  "ONTOLOGYALL" 
#this gets me the probe keys for extracting the other informaiton
#need the probes
probes<-keys(rat2302.db, keytype="PROBEID")
#the entrezis will allowe cross db referencing 
#I think that the go/ontollogy are the right ones.....other option is "all"
test<-select(rat2302.db, keys=probes, columns=c("ENTREZID","GO", "ONTOLOGY", "UNIPROT"), keytype="PROBEID")
fn<-paste("Files/rat2302.probe.entrez.go_", gsub('-','',Sys.Date()),".txt",sep="")
baseDir<-"../"
#write.table(test, file =file.path(baseDir, fn) , sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
#To get in the GO info
library(GO.db)
GO<-unique(select(GO.db, keys=test$GO, columns=c("TERM", "DEFINITION"), keytype="GOID"))
probeDesc<-merge(test, GO, by.x="GO", by.y="GOID", all.x=TRUE)
#to try and keep same ordering as befe
probeDesc<-probeDesc[,c(colnames(test), "TERM",'DEFINITION')]

write.table(probeDesc, file =file.path(baseDir, fn) , sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
