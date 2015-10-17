# TGGATESProc
Code for processing the Rat datafiles from the open TGGATES project

The array data are processed by experiment with the AffyAnalysis.R code. 
Differentially expressed probes are identified in the CombineArrayDE.R scripts.
ProbeAnnotation.R gets the probe level annotation for the microarrays that is used later on
TGGateslabprep.R and TGGatespathologyprep.R preprocesses the clinical chemistry/lab values and the pathology scores (respecitivly) to give a "hit" call needed for a sparse matrix.
