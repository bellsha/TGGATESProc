# TGGATESProc
Code for processing the Rat datafiles from the open TGGATES project

The array data are processed by experiment with the AffyAnalysis.R code. 
Differentially expressed probes are identified in the CombineArrayDE.R scripts.
ProbeAnnotation.R gets the probe level annotation for the microarrays that is used later on
TGGateslabprep.R and TGGatespathologyprep.R preprocesses the clinical chemistry/lab values and the pathology scores (respecitivly) to give a "hit" call needed for a sparse matrix.

gene2biospace.R https://github.com/bellsha/TGGATESProc/blob/master/gene2biospace.R is used to go from the dataframe of differentailly expressed probes (output of CombineArrayDE, https://github.com/bellsha/TGGATESProc/blob/master/CombineArrayDE.R), and using probe to entreze gene id mappings (ProbeAnnotation.R, https://github.com/bellsha/TGGATESProc/blob/master/ProbeAnnotation.R) and rat reactome pathway to uniprot protien mapping (ReactomeCalssv2.R, https://github.com/bellsha/Reactome2Network/blob/master/ReactomeClassv2.R) generate a table of enriched Reactome pathways for each chemical treatment, using the hypergeometric distribution https://github.com/bellsha/cpAOP/blob/master/HyperEnrich.R. Output file used in the networks: https://github.com/bellsha/cpAOP_Supplementary/blob/master/PathEnrichReactomeDataFrame.txt
