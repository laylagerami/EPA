This folder contains the current results for the small subset of data. 
Each signature was run through CausalR and CARNIVAL
CARNIVAL results are pending...

Cytoscape_Nets.pptx -> Tutorial showing how to visualise networks 
Metadata -> Compound metadata e.g. sig id to pert id mapping
Network_Data -> Networks used in the analysis
Results -> CausalR -> Contains several things: 
1. network_enrichment_correlations.csv -> Correlation of each experiment based on their enrichment in the provided cell stress gene sets
2. network_jaccards.csv -> Jaccard coefficients of each experiment based on the nodes in their reconstructed networks
3. Folders for each experiment which contain:
	a. *_enrichment.csv -> Enrichment results for particular gene set
	b. corExplainedNodes.. -> Subnetwork outputs from CausalR. One subnetwork for each key driver protein (subnetwork encompasses all concordant interactions linking the key driver to the TFs)
	c. full_net.siff -> Concatenated subnetworks for each experiment into one subnetwork
	d. ResultsTable* -> Ranked list of nodes output from CausalR
4. Scripts -> Various scripts, most important is the Results_Analysis.Rmd which contains some preliminary analysis
5. Transcriptomics_Data -> Raw and processed LINCS data
