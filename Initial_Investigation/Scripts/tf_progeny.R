# Script to get TF and PROGENy scores from a gene expression matrix

# 1. Script to prepare DoRoTHea TF activities

# Input: Matrix of gene expression data (Rows are genes (Entrez), columns are compounds)
# Output: 1. .txt file for each signature found in the input matrix
#         2. Folder for each compound (compoundname_measurements)
#         with TF activities (meas_50.txt) as UniProt ID

# Import packages

library(CARNIVAL)
library(org.Hs.eg.db)
library(foreach)
library(doParallel)
library(plyr)
library(data.table)

# Initialise cluster
n = 2 # change to number of cores needed
myCluster <- makeCluster(n, type="FORK",outfile="") 
registerDoParallel(myCluster)

# Load files for dorothea
file.copy(from=system.file("dorothea_TF_mapping.csv",package="CARNIVAL"),to=getwd(),overwrite=TRUE)
load(file = system.file("BEST_viperRegulon.rdata",package="CARNIVAL"))
map<-read.csv("dorothea_TF_mapping.csv")

#Open matrix 
gexfile = "../Transcriptomics_Data/processed_signatures/speripone_matrix.csv"
gex_df = fread(gexfile,header=TRUE,sep=",") # First row is python header
gex_df = as.data.frame(gex_df) # Change into df
compound_names = names(gex_df) # Column names are compounds/conditions
compound_names = compound_names[!compound_names %in% 'V1'] # Get rid of 'V1'

#Now change gene ids to gene symbols using metadata
gene_info = fread('gene_info.csv',header=TRUE) # Import metadata
gene_info = as.data.frame(gene_info) # Read as df

converted = merge(gex_df,gene_info,by.x='V1',by.y='pr_gene_id') # Map to gene symbol
converted_symbols = converted$pr_gene_symbol # Extract symbols 
gex_df$V1 = converted_symbols # Make row names into symbols

# Run
output <- 
  foreach(i = compound_names) %dopar% { # Loop over every compound
    print(i)
    fname = paste0("../Transcriptomics_Data/tf_progeny/",i,".txt") # Set the file name for the signature file
    compound_sig = cbind(gex_df[1],gex_df[i]) # Take the gene names and corresponding measurements
    write.table(compound_sig,fname,sep="\t",quote=F,row.names=F,col.names=T,append=T) # Write
    df = read.table(fname,sep="\t",header=T,row.names=1) # Read in the signature .txt file
    TF_genesymbol<-try( # run DoRothEA, confidence levels A, B and C
      runDoRothEA(df, regulon=viper_regulon, confidence_level=c('A','B','C')),
      silent = T
    )
    if(inherits(TF_genesymbol,"try-error")){ # If there is an error for some reason, skip the compound
      next
    }
    #TF_uniprot<-GeneSymbol2Uniprot(TF_genesymbol, map, 1, 2) # Map to UniProt
    folder = paste0("../Transcriptomics_Data/tf_progeny/",i,"_measurements/") # Set folder name
    generate_measfile(measurements=TF_genesymbol, topnumber=50, write2folder=folder) # Write TF activities to folder
  }
# Stop the cluster
stopCluster(myCluster)

# 2. Script to prepare PROGENy pathway scores

# Input: Matrix of gene expression data (Rows are genes (Entrez), columns are compounds), and
#        The .txt file for each signature (From prepare_input_parellel.R)
#         And a measurement folder for each compound (From prepare_input_parellel.R)
# Output: PROGEny pathway weights .txt in each compound's measurement folder

# Import packages
library(CARNIVAL)
library(org.Hs.eg.db)
library(foreach)
library(doParallel)
library(plyr)
library(data.table)

# set n to number of cores
n = 2
myCluster <- makeCluster(n, type="FORK",outfile="")
registerDoParallel(myCluster)

# Load files for progeny
file.copy(from=system.file("model_NatComm+14_human.csv",package="CARNIVAL"),to=getwd(),overwrite=TRUE)
weight_matrix<-read.csv("model_NatComm+14_human.csv")

#Open matrix 
gexfile = "../Transcriptomics_Data/processed_signatures/speripone_matrix.csv"
gex_df = fread(gexfile,header=TRUE,sep=",") # First row is python header
gex_df = as.data.frame(gex_df) # Change into df
compound_names = names(gex_df) # Column names are compounds/conditions
compound_names = compound_names[!compound_names %in% 'V1'] # Get rid of 'V1'

# test for e.g. 10 compounds
# compound_names = compound_names[1:10]

output <-
  foreach(i = compound_names) %dopar% {  # loop over each compound
    print(i)
    fname <- paste0("../Transcriptomics_Data/tf_progeny/",i, ".txt") # get signature .txt
    df = read.table(fname,sep="\t",header=TRUE,row.names=1) # Read back in file
    df_genenames <- data.frame('gene'=rownames(df),df) # make df with rownames = gene symbols
    
    #Run progeny 
    pathway_scores <- try(
      runPROGENy(df_genenames,weight_matrix,z_scores=F),
      silent = T)
    if(inherits(pathway_scores,"try-error")){ # if it fails then skip
      next
    }
    
    #Generate input files
    folder = paste0("../Transcriptomics_Data/tf_progeny/",i,"_measurements/scores_") # get folder name
    scores <- rbind(rownames(pathway_scores),pathway_scores[,1]) # put into correct format
    write.table(scores,paste0(folder,i,".txt"),col.names=F,row.names=F,quote=F,sep='\t') # save
    
  }




stopCluster(myCluster)


