# Script to run CARNIVAL
library(CARNIVAL)

# Get Network
netfilehepg2 = read.csv("../Network_Data/omnipath_full_formatted_hepg2.csv")
netfilefull = read.csv("../Network_Data/omnipath_full_formatted.csv")

colnames(netfilefull) = c('source','Interaction','target')
colnames(netfilehepg2) = c('source','Interaction','target')


netfilehepg2$Interaction = ifelse(netfilehepg2$Interaction=="inhibits",-1,1)
netfilefull$Interaction = ifelse(netfilefull$Interaction=="inhibits",-1,1)

write.table(netfilefull,"../Network_Data/omnipath_full_carnival.sif",sep="\t",quote = F,row.names = F)
write.table(netfilehepg2,"../Network_Data/omnipath_hepg2_carnival.sif",sep="\t",quote = F,row.names = F)

netall = "../Network_Data/omnipath_full_carnival.sif"
nethepg2 = "../Network_Data/omnipath_hepg2_carnival.sif"

# get Targets
targets = as.character(read.csv("../Metadata/spiperone_targets_RH.txt",sep="\t")$Target)
targets = unlist(strsplit(targets,", "))

# List Dirs
dirs = list.dirs("../Transcriptomics_Data/tf_progeny",full.names = T,recursive = T)
dirs = dirs[2:8] # get rid of first dir (root dir)
dirs = dirs[[1]]
dir.create("../Results/Test")
# Iterate over Dirs
for(dir in dirs){
  # get both files
  files = list.files(dir,full.names = T)
  
  # get cond
  cond = strsplit(strsplit(files[1],"/")[[1]][4],"_measurements")[[1]][1]
  
  # extract tf and progeny file
  tf = read.table(files[grepl("meas_",files)],header = T)
  progeny = read.table(files[grepl("scores_",files)],header=T)
  
  # run CARNIVAL (Inv and normal with both networks)
  
  #INV, FULL
  r1 = runCARNIVAL(weightObj = progeny,
                   measObj = tf, 
                   netObj = netfilefull,
                   dir_name="../Results/Test")
  
  #INV, CONS
  
  #NORMAL,FULL
  
  #NORMAL,CONS
}
