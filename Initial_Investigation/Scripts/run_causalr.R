# Script to run CausalR
library(CausalR)
library(org.Hs.eg.db)
library(dplyr)

rank_hypotheses <- function(ccg,expdata,delta) {
  hypothesis <- RankTheHypotheses(ccg, expdata,delta,correctPredictionsThreshold=1,doParallel=TRUE,numCores=32,writeFile=FALSE)
  hypothesis <- hypothesis[complete.cases(hypothesis), ]
  hypothesis <- data.frame(uniprot = row.names(hypothesis), hypothesis)
  #hypothesis(row.names) <- gsymbol
  return(hypothesis[!(hypothesis$'p.value'>0.05),])
  #return(hypothesis)
}

rank_hypotheses_5 <- function(ccg,expdata,delta,results_dir) {
  hypothesis <- RankTheHypotheses(ccg, expdata,delta,correctPredictionsThreshold=1,doParallel=TRUE,numCores=32,writeFile=TRUE,outputDir=results_dir)
  hypothesis <- hypothesis[complete.cases(hypothesis), ]
  hypothesis <- data.frame(uniprot = row.names(hypothesis), hypothesis)
  #hypothesis(row.names) <- gsymbol
  return(hypothesis[!(hypothesis$'p.value'>0.05),])
  #return(hypothesis)
}

# Function to run calcs
run_causalr = function(tf_file,ccg,results_dir){
 
  expData <- ReadExperimentalData(tf_file,ccg)
  
  run1 <- rank_hypotheses(ccg,expData,1)
  
  colnames(run1)[3:9] <- paste("PL1_", colnames(run1[,c(3:9)]), sep = "")
  
  run2 <- rank_hypotheses(ccg,expData,2)
  
  colnames(run2)[3:9] <- paste("PL2_", colnames(run2[,c(3:9)]), sep = "")
  
  run3 <- rank_hypotheses(ccg,expData,3)
  
  colnames(run3)[3:9] <- paste("PL3_", colnames(run3[,c(3:9)]), sep = "")
  
  run4 <- rank_hypotheses(ccg,expData,4)
  
  colnames(run4)[3:9] <- paste("PL4_", colnames(run4[,c(3:9)]), sep = "")
  
  run5 <- rank_hypotheses_5(ccg,expData,5,results_dir)
  
  colnames(run5)[3:9] <- paste("PL5_", colnames(run5[,c(3:9)]), sep = "")
  
  big_run <- setNames(data.frame(matrix(ncol = 30, nrow = 0)), c("PL1_Score","PL1_Correct","PL1_Incorrect","PL1_Ambiguous","PL1_p.value","PL1_Enrichment.p.value","PL2_Score","PL2_Correct","PL2_Incorrect","PL2_Ambiguous","PL2_p.value","PL2_Enrichment.p.value","PL3_Score","PL3_Correct","PL3_Incorrect","PL3_Ambiguous","PL3_p.value","PL3_Enrichment.p.value","PL4_Score","PL4_Correct","PL4_Incorrect","PL4_Ambiguous","PL4_p.value","PL4_Enrichment.p.value","PL5_Score","PL5_Correct","PL5_Incorrect","PL5_Ambiguous","PL5_p.value","PL5_Enrichment.p.value"))
  
  all_sig_hyp_2 <- c(row.names(run1),row.names(run2),row.names(run3),row.names(run4),row.names(run5))
  all_sig_hyp_2 <- unique(all_sig_hyp_2)
  for (hyp in all_sig_hyp_2){
    big_run[hyp,"PL1_Score"]<-run1[hyp,"PL1_Score"]
    big_run[hyp,"PL1_Correct"]<-run1[hyp,"PL1_Correct"]
    big_run[hyp,"PL1_Incorrect"]<-run1[hyp,"PL1_Incorrect"]
    big_run[hyp,"PL1_Ambiguous"]<-run1[hyp,"PL1_Ambiguous"]
    big_run[hyp,"PL1_p.value"]<-run1[hyp,"PL1_p.value"]
    big_run[hyp,"PL1_Enrichment.p.value"]<-run1[hyp,"PL1_Enrichment.p.value"]
    big_run[hyp,"PL2_Score"]<-run2[hyp,"PL2_Score"]
    big_run[hyp,"PL2_Correct"]<-run2[hyp,"PL2_Correct"]
    big_run[hyp,"PL2_Incorrect"]<-run2[hyp,"PL2_Incorrect"]
    big_run[hyp,"PL2_Ambiguous"]<-run2[hyp,"PL2_Ambiguous"]
    big_run[hyp,"PL2_p.value"]<-run2[hyp,"PL2_p.value"]
    big_run[hyp,"PL2_Enrichment.p.value"]<-run2[hyp,"PL2_Enrichment.p.value"]
    big_run[hyp,"PL3_Score"]<-run3[hyp,"PL3_Score"]
    big_run[hyp,"PL3_Correct"]<-run3[hyp,"PL3_Correct"]
    big_run[hyp,"PL3_Incorrect"]<-run3[hyp,"PL3_Incorrect"]
    big_run[hyp,"PL3_Ambiguous"]<-run3[hyp,"PL3_Ambiguous"]
    big_run[hyp,"PL3_p.value"]<-run3[hyp,"PL3_p.value"]
    big_run[hyp,"PL3_Enrichment.p.value"]<-run3[hyp,"PL3_Enrichment.p.value"]
    big_run[hyp,"PL4_Score"]<-run4[hyp,"PL4_Score"]
    big_run[hyp,"PL4_Correct"]<-run4[hyp,"PL4_Correct"]
    big_run[hyp,"PL4_Incorrect"]<-run4[hyp,"PL4_Incorrect"]
    big_run[hyp,"PL4_Ambiguous"]<-run4[hyp,"PL4_Ambiguous"]
    big_run[hyp,"PL4_p.value"]<-run4[hyp,"PL4_p.value"]
    big_run[hyp,"PL4_Enrichment.p.value"]<-run4[hyp,"PL4_Enrichment.p.value"]
    big_run[hyp,"PL5_Score"]<-run5[hyp,"PL5_Score"]
    big_run[hyp,"PL5_Correct"]<-run5[hyp,"PL5_Correct"]
    big_run[hyp,"PL5_Incorrect"]<-run5[hyp,"PL5_Incorrect"]
    big_run[hyp,"PL5_Ambiguous"]<-run5[hyp,"PL5_Ambiguous"]
    big_run[hyp,"PL5_p.value"]<-run5[hyp,"PL5_p.value"]
    big_run[hyp,"PL5_Enrichment.p.value"]<-run5[hyp,"PL5_Enrichment.p.value"]
    
    k <- 0
    if (!is.na(big_run[hyp,"PL1_Score"])){
      k <- k + 1
    }
    if (!is.na(big_run[hyp,"PL2_Score"])){
      k <- k + 1
    }
    if (!is.na(big_run[hyp,"PL3_Score"])){
      k <- k + 1
    }
    if (!is.na(big_run[hyp,"PL4_Score"])){
      k <- k + 1
    }
    if (!is.na(big_run[hyp,"PL5_Score"])){
      k <- k + 1
    }
    big_run[hyp,"Total"] <- k
  }
  
  
  #Write files
  top_nodes <- big_run[which(big_run$Total == max(big_run$Total)), ]
  biggest_pl <- max(big_run$Total)
  if(biggest_pl == 0){
    print("NO")
  }
  list_top_nodes <- rownames(top_nodes)
  
  for(node in list_top_nodes){
    if (grepl('+',node,fixed=TRUE) == TRUE){
      protein <- unlist(strsplit(node,'+',fixed=TRUE))
      sign <- '+1'
    } else { 
      protein <- unlist(strsplit(node,'-',fixed=TRUE))
      sign <- '-1'
    }
    WriteExplainedNodesToSifFile(outputDir=results_dir,protein,sign,ccg,expData,delta=biggest_pl,correctlyExplainedOnly = TRUE)
  }
}

# Get Network
netfilehepg2 = read.csv("../Network_Data/omnipath_full_formatted_hepg2.csv")
netfilefull = read.csv("../Network_Data/omnipath_full_formatted.csv")

colnames(netfilefull) = c('source','Interaction','target')
colnames(netfilehepg2) = c('source','Interaction','target')

netfilehepg2$Interaction = ifelse(netfilehepg2$Interaction=="inhibits","Inhibits","Activates")
netfilefull$Interaction = ifelse(netfilefull$Interaction=="inhibits","Inhibits","Activates")

write.table(netfilefull,"../Network_Data/omnipath_full_causalr.sif",sep="\t",quote = F,row.names = F,col.names = F)
write.table(netfilehepg2,"../Network_Data/omnipath_hepg2_causalr.sif",sep="\t",quote = F,row.names = F,col.names = F)

netall = CreateCCG("../Network_Data/omnipath_full_causalr.sif")
nethepg2 = CreateCCG("../Network_Data/omnipath_hepg2_causalr.sif")

# List Dirs
dirs = list.dirs("../Transcriptomics_Data/tf_progeny",full.names = T,recursive = T)
dirs = dirs[2:8] # get rid of first dir (root dir)
dirs = dirs[[1]]
# Iterate over Dirs
for(dir in dirs){
  # get both files
  files = list.files(dir,full.names = T)

  # get cond
  cond = strsplit(strsplit(files[1],"/")[[1]][4],"_measurements")[[1]][1]
  root_dir = paste0("../Results/CausalR/",cond)
  dir.create("../Results/CausalR")
  dir.create(root_dir)
  
  # extract tf 
  tf_activities = data.frame(t(read.table(files[grepl("meas_",files)],header=T)))
  tf_activities <- as.matrix(tf_activities[,1,drop=FALSE]) # convert factor to numerical
  tf_activities[,1] <- ifelse(tf_activities[,1] > 0,1,-1)
  tf_file = paste0("../Transcriptomics_Data/tf_progeny/",cond,"_tf_disc.txt")
  write.table(tf_activities,file=tf_file,col.names=FALSE,sep="\t",quote=FALSE)
  
  #FULL
  dir_name = paste0(root_dir,"/full")
  dir.create(dir_name)
  run_causalr(tf_file=tf_file,ccg=netall,results_dir=dir_name)
  
  # CONS
  dir_name = paste0(root_dir,"/cons")
  dir.create(dir_name)
  run_causalr(tf_file=tf_file,ccg=nethepg2,results_dir=dir_name)
 
}
