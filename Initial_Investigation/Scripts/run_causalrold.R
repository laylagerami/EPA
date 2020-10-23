library(CARNIVAL)
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

ccg <- CreateCCG("Ex2_network_SBV_Omnipath_causalr.sif")
#cell_line_folders <- list.dirs(recursive=FALSE)
#for(cell_line in cell_line_folders){
dir.create(file.path("RESULTS_CausalR"),showWarnings = FALSE)
compound_folders = list.dirs(recursive=FALSE) #path=cell_line,
print(compound_folders) 
  for(compound in compound_folders){
    drug = unlist(strsplit(unlist(strsplit(compound,"/"))[2],"_"))[1]
    print(drug)
    if(drug == "RESULTS"){
      next
    }
    dir.create(file.path(paste0("RESULTS_CausalR/",drug)))
    results_dir = paste0("RESULTS_CausalR/",drug)
    tf_activities = data.frame(t(read.table(list.files(path=compound,pattern="_50.txt",full.names=TRUE))))
    #tf_activities$gene_name = mapIds(org.Hs.eg.db,
    #                                 keys=as.character(tf_activities$X1),
    #                                 keytype="UNIPROT",
    #                                 column="SYMBOL")
    #tf_activities <- tf_activities %>% distinct(gene_name, .keep_all=TRUE) # remove duplicate gene names
    rownames(tf_activities) <- tf_activities$X1
    tf_activities <- as.matrix(tf_activities[,2,drop=FALSE]) # convert factor to numerical
    tf_activities[,1] <- ifelse(tf_activities[,1] > 0,1,-1)
    write.table(tf_activities,file=paste0(compound,"/tf_disc.txt"),col.names=FALSE,sep="\t",quote=FALSE)
    expData <- ReadExperimentalData(paste0(compound,"/tf_disc.txt"),ccg)
    
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


