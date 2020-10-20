# Script to run CARNIVAL
library(CARNIVAL)


assignPROGENyScores <- function (progeny = progeny, progenyMembers = progenyMembers, 
                                 id = "gene", access_idx = 1) 
{
  if (id == "uniprot") {
    idx <- which(names(progenyMembers) == "uniprot")
    progenyMembers <- progenyMembers[[idx]]
  }
  else {
    idx <- which(names(progenyMembers) == "gene")
    progenyMembers <- progenyMembers[[idx]]
  }
  members <- matrix(data = , nrow = 1, ncol = 2)
  pathways <- colnames(progeny)
  ctrl <- intersect(x = access_idx, y = 1:nrow(progeny))
  if (length(ctrl) == 0) {
    stop("The indeces you inserted do not correspond to \n              the number of rows/samples")
  }
  for (ii in 1:length(pathways)) {
    mm <- progenyMembers[[which(names(progenyMembers) == 
                                  pathways[ii])]]
    for (jj in 1:length(mm)) {
      members <- rbind(members, c(pathways[ii], mm[jj]))
    }
  }
  members <- members[-1, ]
  scores <- matrix(data = , nrow = nrow(progeny), ncol = nrow(members))
  colnames(scores) <- members[, 2]
  rownames(scores) <- rownames(progeny)
  members <- unique(members)
  for (i in 1:ncol(scores)) {
    for (j in 1:nrow(scores)) {
      scores[j, i] <- as.numeric(progeny[j, members[which(members[, 
                                                                  2] == colnames(scores)[i]), 1]])
    }
  }
  pxList <- list()
  for (ii in 1:length(access_idx)) {
    pxList[[length(pxList) + 1]] <- as.data.frame(t(as.matrix(scores[access_idx[ii], 
                                                                     ])))
  }
  names(pxList) <- rownames(progeny)[ctrl]
  return(pxList)
}


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
targets = t(data.frame(targets))
targets = data.frame(rbind(targets,rep(1,ncol(targets))))
colnames(targets) = as.character(unlist(targets[1,]))
targets = targets[-c(1),]
rownames(targets) = NULL

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
  root_dir = paste0("../Results/CARNIVAL/",cond)
  dir.create("../Results/CARNIVAL")
  dir.create(root_dir)
  
  # extract tf and progeny file
  load(file=system.file("progenyMembers.RData",package="CARNIVAL"))
  
  measObj = read.table(files[grepl("meas_",files)],header=T)
  weightObj = read.table(files[grepl("scores_",files)],header=T)

  progenylist = assignPROGENyScores(progeny = weightObj, 
                                    progenyMembers = progenyMembers, 
                                    id = "gene", 
                                    access_idx = 1)
  # run CARNIVAL (Inv and normal with both networks)
  
  #INV, FULL
  dir_name = paste0(root_dir,"/inv_full")
  dir.create(dir_name)
 # r1 = runCARNIVAL(inputObj=NULL, # no targets
#                   weightObj = progenylist$`1`,
#                   measObj = measObj, 
#                   netObj = netfilefull,
#                   solverPath = "../../../ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/cplex",
#                   threads = 10,
#                   solver="cplex", # need this for it to work
#                   dir_name=dir_name)
  
  #INV, CONS
  dir_name = paste0(root_dir,"/inv_cons")
  dir.create(dir_name)
 # r2 = runCARNIVAL(inputObj=NULL, # no targets
  #                 weightObj = progenylist$`1`,
   #                measObj = measObj, 
    #               netObj = netfilehepg2,
       #            threads = 10,
     ##              solverPath = "../../../ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/cplex",
    #               solver="cplex", # need this for it to work
    #               dir_name=dir_name)
  
  #NORMAL,FULL
  dir_name = paste0(root_dir,"/std_full")
  dir.create(dir_name)
  r3 = runCARNIVAL(inputObj=targets, 
                   weightObj = progenylist$`1`,
                   measObj = measObj, 
                   netObj = netfilefull,
                   solverPath = "../../../ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/cplex",
                   threads = 10,
                   solver="cplex", # need this for it to work
                   dir_name=dir_name)
  
  #NORMAL,CONS
  dir_name = paste0(root_dir,"/std_cons")
  dir.create(dir_name)
#  r4 = runCARNIVAL(inputObj=targets, 
  #                 weightObj = progenylist$`1`,
   #                measObj = measObj, 
   #                netObj = netfilehepg2,
   #                solverPath = "../../../ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/cplex",
      #             threads = 10,
      #             solver="cplex", # need this for it to work
      #             dir_name=dir_name)
}
