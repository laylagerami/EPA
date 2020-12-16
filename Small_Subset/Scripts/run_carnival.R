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

netfilefull = read.csv("../Network_Data/omnipath_full_formatted.csv")

colnames(netfilefull) = c('source','Interaction','target')


netfilefull$Interaction = ifelse(netfilefull$Interaction=="inhibits",-1,1)

write.table(netfilefull,"../Network_Data/omnipath_full_carnival.sif",sep="\t",quote = F,row.names = F)

netall = "../Network_Data/omnipath_full_carnival.sif"


# List Dirs
dirs = list.dirs("../Transcriptomics_Data/tf_progeny",full.names = T,recursive = T)
dirs = dirs[2:50] # get rid of first dir (root dir)

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
  r1 = runCARNIVAL(inputObj=NULL, # no targets
                   weightObj = progenylist$`1`,
                   measObj = measObj, 
                   netObj = netfilefull,
                   solverPath = "../../../ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/cplex",
                   threads = 10,
                   solver="cplex", # need this for it to work
                   dir_name=root_dir)
  
}
