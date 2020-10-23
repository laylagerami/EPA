# analyse causal reasoning results
# also extract networks as .sif files with annotations for directionality and cuasl drivers
# for later e.g. Cytoscape visualisation
library(ReactomePA)

# load CausalR
causalr_dir = "../Results/CausalR"

# list the conditions
conditions = list.dirs(path=causalr_dir,recursive = F)

# loop over conditions, extract nodes, store them
per_cond_nodes = list()
conditions_etc = list()
for(cond in conditions){
  
  # extract experimental condition
  time_dose = strsplit(cond,"/")[[1]][4]
  
  # list output tpes
  output_types = list.dirs(path=cond,recursive = F) # i.e. full or constrained network
  
  # loop over output types
  for(output in output_types){
    
    # full or constrained net?
    network_type = strsplit(output,"/")[[1]][5]
    
    # now list .sif in this dir
    files = list.files(path=output,full.names = T,pattern="*.sif")
    files = files[!files=="all.sif"]
    
    # now list .txt in this dir
    txtfiles = list.files(path=output,full.names = T,pattern="*.anno.txt")
    txtfiles = txtfiles[!txtfiles=="node_annotations.txt"]
    txtfiles = txtfiles[!txtfiles=="main_nodes.txt"]
    
    # concatenate the .sif files
    all_files = list() # extract networks
    all_main_nodes = list() # extract driver
    for(sif in files){
      #get the main node from file name
      main_node = strsplit(strsplit(sif,"-")[[1]][5],".sif")[[1]][1]
      
      # omit directionality
      main_node = gsub(pattern="+",replacement="",main_node,fixed = T)
      main_node = gsub(pattern="-",replacement="",main_node,fixed=T)
      
      # read the net if possible
      net = try(read.table(sif,sep="\t"))
      if(!inherits(net,"try-error")){
        all_files[[sif]] = net
      }
      
      # add to list of nets and main nodes
      all_main_nodes[[sif]] = main_node 
    }
    
    # do the same with node annotations
    all_node_annotations = list()
    for(txt in txtfiles){
      # read the table
      tab = try(read.table(txt,sep="\t",skip = 1))
      
      # add to list of nodes annotations
      if(!inherits(tab,"try-error")){
        all_node_annotations[[txt]] = tab
      }
    }
    
    # make df with main nodes
    all_main_nodes = unname(all_main_nodes)
    main_nodes_df = data.frame(t(data.frame(all_main_nodes)))
    main_nodes_df$main = rep(1,nrow(main_nodes_df))
    rownames(main_nodes_df) = NULL
    
    # make full network file
    all_files = unname(all_files)
    all_sif = do.call("rbind",all_files)
    # remove duplicate rows from network
    deduped.data <- unique(all_sif)
    
    # make full anno
    all_files_anno = unname(all_node_annotations)
    all_files_anno = do.call("rbind",all_files_anno)
    # remove duplicate rows from network
    deduped.anno <- unique(all_files_anno)
    
    # save files
    fnet = paste0(output,"/all.sif")
    fnodes = paste0(output,"/main_nodes.txt")
    fstate = paste0(output,"/node_annotations.txt")
    
    write.table(deduped.data,fnet,sep="\t",quote=F,row.names = F)
    write.table(deduped.anno,fstate,sep="\t",quote=F,row.names=F)
    write.table(main_nodes_df,fnodes,sep="\t",quote=F,row.names=F)
    
    # now save the nodelist to a list
    all_nodes_in_net = unique(c(as.character(deduped.data$V1),as.character(deduped.data$V3)))
    
    # create a describing string with time,dose,network type
    data_desc = paste0(time_dose,"_",network_type)
    
    # save
    per_cond_nodes[[output]] = all_nodes_in_net
    conditions_etc[[output]] = data_desc
  }
}

# list of nodes for each condition
conditions_etc = unname(unlist(conditions_etc))
names(per_cond_nodes) = conditions_etc

# split into full and constrained for the network analysis
all_full_nodes = per_cond_nodes[grepl("full",names(per_cond_nodes))]
all_cons_nodes = per_cond_nodes[grepl("cons",names(per_cond_nodes))]

# get network nodes for full and constrained for enrichment background
