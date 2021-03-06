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
full_net = read.table("../Network_Data/omnipath_full_causalr.sif",sep="\t")
cons_net = read.table("../Network_Data/omnipath_hepg2_causalr.sif",sep="\t")

full_net_nodes =  unique(c(as.character(full_net$V1),as.character(full_net$V3)))
cons_net_nodes =  unique(c(as.character(cons_net$V1),as.character(cons_net$V3)))

# perform the clusterprofiler enrichment
library(clusterProfiler)

# first we need to change all gene symbols to entrez IDs
library(org.Hs.eg.db)
convert = function(symbols){
  entrez = AnnotationDbi::select(org.Hs.eg.db,
                                 keys = symbols,
                                 columns= c("SYMBOL","ENTREZID"),
                                 keytype="SYMBOL")
  ids = unique(as.character(entrez$ENTREZID))
  ids = ids[!is.na(ids)]
  return(ids)
}

# convert networks
full_net_nodes_entrez = convert(full_net_nodes)
cons_net_nodes_entrez = convert(cons_net_nodes)

# convert output nodes
all_full_nodes_entrez = lapply(all_full_nodes, convert)
all_cons_nodes_entrez = lapply(all_cons_nodes, convert)

# do the enrichment together, and then separately
test <- compareCluster(geneCluster = all_full_nodes_entrez, fun = "enrichPathway",universe=full_net_nodes_entrez)
library(enrichplot)
dotplot(test)

test2 <- compareCluster(geneCluster = all_cons_nodes_entrez, fun = "enrichPathway",universe=cons_net_nodes_entrez)

# save
p = dotplot(test2)
p + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle("HEPG2 Constrained Network")
ggsave("../Results/CausalR/dotplot_constrained.png",width=11,height=7)

e = emapplot(test2)
ggsave("../Results/CausalR/emapplot_constrained.png",e,width=11,height=7)

p = dotplot(test)
p + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle("Full Network")
ggsave("../Results/CausalR/dotplot_full.png",width=11,height=7)
e = emapplot(test)
ggsave("../Results/CausalR/emapplot_full.png",e,width=11,height=7)
# separate enrichment, save as .txt file, then plot the interesting ones?
# full net
for(i in (1:length(all_full_nodes_entrez))){
  
  # extract condition
  condition = names(all_full_nodes_entrez[i])
  
  # get nodes
  nodes_entrez = all_full_nodes_entrez[[i]]
  
  # perform enrichment
  enrichment = enrichPathway(nodes_entrez,universe=full_net_nodes_entrez)
  
  # turn into table
  tab = data.frame(enrichment)
  
  # save pathway table
  fname = paste0("../Results/CausalR/per_condition_pathways/",condition,".txt")
  write.table(tab,fname,sep="\t",quote=F,row.names=F)
  
  # do plots
  # convert to gene symbol
  edox <- setReadable(enrichment, 'org.Hs.eg.db', 'ENTREZID')
  p = dotplot(edox,showCategory=20)
  fname = paste0("../Results/CausalR/per_condition_plots/",condition,"_dotplot.png")
  ggsave(fname,p,width=11,height=7)
  
  p2 <- cnetplot(edox, node_label="all")
  fname = paste0("../Results/CausalR/per_condition_plots/",condition,"_cnetplot.png")
  ggsave(fname,p2,width=11,height=7)
  
  p3 <- emapplot(edox)
  fname = paste0("../Results/CausalR/per_condition_plots/",condition,"_emapplot.png")
  ggsave(fname,p3,width=11,height=7)
}

#cons net
for(i in (1:length(all_cons_nodes_entrez))){
  
  # extract condition
  condition = names(all_cons_nodes_entrez[i])
  
  # get nodes
  nodes_entrez = all_cons_nodes_entrez[[i]]
  
  # perform enrichment
  enrichment = enrichPathway(nodes_entrez,universe=cons_net_nodes_entrez)
  
  # turn into table
  tab = data.frame(enrichment)
  
  # save pathway table
  fname = paste0("../Results/CausalR/per_condition_pathways/",condition,".txt")
  write.table(tab,fname,sep="\t",quote=F,row.names=F)
  
  # do plots
  # convert to gene symbol
  edox <- setReadable(enrichment, 'org.Hs.eg.db', 'ENTREZID')
  p = dotplot(edox,showCategory=20)
  fname = paste0("../Results/CausalR/per_condition_plots/",condition,"_dotplot.png")
  ggsave(fname,p,width=11,height=7)
  
  p2 <- cnetplot(edox, node_label="all")
  fname = paste0("../Results/CausalR/per_condition_plots/",condition,"_cnetplot.png")
  ggsave(fname,p2,width=11,height=7)
  
  p3 <- emapplot(edox)
  fname = paste0("../Results/CausalR/per_condition_plots/",condition,"_emapplot.png")
  ggsave(fname,p3,width=11,height=7)
}

# similarity between all networks (networks and pathways)
# concatenate the conditions
all_conds_nodes_entrez = c(all_full_nodes_entrez,all_cons_nodes_entrez)

# create a function
overlapcoefficient = function(i,j){
  oc = length(intersect(i,j))/min(length(i),length(j))
  return(oc)
}

# do the loop
rows = data.frame()
for(i in (1:length(all_conds_nodes_entrez))){
  for(j in (1:length(all_conds_nodes_entrez))){
    condi = names(all_conds_nodes_entrez[i])
    condj = names(all_conds_nodes_entrez[j])
    
    nodesi = all_conds_nodes_entrez[[i]]
    nodesj = all_conds_nodes_entrez[[j]]
    
    oc = overlapcoefficient(i=nodesi,j=nodesj)
    
    row = data.frame(t(data.frame(c(condi,condj,oc))))
    rows = rbind(rows,row)
  }
}

colnames(rows)= c("cond_1","cond_2","OC")
rows$OC = as.numeric(as.character(rows$OC))
ggplot(rows,aes(x=cond_1,y=cond_2,fill=OC)) +
  geom_tile()+
  geom_text(aes(label=round(OC,2)))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("../Results/CausalR/network_overlap_coefficient.png",width=8,height=6)
# time series/dose concordance?

# replicate concordance?