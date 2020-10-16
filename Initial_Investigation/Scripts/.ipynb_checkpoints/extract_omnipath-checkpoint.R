library(OmnipathR)

# download
interactions <- import_Omnipath_Interactions()
interactions_directed = subset(interactions,is_directed==1)
interactions_signed = subset(interactions,consensus_stimulation==1 | consensus_inhibition==1)

# convert to igraph
op_g <- interaction_graph(interactions = interactions_signed) #4919 nodes, 20,292 edges

# constrain network
hepg2_prots = read.csv("../Proteomics_Data/E-PROT-20-query-results.tsv",sep="\t")
hepg2_prot_symbols = as.character(hepg2_prots$Gene.Name)

interactions_signed_hepg2 = subset(interactions_signed,source_genesymbol %in% hepg2_prot_symbols | target_genesymbol %in% hepg2_prot_symbols)
op_g_hepg2 <- interaction_graph(interactions = interactions_signed_hepg2) #4087 nodes, 14,259 edges

# save graphs
write.csv(interactions_signed,"../Network_Data/omnipath_full.csv",row.names = F,quote = F)
write.csv(interactions_signed_hepg2,"../Network_Data/omnipath_hepg2.csv",row.names = F, quote=F)

# now convert to proper direction etc
all_rows = list()
for(row in 1:nrow(interactions_signed)){
  rowdat = interactions_signed[row,]
  if(rowdat$consensus_direction==1){
    source = rowdat$source_genesymbol
    target = rowdat$target_genesymbol
  }else{ # if the consensus_direction = 1 then the source and target are flipped round
    source = rowdat$target_genesymbol
    target = rowdat$source_genesymbol
  }
  if(rowdat$consensus_stimulation==1){
    interaction = "activates"
  }else{
    interaction = "inhibits"
  }
  newrow = t(data.frame(c(source,interaction,target)))
  all_rows[[row]] = newrow
}

all_omnipath_formatted = do.call("rbind",all_rows)
rownames(all_omnipath_formatted) = NULL
all_omnipath_formatted = data.frame(all_omnipath_formatted)
colnames(all_omnipath_formatted) = c("source","interaction","target")

write.csv(all_omnipath_formatted,"../Network_Data/omnipath_full_formatted.csv",row.names = F,quote = F)

# do the same for constrained network
# now convert to proper direction etc
all_rows = list()
for(row in 1:nrow(interactions_signed_hepg2)){
  rowdat = interactions_signed[row,]
  if(rowdat$consensus_direction==1){
    source = rowdat$source_genesymbol
    target = rowdat$target_genesymbol
  }else{ # if the consensus_direction = 1 then the source and target are flipped round
    source = rowdat$target_genesymbol
    target = rowdat$source_genesymbol
  }
  if(rowdat$consensus_stimulation==1){
    interaction = "activates"
  }else{
    interaction = "inhibits"
  }
  newrow = t(data.frame(c(source,interaction,target)))
  all_rows[[row]] = newrow
}

all_omnipath_formatted = do.call("rbind",all_rows)
rownames(all_omnipath_formatted) = NULL
all_omnipath_formatted = data.frame(all_omnipath_formatted)
colnames(all_omnipath_formatted) = c("source","interaction","target")

write.csv(all_omnipath_formatted,"../Network_Data/omnipath_full_formatted_hepg2.csv",row.names = F,quote = F)
