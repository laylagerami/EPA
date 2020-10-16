# script to perform exploratory analysis of the data

library(RColorBrewer)
library(gplots)
library(ggplot2)
# first on gene-expression level
# read in data
gene_data = read.csv("../Transcriptomics_Data/processed_signatures/speripone_matrix.csv",row.names = 1)
gene_matrix = as.matrix(gene_data)

# do initial heatmap to see rough clustering 
png(file="../Results/Exploratory_Plots/gene_heatmap.png",res=300,width = 350, height = 225, units='mm')
ht = Heatmap(gene_matrix,show_row_names=F)
draw(ht)
dev.off()

# now do pairwise correlation
library(plyr)
library(tidyr)
d2 <- gene_data %>% 
  as.matrix %>%
  cor %>%
  as.data.frame %>%
  rownames_to_column(var = 'var1') %>%
  gather(var2, value, -var1)

# remove duplicates
d2 %>%
  mutate(var_order = paste(var1, var2) %>%
           strsplit(split = ' ') %>%
           map_chr( ~ sort(.x) %>% 
                      paste(collapse = ' '))) %>%
  mutate(cnt = 1) %>%
  group_by(var_order) %>%
  mutate(cumsum = cumsum(cnt)) %>%
  filter(cumsum != 2) %>%
  ungroup %>%
  select(-var_order, -cnt, -cumsum)

ggplot(d2, aes(var1, var2, fill= value)) + 
  geom_tile()+
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(name = "Correlation",
                      low = "#FFFFFF",
                      high = "#012345")
ggsave("../Results/Exploratory_Plots/gene_corr_mat.png",width=8,height=5)

# now on the TF-level
# read in data
tf_data = list.files(path="../Transcriptomics_Data/tf_progeny/",recursive = T,pattern="meas",full.names = T)

# loop read
list_of_tf_scores = list()
list_of_tfs = list()
names = c("24h_10uM","24h_3uM","6h_10uM_1","6h_10uM_2","6h_10uM_3","6h_10uM","6h_3uM")
for(i in 1:length(tf_data)){
  f = tf_data[[i]]
  read_file = data.frame(t(read.table(f,sep="\t")))
  colnames(read_file) = c("X1",names[[i]])
  read_file[,2] = as.numeric(as.character(read_file[,2]))
  tfs = as.character(read_file[,1])
  list_of_tfs[[i]] = tfs
  list_of_tf_scores[[i]] = read_file
}
names(list_of_tf_scores) = names
names(list_of_tfs) = names

# merge tf scores
all_scores = Reduce(function(...) merge(..., by="X1", all=TRUE), list_of_tf_scores)
tfs = all_scores$X1
all_scores$X1 = NULL
all_scores <- sapply( all_scores, as.numeric )
all_scores = data.frame(all_scores)
rownames(all_scores) = tfs
all_scores = as.matrix(all_scores)
M1 <- as(all_scores, "dgCMatrix")
library(qlcMatrix) 
corSparse(M1, Y = NULL, cov = FALSE)

