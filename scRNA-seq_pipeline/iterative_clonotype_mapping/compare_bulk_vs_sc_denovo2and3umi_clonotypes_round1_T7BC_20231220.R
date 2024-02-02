
library(tidyverse)
library(igraph)
library(stringr)
library(DescTools)
library(patchwork)

df_denovo2 <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/denovo_clonotypes_round1_2umi_table_BC_pairs_20231220.txt",
                        header=TRUE)

df_denovo2 <- df_denovo2 %>% transform(consensus_clone_id=paste0("sc2_",consensus_clone_id),
                                      modality_generation="single_cell2")

df_denovo3 <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/denovo_clonotypes_round1_3umi_table_BC_pairs_20231220.txt",
                         header=TRUE)
df_denovo3 <- df_denovo3 %>% transform(consensus_clone_id=paste0("sc3_",consensus_clone_id),
                                      modality_generation="single_cell3")

df_bulk <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq055SP_T7_shuffle_10X/cloneBCs/table_consensus_clones_w_BCs_20231203.txt",
                        header=TRUE)

df_bulk2 <- df_bulk %>% transform(BC_pair=paste0(real_barcode1,"-",real_barcode2),
                                  consensus_clone_id=paste0("b_",consensus_clone_id),
                                  modality_generation="bulk") %>%
  select(-real_barcode1,-real_barcode2)

# df_joined <- rbind(df_bulk2,df_denovo2,df_denovo3)
df_joined <- rbind(df_denovo2,df_denovo3)


##


clones_id <- df_joined %>% pull(consensus_clone_id) %>% unique()

df_pairs_clones <- CombPairs(clones_id)

for (id in seq(dim(df_pairs_clones)[1])){
  
  if ((id %% 100)==0){
    print(id)
  }
  
  BC_clone1 <- df_joined %>% filter(consensus_clone_id==df_pairs_clones$X1[id]) %>% pull(BC_pair) 
  BC_clone2 <- df_joined %>% filter(consensus_clone_id==df_pairs_clones$X2[id]) %>% pull(BC_pair) 
  
  jaccard_pair_BC <- length(intersect(BC_clone1,BC_clone2))/length(union(BC_clone1,BC_clone2))

  df_pairs_clones$jaccard_pair_BC[id] <- jaccard_pair_BC 

}

dev.new()
dev.set(4)
ggplot(df_pairs_clones) + stat_bin(aes(x=jaccard_pair_BC),color="firebrick",geom="step")+
  # stat_bin(aes(x=jaccard_pair_BC2),color="blue",geom="step")+
  scale_y_log10()

df_similar <- df_pairs_clones %>% filter(jaccard_pair_BC>0.5)
df_pairs_clones %>% filter(jaccard_pair_BC>0.05)

dim(df_similar)
df_similar %>% pull(X1) %>% unique %>% length
df_similar %>% pull(X2) %>% unique %>% length
df_similar %>% pull(jaccard_pair_BC) %>% median()
df_similar %>% pull(jaccard_pair_BC) %>% mean()

# heuristic: retain only the bulk if joined. 
df_joined_filtered <- df_joined %>% filter(!(consensus_clone_id %in% df_similar$X2))
df_joined_filtered$consensus_clone_id %>% unique() %>% length()

# write.table(df_joined_filtered,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/denovo2umi_w_bulk_joined_round1_clonotypes_BC_pairs_20231220.txt",
#                         row.names=TRUE,quote=FALSE,sep="\t")

write.table(df_joined_filtered,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/denovo3umi_w_bulk_joined_round1_clonotypes_BC_pairs_20231220.txt",
            row.names=TRUE,quote=FALSE,sep="\t")

# Convert the data frame to a matrix

jac_mat_BC_v2 <- matrix(data=NA,nrow=length(clones_id),ncol=length(clones_id))
rownames(jac_mat_BC_v2) <- clones_id
colnames(jac_mat_BC_v2) <- clones_id
for (id in seq(dim(df_pairs_clones)[1])){
  jac_mat_BC_v2[df_pairs_clones$X1[id],df_pairs_clones$X2[id]] <- df_pairs_clones$jaccard_pair_BC[id]
}


dist_thresh <- 0.1
conn_mat <- jac_mat_BC_v2>=dist_thresh
g  <- graph.adjacency(conn_mat)
clu <- components(g)

graph <- graph_from_adjacency_matrix(conn_mat, mode = 'undirected')
dev.set(4)
plot(graph,vertex.size=2,vertex.label.cex=0.5)


df_metad_denovo <- df_denovo2 %>% group_by(consensus_clone_id) %>% summarize(MOI=length(BC_pair)) %>% 
  transform(is_matched= (consensus_clone_id %in% df_similar$X2))

df_metad_bulk <- df_bulk2 %>% group_by(consensus_clone_id) %>% summarize(MOI=length(BC_pair)) %>% 
  transform(is_matched= (consensus_clone_id %in% df_similar$X1))

ggplot(df_metad_denovo) + stat_bin(aes(x=MOI,color=is_matched),geom="step",position="identity")+
  scale_x_log10()

plt1 <- ggplot(df_metad_denovo) + geom_step(aes(x=MOI,color=is_matched),stat="ecdf")+
  scale_x_log10()+labs(title="de novo")

plt2 <- ggplot(df_metad_bulk) + geom_step(aes(x=MOI,color=is_matched),stat="ecdf")+
  scale_x_log10()+labs(title="bulk")


plt1+plt2
