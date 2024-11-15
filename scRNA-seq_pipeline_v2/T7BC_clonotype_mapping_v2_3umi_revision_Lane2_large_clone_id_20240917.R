
library(Seurat)
library(Matrix)
library(tidyverse)

df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling/valid_T7BC_Lane2_20240907.txt.gz",
                            header=TRUE)


set.seed(123)

# update the BC set
BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split(regex("_att.__"))
# BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split("_CS")

df_T7BC_valid$BC_first <- BC_all %>% lapply("[[",1) %>% unlist()
df_T7BC_valid$BC_second <- BC_all %>% lapply("[[",2) %>% unlist()
df_T7BC_valid$first_BC_CS <- ifelse(str_detect(df_T7BC_valid$BC_first,"CS1"),"CS1","CS2")
df_T7BC_valid$CS1_BC <- NA
df_T7BC_valid$CS1_BC[df_T7BC_valid$first_BC_CS=="CS1"] <- df_T7BC_valid$BC_first[df_T7BC_valid$first_BC_CS=="CS1"]
df_T7BC_valid$CS1_BC[df_T7BC_valid$first_BC_CS=="CS2"] <- df_T7BC_valid$BC_second[df_T7BC_valid$first_BC_CS=="CS2"]

df_T7BC_valid$CS2_BC <- NA
df_T7BC_valid$CS2_BC[df_T7BC_valid$first_BC_CS=="CS2"] <- df_T7BC_valid$BC_first[df_T7BC_valid$first_BC_CS=="CS2"]
df_T7BC_valid$CS2_BC[df_T7BC_valid$first_BC_CS=="CS1"] <- df_T7BC_valid$BC_second[df_T7BC_valid$first_BC_CS=="CS1"]

df_T7BC_valid2 <- df_T7BC_valid %>% select(-BC,-BC_first,-BC_second) %>% transform(BC_pair=paste0(CS1_BC,"__",CS2_BC)) %>% select(-CS1_BC,-CS2_BC)

# df_T7BC_valid$capture_direction <- paste0("CS",BC_all %>% lapply("[[",2) %>% unlist())
# df_T7BC_valid <- df_T7BC_valid %>% select(-BC)

df_T7BC_valid3 <- df_T7BC_valid2 %>% group_by(cBC_pdT,BC_pair) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC))
df_T7BC_valid3 <- df_T7BC_valid3 %>% filter(UMIs_BC>=3)

df_BC_UMI_per_cell <- df_T7BC_valid3 %>% group_by(cBC_pdT) %>% summarize(sum_T7_UMI_per_cell=sum(UMIs_BC),MOI=length(BC_pair))
df_BC_pairs_sum_UMI <- df_T7BC_valid3 %>% group_by(BC_pair) %>% summarize(sum_UMI=sum(UMIs_BC))

thresh_T7umi_per_cell <- 10
thresh_T7_per_BC_pair <- 10
dev.new()
dev.set(4)
ggplot(df_BC_UMI_per_cell) + stat_bin(aes(x=sum_T7_UMI_per_cell),geom="step")+
  scale_x_log10()+scale_y_log10()+
  geom_vline(xintercept=thresh_T7umi_per_cell,color="red",linetype="dashed")

ggplot(df_BC_pairs_sum_UMI) + stat_bin(aes(x=sum_UMI),geom="step")+
  scale_x_log10()+scale_y_log10()+
  geom_vline(xintercept=thresh_T7_per_BC_pair,color="red",linetype="dashed")

# TAATACAGACACGCCTCCGG_CS1__TTTGTTGACCAGATCGCTAA_CS2
#
#CACTACGGCGCCCTAGTACT-CS1--ACGTTCGTCTGAGTCAAATG-CS2                36
#GGTTACCCAGACCAATAAAA-CS1--GTATAGCATTAGTCATCAAC-CS2                27
# df_BC_pairs_sum_UMI %>% filter(BC_pair=="TAATACAGACACGCCTCCGG_CS1__TTTGTTGACCAGATCGCTAA_CS2")
# df_BC_pairs_sum_UMI %>% filter(BC_pair=="CACTACGGCGCCCTAGTACT_CS1__ACGTTCGTCTGAGTCAAATG_CS2")
# df_BC_pairs_sum_UMI %>% filter(BC_pair=="GGTTACCCAGACCAATAAAA_CS1__GTATAGCATTAGTCATCAAC_CS2")
# df_BC_pairs_sum_UMI %>% filter(sum_UMI>1.5e4)
# df_BC_pairs_sum_UMI %>% arrange(desc(sum_UMI)) %>% head(n=10)



hi_BC_pairs <- df_BC_pairs_sum_UMI %>% filter(sum_UMI>thresh_T7_per_BC_pair ) %>% pull(BC_pair)
high_cells <- df_BC_UMI_per_cell %>% filter(sum_T7_UMI_per_cell>thresh_T7umi_per_cell) %>% pull(cBC_pdT)

# # iterative exclusion of high-represented cells:
# hi_BC_pairs2 <- hi_BC_pairs[!(hi_BC_pairs %in% str_replace_all(high_rep_BCs_excluded_v2,"-","_"))]
# high_cells2 <- high_cells[!(high_cells %in% high_rep_cells_excluded_v2)]

length(hi_BC_pairs)
length(high_cells)



# df_T7BC_valid_hi <- df_T7BC_valid2 %>% filter(BC_pair %in% hi_BC_pairs) %>% 
#   filter(cBC_pdT %in% high_cells)
df_T7BC_valid_hi <- df_T7BC_valid3 %>% filter(BC_pair %in% hi_BC_pairs) %>% 
  filter(cBC_pdT %in% high_cells)


# generate T7BC matrix on filtered cells and BC pairs: 
cBC_list <- df_T7BC_valid_hi %>% pull(cBC_pdT) %>% unique()
cBC_list <- data.frame(cBC_id=seq(1:length(cBC_list)),cBC_pdT=cBC_list)
T7BC_list <- df_T7BC_valid_hi %>% pull(BC_pair) %>% unique()
T7BC_list <- data.frame(T7BC_id=seq(1:length(T7BC_list)),BC_pair=T7BC_list)



# cBC_list_2 <- cBC_list %>% filter(!(cBC_pdT %in% high_rep_cells_excluded_v2))
# T7BC_list_2 <- T7BC_list %>% filter(!(BC_pair %in% str_replace_all(high_rep_BCs_excluded_v2,"-","_")))

df_T7BC_valid_hi2 <-  df_T7BC_valid_hi %>% left_join(cBC_list) %>% left_join(T7BC_list)
# df_T7BC_valid_hi2 <-  df_T7BC_valid_hi %>% left_join(cBC_list_2) %>% left_join(T7BC_list_2)

T7BC_mat <- sparseMatrix(j=df_T7BC_valid_hi2$T7BC_id,
                         i=df_T7BC_valid_hi2$cBC_id,
                         x=df_T7BC_valid_hi2$UMIs_BC)
rownames(T7BC_mat) <- cBC_list$cBC_pdT
colnames(T7BC_mat) <- T7BC_list$BC_pair

BC_obj <- CreateSeuratObject(counts = t(T7BC_mat), min.cells = 0, min.features = 0)


# dimensional reduction
all_BCs <- rownames(BC_obj)
top_pc <- 100
n_neighbors <- 10
cluster_res <- 1
umap_NN <- 10

BC_obj_sub <- BC_obj
BC_obj_sub <- NormalizeData(BC_obj_sub, normalization.method = "RC", scale.factor = 10000)
BC_obj_sub <- FindVariableFeatures(BC_obj_sub, selection.method = "vst", nfeatures = length(all_BCs), verbose = TRUE)
BC_obj_sub <- ScaleData(BC_obj_sub, features = all_BCs)

BC_obj_sub <- RunPCA(BC_obj_sub, 
                     features = VariableFeatures(object = BC_obj_sub), 
                     verbose = FALSE, 
                     npcs = top_pc)

BC_obj_sub <- FindNeighbors(BC_obj_sub, dims = 1:top_pc, k.param = n_neighbors)
BC_obj_sub <- FindClusters(BC_obj_sub, resolution = cluster_res) 
BC_obj_sub <- RunUMAP(BC_obj_sub, dims = 1:top_pc/2, n.neighbors = umap_NN, seed.use = 42) 



dev.new()
dev.set(4)
DimPlot(BC_obj_sub,
        reduction = "umap", label=TRUE,
        label.box = FALSE,label.size =1.5,
        pt.size =0.5)+theme(legend.position="none")#+labs(title=sample_cell_oi)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# studying count distribution in the T7 BC per clone
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

seurat_clusters <- 0:64
labels <- c(1E-4,1E-3,1E-2,1E-1,1)
lims <- c(1E-4,1)

UMI_thresh <- 1
thresh_frac_cBC <- 0.1
thresh_min_cBC <- 3
FC_cBC_frac_thresh <- 1.5

list_core_T7BC <- vector(mode = "list", length = length(seurat_clusters))
names(list_core_T7BC) <- seurat_clusters
max_frac_per_cluster <- vector(mode= "numeric", length=length(seurat_clusters))
names(max_frac_per_cluster) <- seurat_clusters

# cluster_oi <- 35
dir_out <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling/cluster_UMAP_w_distribution_T7BC_representation_per_cell_3umi"
for (cluster_oi in seurat_clusters){
  
  
  print(cluster_oi)
  seur_obj_cluster_oi <- subset(BC_obj_sub, subset =( seurat_clusters==cluster_oi) )
  
  count_mat_oi <- seur_obj_cluster_oi@assays$RNA$counts
  
  # sum_UMI_BC <- sort(rowSums(count_mat_oi),decreasing = TRUE)
  
  # for each pBC, determine fraction of cBC with >10 UMI to that pBC
  n_cell_cluster <- dim(count_mat_oi)[2]
  frac_cBC <- (rowSums(count_mat_oi>UMI_thresh,na.rm=TRUE))/n_cell_cluster
  
  # Addition: mean UMI
  mean_UMI_per_cell <- rowMeans(count_mat_oi)
  dev.set(4)
  plt_pdf <- ggplot()+stat_bin(aes(x=mean_UMI_per_cell),geom="step")+scale_x_log10()+scale_y_log10()
  plt_cdf <- ggplot()+geom_step(aes(x=mean_UMI_per_cell,y=log(1-..y..)),stat="ecdf")+scale_x_log10()
  
  
  # finding per cluster threshold
  sorted_frac_cBC <- sort(frac_cBC,decreasing=TRUE)
  FC_frac_cBC <- sorted_frac_cBC[1:(length(sorted_frac_cBC)-1)]/sorted_frac_cBC[2:length(sorted_frac_cBC)]
  
  frac_cBC_cut <- sorted_frac_cBC[min(which(FC_frac_cBC>FC_cBC_frac_thresh))]
  
  max_frac_per_cluster_oi <- max(frac_cBC)
  max_frac_per_cluster[[as.character(cluster_oi)]] <- max_frac_per_cluster_oi
  
  ecdf_frac_cBC <- ecdf(frac_cBC)
  xplot <- (1:length(frac_cBC))/length(frac_cBC)
  
  final_thresh_frac_cBC <- max(thresh_min_cBC/n_cell_cluster,thresh_frac_cBC,frac_cBC_cut)
  
  cluster_core_eBC <- names(which(frac_cBC>=final_thresh_frac_cBC))
  list_core_T7BC[[as.character(cluster_oi)]] <- cluster_core_eBC
  # 
  core_pBC_plot <- ggplot()+geom_step(aes(x=xplot,y=1-ecdf_frac_cBC(xplot)))+
    geom_vline(xintercept=final_thresh_frac_cBC,color="red",linetype="dashed")+
    labs(y="1-CDF (top % eBC)",
         x=sprintf("Fraction of cell BCs in cluster with >%d UMI to T7-BC",UMI_thresh),title=sprintf("Cluster %d",cluster_oi))+
    scale_x_log10(breaks = labels,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lims)+
    scale_y_log10(breaks = labels,
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = lims)+annotation_logticks()
  
  cBC_oi <- Cells(BC_obj_sub)[BC_obj_sub$seurat_clusters==cluster_oi]
  Dim_plt <- DimPlot(BC_obj_sub, reduction = "umap", cells.highlight = cBC_oi,
                     sizes.highlight=0.3)+ggtitle(sprintf("cluster %s",cluster_oi))
  
  # Dim_plt+core_pBC_plot
  
  # core_pBC_plot
  fig_name <- sprintf("%s/clusters_T7BCs_id_iteration_%d_20240917.pdf",
                      dir_out,
                      cluster_oi)
  pdf(fig_name,width=10,height=10)
  print(Dim_plt+core_pBC_plot+plt_pdf+plt_cdf)
  dev.off()
  # 
  
}


# highly represented clones:
clusters_high_representation <- c(0,1)
high_rep_cells_excluded_v2 <- BC_obj_sub@meta.data %>% filter(seurat_clusters %in% clusters_high_representation) #%>% transform( rownames()
high_rep_BCs_excluded_v2 <- rbind(data.frame(clone_id="0",BC=list_core_T7BC[c("0")] %>% unlist()),data.frame(clone_id="1",BC=list_core_T7BC[c("1")] %>% unlist()))

write.table(high_rep_cells_excluded_v2 %>% transform(cellBC=rownames(high_rep_cells_excluded_v2)),"/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling/large_clones_cellBC_Lane2_20240917.txt",
            sep="\t",row.names=FALSE,quote=FALSE)
write.table(high_rep_BCs_excluded_v2,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling/large_clones_shuffleBC_Lane2_20240917.txt",
            sep="\t",row.names=FALSE,quote=FALSE)
# 
# # threshold some base quality on cluster
# dev.set(4)
# ggplot()+stat_bin(aes(x=max_frac_per_cluster),bins=20)+
#   geom_vline(xintercept=0.5,color="red")
# 
# MOI_per_cluster <- lengths(list_core_T7BC) 
# 
# df_MOI_per_cluster <- data.frame(MOI_per_cluster) 
# df_MOI_per_cluster <- df_MOI_per_cluster %>%  transform(seurat_clusters=rownames(df_MOI_per_cluster))
# 
# cells_per_cluster <- BC_obj_sub@meta.data %>% group_by(seurat_clusters) %>% summarize(cells_per_cluster=length(nCount_RNA)) %>% data.frame() %>% left_join(df_MOI_per_cluster)
# 
# ggplot()+stat_bin(aes(x=MOI_per_cluster),bins=20)+
#   geom_vline(xintercept=7,color="red")+scale_x_log10()
# 
# 
# # # manual curation based on bimodality of the map
# list_core_T7BC2 <- list_core_T7BC[max_frac_per_cluster>0.5 & MOI_per_cluster>=2]
# # list_core_T7BC2[["113"]] <- c()
# # list_core_T7BC2[["85"]] <- c()
# 
# 
# # # # # # # # # # # # # # # # # # # # # 
# # assess uniformity of cluster
# # # # # # # # # # # # # # # # # # # # # 
# library(ggplotify)
# 
# for (cluster_oi in seurat_clusters){
#   print(cluster_oi)
#   seur_obj_cluster_oi <- subset(BC_obj_sub, subset =( seurat_clusters==cluster_oi) )
#   BC_count_mat <- seur_obj_cluster_oi@assays$RNA$counts
#   if (length(list_core_T7BC2[[as.character(cluster_oi)]])>1){
#     BC_count_mat_core <- BC_count_mat[list_core_T7BC2[[as.character(cluster_oi)]],]
#     # plt_hclust <- heatmap(as.matrix(BC_count_mat_core))
#     
#     
#     # core_pBC_plot
#     fig_name <- sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling/heatmap_T7BC_count_across_assigned_cells_3umi/heatmap_T7BC_counts_cluster_%d_20240909.pdf",
#                         cluster_oi)
#     pdf(fig_name,width=5,height=5)
#     
#     
#     heatmap(as.matrix(t(BC_count_mat_core)),
#             #scale="none",
#             margins = c(10, 10),
#             cexRow = 0.3,
#             cexCol=0.3)
#     # print(plt_hclust)
#     dev.off()
#     
#   }
# }
# 
# 
# 
# 
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # Jaccard index similarity analysis to flag identical clones
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 
# clones_id <- names(list_core_T7BC2)
# 
# # clones_id <- paste0("seurat",names(list_core_T7BC2))
# # names(list_core_T7BC2) <- clones_id
# 
# df_pairs_clones <- CombPairs(clones_id)
# 
# for (id in seq(dim(df_pairs_clones)[1])){
#   
#   
#   if ((id %% 100)==0){
#     print(id)
#   }
#   
#   BC_clone1 <- list_core_T7BC2[[df_pairs_clones$X1[id]]] #df_all_clones %>% filter(clone_id==df_pairs_clones$X1[id]) %>% pull(real_barcode1)
#   BC_clone2 <- list_core_T7BC2[[df_pairs_clones$X2[id]]] #df_all_clones %>% filter(clone_id==df_pairs_clones$X2[id]) %>% pull(real_barcode1)
#   
#   jaccard_pair_BC <- length(intersect(BC_clone1,BC_clone2))/length(union(BC_clone1,BC_clone2))
# 
#   df_pairs_clones$jaccard_pair_BC[id] <- jaccard_pair_BC 
# 
# }
# 
# dev.set(4)
# ggplot(df_pairs_clones) + stat_bin(aes(x=jaccard_pair_BC),color="firebrick",geom="step")+
#   # stat_bin(aes(x=jaccard_pair_BC2),color="blue",geom="step")+
#   scale_y_log10()
# 
# df_pairs_clones %>% filter(jaccard_pair_BC>0.2)
# 
# df_pairs_clones %>% filter(jaccard_pair_BC>0.05)
# 
# 
# # spot check weird ones:
# list_core_T7BC2[["32"]]
# list_core_T7BC2[["38"]]
# intersect(list_core_T7BC2[["47"]],list_core_T7BC2[["48"]])
# 
# # Convert the data frame to a matrix
# 
# jac_mat_BC_v2 <- matrix(data=NA,nrow=length(clones_id),ncol=length(clones_id))
# rownames(jac_mat_BC_v2) <- clones_id
# colnames(jac_mat_BC_v2) <- clones_id
# for (id in seq(dim(df_pairs_clones)[1])){
#   jac_mat_BC_v2[df_pairs_clones$X1[id],df_pairs_clones$X2[id]] <- df_pairs_clones$jaccard_pair_BC[id]
# }
# 
# 
# dist_thresh <- 0.1
# conn_mat <- jac_mat_BC_v2>=dist_thresh
# g  <- graph_from_adjacency_matrix(conn_mat)
# clu <- igraph:::components(g)
# 
# graph <- graph_from_adjacency_matrix(conn_mat, mode = 'undirected')
# dev.set(4)
# plot(graph,vertex.size=3)
# 
# list_core_T7BC2[["27"]]
# list_core_T7BC2[["36"]]
# list_core_T7BC2[["40"]] 
# intersect(list_core_T7BC2[["27"]],list_core_T7BC2[["36"]])
# intersect(list_core_T7BC2[["27"]],list_core_T7BC2[["40"]])
# intersect(list_core_T7BC2[["27"]],list_core_T7BC2[["51"]])
# intersect(list_core_T7BC2[["40"]],list_core_T7BC2[["51"]])
# intersect(list_core_T7BC2[["36"]],list_core_T7BC2[["51"]])
# intersect(list_core_T7BC2[["19"]],list_core_T7BC2[["55"]])
# intersect(list_core_T7BC2[["12"]],list_core_T7BC2[["16"]])
# 
# 
# # identification of barcodes found in all clones?
# all_BCs <- list_core_T7BC2 %>% unlist() %>% unique()
# n_cluster_per_BC <- c()
# for (id in seq(length(all_BCs))){
#   BC_oi <- all_BCs[id]
#   
#   n_clonotypes_oi <- lapply(list_core_T7BC2, function(x) BC_oi %in% x) %>% unlist() %>% sum()
#   n_cluster_per_BC <- c(n_cluster_per_BC,n_clonotypes_oi)
# }
# 
# all_BCs[n_cluster_per_BC>1]
# 
# df_BC_per_cluster <- data.frame(BC=all_BCs,number_clonotypes=n_cluster_per_BC) %>% arrange(desc(number_clonotypes))
# 
# ggplot()+stat_bin(aes(x=n_cluster_per_BC))
# 
# 
# # assess bleedover hypothesis: group cells by cluster and get average UMI representation
# 
# df_cell_clusters <- BC_obj_sub@meta.data
# df_cell_clusters2 <- df_cell_clusters %>% transform(cellBC=rownames(df_cell_clusters))
# 
# df_T7BC_valid_hi2_w_clusters <- df_T7BC_valid_hi2 %>% left_join(df_cell_clusters2 %>% select(cBC_pdT=cellBC,seurat_clusters))
# 
# df_cell_clusters2 %>% select(cBC_pdT=cellBC,seurat_clusters) %>% group_by(seurat_clusters) %>% summarize(n_cells=length(cBC_pdT))
# 
# df_BC_per_cluster <- df_T7BC_valid_hi2_w_clusters %>% group_by(seurat_clusters,BC_pair) %>% 
#   summarize(n_cells=length(unique(cBC_pdT)),
#             mean_UMI_per_cell=mean(UMIs_BC)) %>% left_join(cells_per_cluster %>% select(seurat_clusters,cells_per_cluster))
# 
# df_BC_per_cluster2 <- df_BC_per_cluster %>% transform(fraction_cells_w_shuffleBC=n_cells/cells_per_cluster)
# 
# df_BC_pairs_sum_UMI %>% filter(BC_pair=="TAATACAGACACGCCTCCGG_CS1__TTTGTTGACCAGATCGCTAA_CS2")
# 
# df_BC_pairs_sum_UMI %>% filter(BC_pair=="CACTACGGCGCCCTAGTACT_CS1__ACGTTCGTCTGAGTCAAATG_CS2")
# df_BC_pairs_sum_UMI %>% filter(BC_pair=="GGTTACCCAGACCAATAAAA_CS1__GTATAGCATTAGTCATCAAC_CS2")
# 
# df_BC_per_cluster2 %>% filter(BC_pair=="TAATACAGACACGCCTCCGG_CS1__TTTGTTGACCAGATCGCTAA_CS2") %>% data.frame
# df_BC_per_cluster2 %>% filter(BC_pair=="GGTTACCCAGACCAATAAAA_CS1__GTATAGCATTAGTCATCAAC_CS2") %>% data.frame
# df_BC_per_cluster2 %>% filter(seurat_clusters=="0") %>% arrange(desc(mean_UMI_per_cell))
# 
# 
# clu$membership[clu$membership %in% which(clu$csize>1)] %>% table()
# 
# df_clone_membership_no_doub <- data.frame(clone_id=names(clu$membership),
#                                           conn_clu=clu$membership)
# 
# df_consensus_clone_name <- df_clone_membership_no_doub %>% group_by(conn_clu) %>% summarize(consensus_clone_id=paste0(clone_id,collapse="_"))
# 
# df_clone_membership_no_doub2 <- df_clone_membership_no_doub %>% left_join(df_consensus_clone_name)
# 
# 
# 
# # from heatmap analysis
# problematic_clones <- c("80","82","97")
# df_clone_membership_no_doub3 <- df_clone_membership_no_doub2 %>% filter(!(clone_id %in% problematic_clones))
# 
# consensus_clone_ids <- df_clone_membership_no_doub3 %>% pull(consensus_clone_id) %>% unique()
# df_consensus_clone_seurat <- data.frame()
# for (cclone_id in consensus_clone_ids){
#   
#   print(cclone_id)
#   df_oi <- df_clone_membership_no_doub3 %>% filter(consensus_clone_id==cclone_id)
#   
#   list_BCs_oi <- list()
#   list_BCs_oi <- list_core_T7BC2[df_oi$clone_id]
#   # bla
#   # for (clone_id in (df_oi %>% pull(clone_id))){
#   #   list_BCs_oi[clone_id] <- list_core_T7BC2[[clone_id]]
#   # }
#   print(lengths(list_BCs_oi))
#   
#   intersected_list_BC <- Reduce(intersect, list_BCs_oi)
#   # print(length(intersected_list_BC))
#   
#   df_consensus_clone_seurat <- rbind(df_consensus_clone_seurat,
#                                      data.frame(consensus_clone_id=cclone_id,
#                                                 BC_pair=intersected_list_BC))
# }
# 
# write.table(df_consensus_clone_seurat,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/denovo_clonotypes_round1_3umi_table_BC_pairs_20231220.txt",
#             quote=FALSE,row.names=FALSE,sep="\t")
# 
# 
# df_consensus_clone_seurat %>% group_by(consensus_clone_id) %>% summarize(n_BC=length(BC_pair)) %>% pull(n_BC) %>% quantile
# df_consensus_clone_seurat %>% group_by(consensus_clone_id) %>% summarize(n_BC=length(BC_pair)) %>% pull(n_BC) %>% mean
# 
# 
# 
# 
# 
# 
# # # printing the count table for the clusters outputs
# # date_str <- str_replace_all(Sys.Date(),"-","") 
# # dir_out_clone <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/test_de_novo_clones_count_tables"
# # BC_thresh <- 1
# # 
# # clone_oi <- 115
# # for (clone_oi in seq(1,130)){
# #   print(clone_oi)
# #   cells_clone <- colnames(BC_obj_sub)[which(BC_obj_sub$seurat_clusters==clone_oi)]
# #   
# #   # printing count output
# #   mat_oi <- T7BC_mat[cells_clone,]
# #   sum_BC <- colSums(mat_oi)
# #   
# #   idx_sort <- sort(sum_BC,index.return=TRUE,decreasing=TRUE)
# #   
# #   sorted_mat_oi <- t(mat_oi[,idx_sort$ix]) %>% as.matrix
# #   frac_cells_w_oBC <-  colSums(mat_oi[,idx_sort$ix]>BC_thresh)/length(cells_clone)
# #   
# #   mean_BC_count <-  colSums(mat_oi[,idx_sort$ix])/length(cells_clone)
# #   
# #   sorted_mat_oi_subset <- sorted_mat_oi[mean_BC_count>3/length(cells_clone),]
# #   
# #   n_cells_w_BC <-  colSums(sorted_mat_oi_subset>BC_thresh)
# #   idx_cells <- sort(n_cells_w_BC,index.return=TRUE,decreasing=TRUE)
# #   sorted_mat_oi_subset2 <- sorted_mat_oi_subset[,idx_cells$ix]
# #   
# #   write.table(sorted_mat_oi_subset2,
# #               sprintf("%s/clone_%s_%s.txt",dir_out_clone,clone_oi,date_str),
# #               sep="\t",quote=FALSE)
# # }
