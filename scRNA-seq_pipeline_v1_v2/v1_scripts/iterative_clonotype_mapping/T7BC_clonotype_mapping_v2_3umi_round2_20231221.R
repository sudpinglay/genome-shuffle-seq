
library(Seurat)
library(Matrix)

df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/valid_T7BC_Cre_v3_CS_20231220.txt",
                            header=TRUE)


# filter out the round1 assigned cells! 
df_assignment_round1 <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_v2/round1_assignment_cells_to_denovo3umi_bulk_clonotypes_20231221.txt",
                                   header=TRUE)

low_capture_cells <- df_assignment_round1 %>% filter(assignment_category %in% c("no_capture","low_capture")) %>% pull(cellBC)


BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split("_CS")
df_T7BC_valid$BC_pair <- BC_all %>% lapply("[[",1) %>% unlist()
df_T7BC_valid$capture_direction <- paste0("CS",BC_all %>% lapply("[[",2) %>% unlist())
df_T7BC_valid <- df_T7BC_valid %>% select(-BC)

df_T7BC_valid2 <- df_T7BC_valid %>% group_by(cBC_pdT,BC_pair) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC))
df_T7BC_valid2 <- df_T7BC_valid2 %>% filter(UMIs_BC>=3)
df_T7BC_valid2 <- df_T7BC_valid2 %>% filter(cBC_pdT %in% low_capture_cells)


df_BC_pairs_sum_UMI <- df_T7BC_valid2 %>% group_by(BC_pair) %>% summarize(sum_UMI=sum(UMIs_BC))
df_BC_UMI_per_cell <- df_T7BC_valid2 %>% group_by(cBC_pdT) %>% summarize(sum_T7_UMI_per_cell=sum(UMIs_BC))


thresh_T7umi_per_cell <- 10
thresh_T7_per_BC_pair <- 20
dev.set(4)
ggplot(df_BC_UMI_per_cell) + stat_bin(aes(x=sum_T7_UMI_per_cell),geom="step")+
  scale_x_log10()+scale_y_log10()+
  geom_vline(xintercept=thresh_T7umi_per_cell,color="red",linetype="dashed")

ggplot(df_BC_pairs_sum_UMI) + stat_bin(aes(x=sum_UMI),geom="step")+
  scale_x_log10()+scale_y_log10()+
  geom_vline(xintercept=thresh_T7_per_BC_pair,color="red",linetype="dashed")

hi_BC_pairs <- df_BC_pairs_sum_UMI %>% filter(sum_UMI>thresh_T7_per_BC_pair ) %>% pull(BC_pair)
high_cells <- df_BC_UMI_per_cell %>% filter(sum_T7_UMI_per_cell>thresh_T7umi_per_cell) %>% pull(cBC_pdT)
length(hi_BC_pairs)
length(high_cells)

df_T7BC_valid_hi <- df_T7BC_valid2 %>% filter(BC_pair %in% hi_BC_pairs) %>% 
  filter(cBC_pdT %in% high_cells)


# generate T7BC matrix on filtered cells and BC pairs: 
cBC_list <- df_T7BC_valid_hi %>% pull(cBC_pdT) %>% unique()
cBC_list <- data.frame(cBC_id=seq(1:length(cBC_list)),cBC_pdT=cBC_list)
T7BC_list <- df_T7BC_valid_hi %>% pull(BC_pair) %>% unique()
T7BC_list <- data.frame(T7BC_id=seq(1:length(T7BC_list)),BC_pair=T7BC_list)

df_T7BC_valid_hi2 <-  df_T7BC_valid_hi %>% left_join(cBC_list) %>% left_join(T7BC_list)
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
DimPlot(BC_obj_sub,
        reduction = "umap", label=TRUE,
        label.box = FALSE,label.size =1.5,
        pt.size =0.5)+theme(legend.position="none")#+labs(title=sample_cell_oi)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# studying count distribution in the T7 BC per clone
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

seurat_clusters <- 0:46
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
for (cluster_oi in seurat_clusters){
  
  print(cluster_oi)
  seur_obj_cluster_oi <- subset(BC_obj_sub, subset =( seurat_clusters==cluster_oi) )
  
  count_mat_oi <- seur_obj_cluster_oi@assays$RNA@counts
  
  # sum_UMI_BC <- sort(rowSums(count_mat_oi),decreasing = TRUE)
  
  # for each pBC, determine fraction of cBC with >10 UMI to that pBC
  n_cell_cluster <- dim(count_mat_oi)[2]
  frac_cBC <- (rowSums(count_mat_oi>UMI_thresh,na.rm=TRUE))/n_cell_cluster
  
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
  fig_name <- sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_round2/cluster_UMAP_w_distribution_T7BC_representation_per_cell_3umi/clusters_T7BCs_round2_id_%d_20231221.pdf",
                      cluster_oi)
  pdf(fig_name,width=10,height=5)
  print(Dim_plt+core_pBC_plot)
  dev.off()
  # 
  
}


# threshold some base quality on cluster
dev.set(4)
ggplot()+stat_bin(aes(x=max_frac_per_cluster),bins=20)+
  geom_vline(xintercept=0.5,color="red")

MOI_per_cluster <- lengths(list_core_T7BC)

ggplot()+stat_bin(aes(x=MOI_per_cluster),bins=20)+
  geom_vline(xintercept=7,color="red")+scale_x_log10()


list_core_T7BC2 <- list_core_T7BC[max_frac_per_cluster>0.5 & MOI_per_cluster>=7]
# list_core_T7BC2[["113"]] <- c()
# list_core_T7BC2[["85"]] <- c()


# # # # # # # # # # # # # # # # # # # # 
# assess uniformity of cluster
# # # # # # # # # # # # # # # # # # # # 
library(ggplotify)

for (cluster_oi in seurat_clusters){
  print(cluster_oi)
  seur_obj_cluster_oi <- subset(BC_obj_sub, subset =( seurat_clusters==cluster_oi) )
  BC_count_mat <- seur_obj_cluster_oi@assays$RNA@counts
  if (length(list_core_T7BC2[[as.character(cluster_oi)]])>1){
    BC_count_mat_core <- BC_count_mat[list_core_T7BC2[[as.character(cluster_oi)]],]
    # plt_hclust <- heatmap(as.matrix(BC_count_mat_core))
    
    
    # core_pBC_plot
    fig_name <- sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_round2/heatmap_T7BC_count_across_assigned_cells_3umi/heatmap_T7BC_counts_round2_cluster_%d_20231221.pdf",
                        cluster_oi)
    pdf(fig_name,width=5,height=5)
    
    
    heatmap(as.matrix(t(BC_count_mat_core)),
            #scale="none",
            margins = c(10, 10),
            cexRow = 0.3,
            cexCol=0.3)
    # print(plt_hclust)
    dev.off()
    
  }
}





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Jaccard index similarity analysis to flag identical clones
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

clones_id <- names(list_core_T7BC2)

# clones_id <- paste0("seurat",names(list_core_T7BC2))
# names(list_core_T7BC2) <- clones_id

df_pairs_clones <- CombPairs(clones_id)

for (id in seq(dim(df_pairs_clones)[1])){
  
  
  if ((id %% 100)==0){
    print(id)
  }
  
  BC_clone1 <- list_core_T7BC2[[df_pairs_clones$X1[id]]] #df_all_clones %>% filter(clone_id==df_pairs_clones$X1[id]) %>% pull(real_barcode1)
  BC_clone2 <- list_core_T7BC2[[df_pairs_clones$X2[id]]] #df_all_clones %>% filter(clone_id==df_pairs_clones$X2[id]) %>% pull(real_barcode1)
  
  jaccard_pair_BC <- length(intersect(BC_clone1,BC_clone2))/length(union(BC_clone1,BC_clone2))

  df_pairs_clones$jaccard_pair_BC[id] <- jaccard_pair_BC 

}

dev.set(4)
ggplot(df_pairs_clones) + stat_bin(aes(x=jaccard_pair_BC),color="firebrick",geom="step")+
  # stat_bin(aes(x=jaccard_pair_BC2),color="blue",geom="step")+
  scale_y_log10()

df_pairs_clones %>% filter(jaccard_pair_BC>0.05)

# # spot check weird ones:
# list_core_T7BC2[["31"]]
# list_core_T7BC2[["90"]]
# intersect(list_core_T7BC2[["31"]],list_core_T7BC2[["90"]])

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
plot(graph,vertex.size=3)



clu$membership[clu$membership %in% which(clu$csize>1)] %>% table()






df_clone_membership_no_doub <- data.frame(clone_id=names(clu$membership),
                                          conn_clu=clu$membership)

df_consensus_clone_name <- df_clone_membership_no_doub %>% group_by(conn_clu) %>% summarize(consensus_clone_id=paste0(clone_id,collapse="_"))

df_clone_membership_no_doub2 <- df_clone_membership_no_doub %>% left_join(df_consensus_clone_name)

# from heatmap analysis
problematic_clones <- c("13","31","34")
df_clone_membership_no_doub3 <- df_clone_membership_no_doub2 %>% filter(!(clone_id %in% problematic_clones))

consensus_clone_ids <- df_clone_membership_no_doub3 %>% pull(consensus_clone_id) %>% unique()
df_consensus_clone_seurat <- data.frame()
for (cclone_id in consensus_clone_ids){
  
  print(cclone_id)
  df_oi <- df_clone_membership_no_doub3 %>% filter(consensus_clone_id==cclone_id)
  
  list_BCs_oi <- list()
  list_BCs_oi <- list_core_T7BC2[df_oi$clone_id]
  # bla
  # for (clone_id in (df_oi %>% pull(clone_id))){
  #   list_BCs_oi[clone_id] <- list_core_T7BC2[[clone_id]]
  # }
  print(lengths(list_BCs_oi))
  
  intersected_list_BC <- Reduce(intersect, list_BCs_oi)
  # print(length(intersected_list_BC))
  
  df_consensus_clone_seurat <- rbind(df_consensus_clone_seurat,
                                     data.frame(consensus_clone_id=cclone_id,
                                                BC_pair=intersected_list_BC))
}

write.table(df_consensus_clone_seurat,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_round2/denovo_clonotypes_round2_3umi_table_BC_pairs_20231221.txt",
            quote=FALSE,row.names=FALSE,sep="\t")


df_consensus_clone_seurat %>% group_by(consensus_clone_id) %>% summarize(n_BC=length(BC_pair)) %>% pull(n_BC) %>% quantile
df_consensus_clone_seurat %>% group_by(consensus_clone_id) %>% summarize(n_BC=length(BC_pair)) %>% pull(n_BC) %>% mean






# # printing the count table for the clusters outputs
# date_str <- str_replace_all(Sys.Date(),"-","") 
# dir_out_clone <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/test_de_novo_clones_count_tables"
# BC_thresh <- 1
# 
# clone_oi <- 115
# for (clone_oi in seq(1,130)){
#   print(clone_oi)
#   cells_clone <- colnames(BC_obj_sub)[which(BC_obj_sub$seurat_clusters==clone_oi)]
#   
#   # printing count output
#   mat_oi <- T7BC_mat[cells_clone,]
#   sum_BC <- colSums(mat_oi)
#   
#   idx_sort <- sort(sum_BC,index.return=TRUE,decreasing=TRUE)
#   
#   sorted_mat_oi <- t(mat_oi[,idx_sort$ix]) %>% as.matrix
#   frac_cells_w_oBC <-  colSums(mat_oi[,idx_sort$ix]>BC_thresh)/length(cells_clone)
#   
#   mean_BC_count <-  colSums(mat_oi[,idx_sort$ix])/length(cells_clone)
#   
#   sorted_mat_oi_subset <- sorted_mat_oi[mean_BC_count>3/length(cells_clone),]
#   
#   n_cells_w_BC <-  colSums(sorted_mat_oi_subset>BC_thresh)
#   idx_cells <- sort(n_cells_w_BC,index.return=TRUE,decreasing=TRUE)
#   sorted_mat_oi_subset2 <- sorted_mat_oi_subset[,idx_cells$ix]
#   
#   write.table(sorted_mat_oi_subset2,
#               sprintf("%s/clone_%s_%s.txt",dir_out_clone,clone_oi,date_str),
#               sep="\t",quote=FALSE)
# }
