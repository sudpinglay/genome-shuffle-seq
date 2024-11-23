
library(tidyverse)
library(stringr)
library(Matrix)
library(Biostrings)


# read T7BC data

sample_oi <- "Cre"
df_T7BC_CS <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/processed_out/Cre_CS_get_R2_2BC_umi_cleaned_20231206.txt.gz",
                      header=TRUE)

sample_oi <- "Parental"
df_T7BC_CS <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/processed_out/Par_CS_get_R2_2BC_umi_cleaned_20231206.txt.gz",
                         header=TRUE)

df_BC2 <- df_T7BC_CS %>% select(cBC,BC=mBC,reads_BC=n_reads_filtered,UMIs_BC=filtered_corrected_UMIs) #filtered_corrected_UMIs)
# df_BC3 <- df_BC2 %>% transform(CS_end=ifelse(substr(BC,1,4)=="AAGC","CS1","CS2"),BC2=substr(BC,5,24)) %>% select(-BC)

# convert the cellBC to the right one (from the same bead, CS1/CS2 cellBC are different from poly-dT): 
cBC_CS <- df_BC2$cBC
start_cBC <- substr(cBC_CS,1,7)
end_cBC <- substr(cBC_CS,10,18)
complemented_center_cBC <- as.character(complement(DNAStringSet(substr(cBC_CS,8,9))))
cBC_pdT <- paste0(start_cBC,complemented_center_cBC,end_cBC)
df_BC4 <- df_BC2 %>% rename(c("cBC"="cBC_CS")) 
df_BC4$cBC_pdT <- cBC_pdT

df_T7BC_CS2 <- df_BC4


# get set of "bona fide cells" from the GEx criteria 
file_valid_cellBC <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/gex_prepro_v2/metadata_GEx_cellBC_SP_T7loxShuffle_combined_20231220.txt"
df_valid_cellBC_all <- read.table(file_valid_cellBC,sep="\t",header=TRUE) #%>% transform(raw_cellBC_pdT=paste0(raw_cellBC_pdT,))

df_valid_cellBC_all <-  df_valid_cellBC_all %>% 
  transform(raw_cellBC_pdT=str_split(raw_cellBC_pdT,"_") %>% lapply("[[",2) %>% unlist())

# restrict to sample of interest
df_valid_cellBC_all2 <- df_valid_cellBC_all %>% filter(orig.ident==sample_oi)

df_T7BC_CS3 <- df_T7BC_CS2 %>% transform(valid_cBC=(cBC_pdT %in% df_valid_cellBC_all2$raw_cellBC_pdT))

# loose filter on counts
df_T7BC_CS_valid <- df_T7BC_CS3 %>% filter(valid_cBC) %>% select(-valid_cBC) %>%
  filter(UMIs_BC>=2 | reads_BC/UMIs_BC>=8)


write.table(df_T7BC_CS_valid,
            sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/valid_T7BC_%s_v3_CS_20231220.txt",sample_oi),
            sep="\t",row.names=FALSE,quote=FALSE)

dev.set(4)
ggplot(df_T7BC_CS_valid) + geom_bin_2d(aes(x=reads_BC,y=UMIs_BC))+
  scale_x_log10()+scale_y_log10()+coord_fixed()+labs(title=sample_oi)






# df_T7BC_pdT <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/processed_out/Cre_dT_get_R2_2BC_umi_cleaned_20231206.txt.gz",
#                          header=TRUE)
# 
# df_T7BC_pdT2 <- df_T7BC_pdT %>% select(cBC_pdT=cBC,BC=mBC,reads_BC=n_reads_filtered,UMIs_BC=filtered_corrected_UMIs) #filtered_corrected_UMIs)


# df_T7BC_pdT3 <- df_T7BC_pdT2 %>% transform(valid_cBC=(cBC_pdT %in% df_valid_cellBC_all2$raw_cellBC_pdT))


# df_T7BC_both_capture <- rbind(df_T7BC_CS3 %>% select(-valid_cBC,-cBC_CS) %>% transform(capture="CS"),
#                               df_T7BC_pdT3 %>% select(-valid_cBC) %>% transform(capture="pdT"))



# df_T7BC_both_capture_wide <- df_T7BC_both_capture %>% pivot_wider(id_cols=c(BC,cBC_pdT),
# values_from=c(reads_BC,UMIs_BC),
# names_from=capture)





df_T7BC_pdT_valid <- df_T7BC_pdT3 %>% filter(valid_cBC) %>% select(-valid_cBC) %>%
  filter(UMIs_BC>=2 | reads_BC/UMIs_BC>=8)

df_T7BC_pdT3 %>% filter(valid_cBC) %>% select(-valid_cBC) %>%
  filter(UMIs_BC>=2) %>% dim()

write.table(df_T7BC_pdT_valid,
            "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/valid_T7BC_Cre_v2_pdT_20231210.txt",
            sep="\t",row.names=FALSE,quote=FALSE)



df_T7BC_both_capture_wide[is.na(df_T7BC_both_capture_wide)] <- 0

p_count <- 0.3
ggplot(df_T7BC_both_capture_wide) + stat_bin_2d(aes(x=UMIs_BC_CS+p_count,y=UMIs_BC_pdT+p_count))+
  scale_x_log10()+scale_y_log10()+coord_fixed()+geom_abline(linetype="dotted",color="orange")+
  geom_hline(yintercept=2,color="red",linewidth=0.2)+
  geom_vline(xintercept=2,color="red",linewidth=0.2)



write.table(df_T7BC_both_capture_wide %>% arrange(cBC_pdT),"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/T7_BC_both_capture_wide_20231210.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


df_T7BC_both_capture_wide %>% filter(UMIs_BC_CS>1 | UMIs_BC_pdT>1) %>% transform(CS_better_pdT=UMIs_BC_CS>=UMIs_BC_pdT) %>% pull(CS_better_pdT) %>% mean
df_T7BC_both_capture_wide %>% filter(UMIs_BC_CS>1 | UMIs_BC_pdT>1) %>% transform(CS_better_pdT=UMIs_BC_CS<=UMIs_BC_pdT) %>% pull(CS_better_pdT) %>% mean


df_sum_T7_BC_per_cell <- df_T7BC_CS3 %>% transform(noise_BC = (reads_BC==1 & UMIs_BC==1)) %>% group_by(cBC_pdT,valid_cBC,noise_BC) %>%
  summarize(total_T7_UMI=sum(UMIs_BC))

df_sum_T7_BC_per_cell <- df_T7BC_pdT3 %>% transform(noise_BC = (reads_BC==1 & UMIs_BC==1)) %>% group_by(cBC_pdT,valid_cBC,noise_BC) %>%
  summarize(total_T7_UMI=sum(UMIs_BC))

plt_logy <- ggplot(df_sum_T7_BC_per_cell) + stat_bin(aes(x=total_T7_UMI,color=valid_cBC,linetype=noise_BC),geom="step",position="identity")+
  scale_x_log10()+scale_y_log10()
plt_liny <- ggplot(df_sum_T7_BC_per_cell) + stat_bin(aes(x=total_T7_UMI,color=valid_cBC,linetype=noise_BC),geom="step",position="identity")+
  scale_x_log10()
plt_liny+plt_logy


quantile(df_sum_T7_BC_per_cell %>% filter(valid_cBC & !noise_BC) %>% pull(total_T7_UMI))


library(patchwork)
plt_pdT <- ggplot(df_T7BC_pdT3)+stat_bin(aes(x=UMIs_BC,color=valid_cBC),geom="step", position="identity")+
  scale_x_log10()+scale_y_log10()

plt_CS <- ggplot(df_T7BC_CS3)+stat_bin(aes(x=UMIs_BC,color=valid_cBC),geom="step", position="identity")+
  scale_x_log10()+scale_y_log10()

plt_pdT+plt_CS

plt_pdT <- ggplot(df_T7BC_pdT3)+stat_bin_2d(aes(x=reads_BC,y=UMIs_BC))+
  facet_wrap(~valid_cBC)+
  scale_x_log10(limits=c(0.5,1E4))+scale_y_log10(limits=c(0.5,200))

plt_CS <- ggplot(df_T7BC_CS3)+stat_bin_2d(aes(x=reads_BC,y=UMIs_BC))+
  facet_wrap(~valid_cBC)+
  scale_x_log10(limits=c(0.5,1E4))+scale_y_log10(limits=c(0.5,200))
plt_pdT+plt_CS


plt_pdT <- ggplot(df_T7BC_pdT3)+stat_bin_2d(aes(x=UMIs_BC,y=reads_BC/UMIs_BC-1))+
  facet_wrap(~valid_cBC)+
  scale_x_log10(limits=c(0.5,200))+scale_y_log10(limits=c(0.5,5E3))

plt_CS <- ggplot(df_T7BC_CS3)+stat_bin_2d(aes(x=UMIs_BC,y=reads_BC/UMIs_BC-1))+
  facet_wrap(~valid_cBC)+
  scale_x_log10(limits=c(0.5,200))+scale_y_log10(limits=c(0.5,5E3))
plt_pdT+plt_CS


# df_GEx_raw <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq055SP_T7_shuffle_10X/T7_BC_prepro/raw_GEx_data_sample1_Cre_20231121.txt",
#                          header=TRUE)
# 
# UMI_thresh <- 750
# mt_thresh <- c(1,12.5)
# dev.new()
# ggplot(df_GEx_raw)+
#   stat_bin2d(aes(x=nCount_RNA,y=percent_mt,fill=log(..count..),color=log(..count..)),bins=100)+
#   geom_point(data=df_GEx_raw %>% filter(cBC_pdT %in% df_valid_cellBC$cBC_pdT),
#              aes(x=nCount_RNA,y=percent_mt),color="red",size=0.1)+
#   scale_x_log10(limits=c(30,20000))+
#   scale_y_log10(limits=c(0.3,100))+
#   annotation_logticks()+
#   geom_vline(xintercept=UMI_thresh, color="red", linetype="dashed")+
#   geom_hline(yintercept=mt_thresh, color="red", linetype="dashed")





















df_valid_cellBC <-  df_valid_cellBC %>% 
  transform(raw_cellBC_pdT=str_split(raw_cellBC_pdT,"_") %>% lapply("[[",3) %>% unlist())

valid_cBC <- df_T7BC4$cBC_pdT %in% df_valid_cellBC$raw_cellBC_pdT

# threshold on valid cellBC and high count T7BC (might be worth revisiting and being more stringent)
df_T7BC_valid <- df_BC4 %>% filter(cBC_pdT %in% df_valid_cellBC$raw_cellBC_pdT) %>%
  filter(UMIs_BC>=2 | reads_BC/UMIs_BC>=5)

write.table(df_T7BC_valid,
            "/Users/jbl/Documents/UW/Data/sequencing_runs/seq055SP_T7_shuffle_10X/valid_cellBC_T7BC_sample1_Cre_ShuffleLox_20231116.txt",
            sep="\t",row.names=FALSE,quote=FALSE)


df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq055SP_T7_shuffle_10X/valid_cellBC_T7BC_sample1_Cre_ShuffleLox_20231116.txt.gz",
                            header=TRUE)

df_T7BC_sum <- df_T7BC_valid %>% group_by(BC2) %>% summarize(sum_UMI=sum(UMIs_BC))


df_T7BC_hi <- df_T7BC_sum %>% filter(sum_UMI>10)

dev.set(5)
ggplot(df_T7BC_sum) + stat_bin(aes(x=sum_UMI),geom="step")+
  scale_x_log10()+scale_y_log10()


# df_cellBC <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq055SP_T7_shuffle_10X/GEx_prepro/metadata_GEx_cellBC_SP_T7loxShuffle_20231121.txt",
#                         header=TRUE)





permute_id <- 1 #as.numeric(str_content[[1]][3])
# rep_oi <- str_content[[1]][1]
# cell_type_oi <- as.numeric(str_content[[1]][2])


# thresh_T7BC <- 11
# rep_oi <- "repA"
# cell_type_oi <- 0

# valid T7BC (list of T7BC with at least one cell with above threshold): for this don't need all counts. 
# df_T7BC_valid <- df_T7BC %>% filter( (biol_rep %in% c(rep_oi)) & UMIs_T7BC>=thresh_T7BC & p25_valid_T7BC & CRE_class=="devCRE") # & gex_cluster_id==cell_type_oi)

all_cells <- df_T7BC_valid %>% pull(cBC_pdT) %>% unique()

cBC_list <- df_T7BC_valid %>% pull(cBC_pdT) %>% unique()
cBC_list <- data.frame(cBC_id=seq(1:length(cBC_list)),cBC_pdT=cBC_list)
T7BC_list <- df_T7BC_valid %>% pull(BC2) %>% unique()
T7BC_list <- data.frame(T7BC_id=seq(1:length(T7BC_list)),BC2=T7BC_list)

df_T7BC_valid2 <-  df_T7BC_valid %>% left_join(cBC_list) %>% left_join(T7BC_list)

T7BC_mat <- sparseMatrix(j=df_T7BC_valid2$T7BC_id,
                         i=df_T7BC_valid2$cBC_id,
                         x=df_T7BC_valid2$UMIs_BC)

rownames(T7BC_mat) <- cBC_list$cBC_pdT
colnames(T7BC_mat) <- T7BC_list$BC2

library(Seurat)
BC_obj <- CreateSeuratObject(counts = t(T7BC_mat), min.cells = 0, min.features = 0)



# usual dimensional reduction
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


# printing the outputs
date_str <- str_replace_all(Sys.Date(),"-","") 
dir_out_clone <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq055SP_T7_shuffle_10X/test_clone_count_table"
BC_thresh <- 1


clone_oi <- 115

for (clone_oi in seq(1,122)){
  print(clone_oi)
  cells_clone <- colnames(BC_obj_sub)[which(BC_obj_sub$seurat_clusters==clone_oi)]
  
  # printing count output
  mat_oi <- T7BC_mat[cells_clone,]
  sum_BC <- colSums(mat_oi)
  
  idx_sort <- sort(sum_BC,index.return=TRUE,decreasing=TRUE)
  
  sorted_mat_oi <- t(mat_oi[,idx_sort$ix]) %>% as.matrix
  frac_cells_w_oBC <-  colSums(mat_oi[,idx_sort$ix]>BC_thresh)/length(cells_clone)
  
  mean_BC_count <-  colSums(mat_oi[,idx_sort$ix])/length(cells_clone)
  
  sorted_mat_oi_subset <- sorted_mat_oi[mean_BC_count>3/length(cells_clone),]
  
  n_cells_w_BC <-  colSums(sorted_mat_oi_subset>BC_thresh)
  idx_cells <- sort(n_cells_w_BC,index.return=TRUE,decreasing=TRUE)
  sorted_mat_oi_subset2 <- sorted_mat_oi_subset[,idx_cells$ix]
  
  write.table(sorted_mat_oi_subset2,
              sprintf("%s/clone_%s_%s.txt",dir_out_clone,clone_oi,date_str),
              sep="\t",quote=FALSE)
}




heatmap(sorted_mat_oi_subset2,Rowv = NULL,Colv=NULL)











# setwd("/Users/jbl/Documents/UW/projects/clonotype_mapping_v2_20220622")
# 
# # for cluster usage
# T7BC_mat <- readRDS("/net/shendure/vol10/projects/JBL/clonotype_mapping/nobackup/mEB_rep2B_20220704/T7BC_mat_mEB_rep2B_T7BC_thresh_10_20220704.RDS")


# threshold for assignment
thresh_T7BC <- 1
T7BC_mat_binary <- T7BC_mat>=thresh_T7BC

clone_list <- list()
n_T7BC <- dim(T7BC_mat)[2]
n_cells <- dim(T7BC_mat)[1]
p_val_thresh <- 0.05/(n_cells*n_cells/2)
cell_category <- vector(mode="character",length=n_cells)
  
shuffled_cell_ids <- sample(rownames(T7BC_mat),size=length(rownames(T7BC_mat)),replace=FALSE)
names(cell_category) <- shuffled_cell_ids
  
counter <- 1

# n_in_both <- 2
# n_in_clone_not_in_cell <- 0
# n_in_cell_not_in_clone <- 0
# n_in_either <- n_in_both
# fisher.test(rbind(c(n_in_both,n_in_clone_not_in_cell),c(n_in_cell_not_in_clone,n_T7BC-n_in_either)), alternative="greater")

for (cell_id in shuffled_cell_ids){
  
  if (length(clone_list)==0){
    clone_list[cell_id] <- cell_id
    
    # clone_list <- append(clone_list,cell_id)
    cell_category[cell_id] <- "new_clone"
    
  } else {
    
    
    T7BC_list_cell_oi <- which(T7BC_mat_binary[cell_id,])
    n_T7BC_cell_oi <- length(T7BC_list_cell_oi)
    
    member_clones <- c()
    for (clone_oi in names(clone_list)){
      
      T7BC_list_clone <- which(T7BC_mat_binary[clone_oi,])
      
      n_T7BC_clone_oi <- length(T7BC_list_clone)
      
      n_in_both <- length(intersect(T7BC_list_cell_oi,T7BC_list_clone))
      n_in_either <- length(union(T7BC_list_cell_oi,T7BC_list_clone))
      
      # build contingency table
      n_in_clone_not_in_cell <- n_T7BC_clone_oi-n_in_both
      n_in_cell_not_in_clone <- n_T7BC_cell_oi-n_in_both
      
      # one-sided Fisher's exact test
      p_val <- fisher.test(rbind(c(n_in_both,n_in_clone_not_in_cell),c(n_in_cell_not_in_clone,n_T7BC-n_in_either)), 
                           alternative="greater")$p.value
      
      if (p_val<p_val_thresh){
        member_clones <- c(member_clones,clone_oi)
      }
    }
    
    
    if (length(member_clones)==0){
      # no overlap found
      clone_list[cell_id] <- cell_id
      
      cell_category[cell_id] <- "new_clone"
      print(sprintf('%d, cell_id: %s, category: %s',counter, cell_id,"new_clone"))
      
    } else if (length(member_clones)==1){
      # one overlap found: plausible singlet
      clone_list[[member_clones]] <- c(clone_list[[member_clones]],cell_id)
      cell_category[cell_id] <- "singlet"
      
      print(sprintf('%d, cell_id: %s, category: %s',counter,cell_id,"singlet"))

    } else if (length(member_clones)>1){
      cell_category[cell_id] <- "doublet"
      print(sprintf('%d, cell_id: %s, category: %s',counter,cell_id,"doublet"))
      
    }
  }
  counter <- counter+1
}



saveRDS(cell_category,clone_list)

cell_category[cell_category=="doublet"] %>% head()

quantile(lengths(clone_list))
mean(lengths(clone_list)>1)

table(cell_category)



# singling out specific cells

cell_id_oi <- "ATCGGCGGTTATAGAG-1"


cell_id_oi <- "GTGCTTCCACTCCGAG-1"
clone_oi <- "GGGCGTTCAGAGACTG-1"

cell_category[cell_id_oi]

T7BC_list_cell_oi <- which(T7BC_mat_binary[cell_id_oi,])
n_T7BC_cell_oi <- length(T7BC_list_cell_oi)


member_clones <- c()
for (clone_oi in names(clone_list)){
  
  T7BC_list_clone <- which(T7BC_mat_binary[clone_oi,])
  
  n_T7BC_clone_oi <- length(T7BC_list_clone)
  n_T7BC_clone_oi
  
  n_in_both <- length(intersect(T7BC_list_cell_oi,T7BC_list_clone))
  n_in_both
  n_in_either <- length(union(T7BC_list_cell_oi,T7BC_list_clone))
  n_in_either
  
  # build contingency table
  n_in_clone_not_in_cell <- n_T7BC_clone_oi-n_in_both
  n_in_cell_not_in_clone <- n_T7BC_cell_oi-n_in_both
  
  # one-sided Fisher's exact test
  rbind(c(n_in_both,n_in_clone_not_in_cell),c(n_in_cell_not_in_clone,n_T7BC-n_in_either))
  p_val <- fisher.test(rbind(c(n_in_both,n_in_clone_not_in_cell),c(n_in_cell_not_in_clone,n_T7BC-n_in_either)), 
                       alternative="greater")$p.value
  p_val
  
  if (p_val<p_val_thresh){
    member_clones <- c(member_clones,clone_oi)
  }
}


clone_oi <- member_clones[1]
clone_oi <- member_clones[2]
clone_oi <- member_clones[3]
clone_oi <- member_clones[4]


clone_list[lengths(clone_list)>10] %>% tail

cell_list <- c("GGGCGTTCAGAGACTG-1",
               "TACAACGCAATTTCCT-1",
               "CACCGTTCACGATTCA-1",
               "CCTAACCAGAAGCTGC-1",
               "GTGCTTCCACTCCGAG-1")

cell_list <- clone_list[["AAGTGAATCGGTTCAA-1"]]

# library(DescTools)
df_pairs <- CombPairs(cell_list)


general_union_T7_BC <- c()

for (idx in seq(dim(df_pairs)[1])){
  
  cell1 <- df_pairs$X1[idx]
  cell2 <- df_pairs$X2[idx]
  
  print(cell1)
  print(cell2)
  
  T7BC_list_cell_oi <- which(T7BC_mat_binary[cell1,])
  n_T7BC_cell_oi <- length(T7BC_list_cell_oi)
  
  general_union_T7_BC <- c(general_union_T7_BC,T7BC_list_cell_oi)
  
  T7BC_list_clone <- which(T7BC_mat_binary[cell2,])
  
  general_union_T7_BC <- c(general_union_T7_BC,T7BC_list_clone)
  
  n_T7BC_clone_oi <- length(T7BC_list_clone)
  n_T7BC_clone_oi
  
  n_in_both <- length(intersect(T7BC_list_cell_oi,T7BC_list_clone))
  n_in_both
  n_in_either <- length(union(T7BC_list_cell_oi,T7BC_list_clone))
  n_in_either
  
  # build contingency table
  n_in_clone_not_in_cell <- n_T7BC_clone_oi-n_in_both
  n_in_cell_not_in_clone <- n_T7BC_cell_oi-n_in_both
  
  # one-sided Fisher's exact test
  print(rbind(c(n_in_both,n_in_clone_not_in_cell),c(n_in_cell_not_in_clone,n_T7BC-n_in_either)))
  p_val <- fisher.test(rbind(c(n_in_both,n_in_clone_not_in_cell),c(n_in_cell_not_in_clone,n_T7BC-n_in_either)), 
                       alternative="greater")$p.value
  print(p_val)
  
  df_pairs$n_cell1[idx] <- n_T7BC_cell_oi
  df_pairs$n_cell2[idx] <- n_T7BC_clone_oi
  df_pairs$n_both[idx] <- n_in_both
  df_pairs$p_val[idx] <- p_val
  
  # if (p_val<p_val_thresh){
  #   member_clones <- c(member_clones,clone_oi)
  # }
}
df_pairs
unique(general_union_T7_BC)

T7BC_list_clone1 <- which(T7BC_mat_binary[member_clones[1],])







table(cell_category)

# save output:
save("clone_list","cell_category","shuffled_cell_ids",
        file=sprintf("clone_list_Fisher_mEBs_rep_%s_permute_id%d_20220823.RData",rep_oi,permute_id))

