library(Seurat)
library(dplyr)
library(patchwork)
library(sctransform)
library(ggplot2)
library(tidyverse)
library(scales)
library(tictoc)
library(irlba)
# library(schex) # might need to install

library(DropletUtils)


library(pryr)


wd_oi <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_GEx"
setwd(wd_oi)

date_str <- "20240906"

# read raw data
reps <- c("Lane1","Lane2")



# # # PRE FILTERING
for (rep_oi in reps[2]){
  
  print(rep_oi)
  path_10x_mtx <- sprintf("%s/%s",wd_oi,rep_oi)
  gex_mat <- Read10X(data.dir =path_10x_mtx)
  
  
  # create separate Seurat objects
  gex <- CreateSeuratObject(counts = gex_mat, project = rep_oi, min.cells = 3, min.features = 50)
  gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^MT-")
  
  mt_thresh <- c(2,8)
  UMI_thresh <- 1500
  fig_QC <- ggplot()+stat_bin2d(aes(x=gex$nCount_RNA,y=gex$percent_mt,fill=log(..count..),color=log(..count..)),bins=100)+ #
    scale_x_log10(limits=c(30,50000))+
    scale_y_log10(limits=c(0.05,100))+
    annotation_logticks()+
    geom_vline(xintercept=UMI_thresh, color="red", linetype="dashed")+
    geom_hline(yintercept=mt_thresh, color="red", linetype="dashed")
  
  # fig_QC
  
  fig_QC_name <- sprintf("mt_vs_RNA_UMI_gex_%s_%s.pdf",rep_oi,date_str)
  
  pdf(fig_QC_name,width=4,height=4)
  print(fig_QC)
  dev.off()
  
  
  # subsetting on cells
  gex_2 <- subset(gex, subset = (percent_mt < mt_thresh[2]) & (percent_mt > mt_thresh[1]) & (nCount_RNA >  UMI_thresh ) )
  
  print(dim(gex))
  print(dim(gex_2))
  
  write10xCounts(sprintf("%s/%s_filtered_10x_mtx_%s/",wd_oi,rep_oi,date_str),
                 gex_2@assays$RNA$counts,
                 barcodes = colnames(gex_2),
                 gene.id = rownames(gex_2))
  
}

#gex_2 <- subset(gex, subset = (percent_mt < mt_thresh[2]) & (nCount_RNA >  UMI_thresh ) )



# 
# rep_oi <- "sample1_Cre"
# 
# print(rep_oi)
# path_10x_mtx <- sprintf("%s/%s",wd_oi,rep_oi)
# gex_mat <- Read10X(data.dir =path_10x_mtx)
# 
# 
# # create separate Seurat objects
# gex <- CreateSeuratObject(counts = gex_mat, project = rep_oi, min.cells = 3, min.features = 50)
# gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^mt-")
# 
# write.table(gex@meta.data %>% transform(cBC_pdT=row.names(gex@meta.data)),"raw_GEx_data_sample1_Cre_20231121.txt",
#             row.names=FALSE, quote=FALSE, sep="\t")
# 
# mt_thresh <- c(1,12)
# UMI_thresh <- 750
# fig_QC <- ggplot()+stat_bin2d(aes(x=gex$nCount_RNA,y=gex$percent_mt,fill=log(..count..),color=log(..count..)),bins=100)+
#   scale_x_log10(limits=c(30,50000))+
#   scale_y_log10(limits=c(0.05,100))+
#   annotation_logticks()+
#   geom_vline(xintercept=UMI_thresh, color="red", linetype="dashed")+
#   geom_hline(yintercept=mt_thresh, color="red", linetype="dashed")









# read raw data
reps <- c("Lane1","Lane2")
date_str <- "20240906"
# # # READ ALL W/ DOUBLET SCORES
for (rep_oi in reps){
  
  print(rep_oi)
  path_10x_mtx <- sprintf("%s/%s_filtered_10x_mtx_%s",wd_oi,rep_oi,date_str)
  gex_mat <- Read10X(data.dir =path_10x_mtx)
  
  # create separate Seurat objects
  gex <- CreateSeuratObject(counts = gex_mat, project = rep_oi, min.cells = 3, min.features = 50)
  gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^MT-")
  

  file_doublet <- sprintf("%s/%s_filtered_10x_mtx_%s/scrublet_score_scRNA%s_filtered_10x_mtx_%s.txt",
                          wd_oi,rep_oi,date_str, rep_oi,date_str)
  df_doublet <- read.table(file_doublet,sep="\t",header=TRUE)
  
  rownames(df_doublet) <- df_doublet$cell_bc
  gex$doublet_score <- df_doublet[Cells(gex),2]
  
  assign(sprintf("gex_%s",rep_oi),gex)

}


# # # # threshold on scrublet score
scrublet_thresh <- 0.4
gex_sample1_2 <- subset(gex_Lane1, subset = doublet_score<scrublet_thresh) 
gex_sample2_2 <- subset(gex_Lane2, subset = doublet_score<scrublet_thresh) 



# # # MERGE ALL IN ONE SEURAT OBJECT

Idents(gex_sample1_2) <- "Lane1"
Idents(gex_sample2_2) <- "Lane2"

saveRDS(gex_sample2_2,"seurat_obj_GEx_filtered_SP_shuffleBC_revision_Lane2_raw_20240906.RDS")

gex_all <- merge(gex_sample1_2, y = c(gex_sample2_2),
                       add.cell.ids = reps)

df_metad_GEx <- gex_all@meta.data 
df_metad_GEx$raw_cellBC_pdT <- row.names(df_metad_GEx)

setwd(wd_oi)
saveRDS(gex_all,"seurat_obj_GEx_filtered_SP_shuffleBC_revision_20240906.RDS")

write.table(df_metad_GEx,"metadata_GEx_cellBC__SP_shuffleBC_revision_20240906.txt",
            quote=FALSE,row.names=FALSE,sep="\t")

table(df_metad_GEx$orig.ident)


df_metad_GEx <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_GEx/metadata_GEx_cellBC_SP_shuffleBC_revision_20240906.txt",header=TRUE)


# saveRDS(gex_all,"gex_all_reps_no_dim_red_20220525.RDS")



# normalization/dimensional reduction
gex_3 <- gex_all
gex_3 <- NormalizeData(gex_3, normalization.method = "LogNormalize", scale.factor = 10000)
gex_3 <- FindVariableFeatures(gex_3, selection.method = "vst", nfeatures = 1000, verbose = TRUE)
all.genes <- rownames(gex_3)
gex_3 <- ScaleData(gex_3, features = all.genes)
gex_3 <- RunPCA(gex_3, features = VariableFeatures(object = gex_3), verbose = FALSE, npcs = 100)


gex_4 <- gex_3
top_pc <- 50 # before was 40
gex_4 <- FindNeighbors(gex_4, dims = 1:top_pc)
gex_4 <- FindClusters(gex_4, resolution = 0.2) # before was 0.5
gex_4 <- RunUMAP(gex_4, dims = 1:top_pc, n.neighbors = 50, seed.use = 42)

plt0 <- DimPlot(gex_4, reduction = "umap", label=TRUE)
plt1 <- FeaturePlot(object = gex_4,  features = 'nCount_RNA',  pt.size = 0.1,  max.cutoff = 'q99',  ncol = 1,  order = TRUE)
plt2 <- DimPlot(object = gex_4,  group.by = 'orig.ident')
plt0+plt1+plt2

dev.set(4)
DimPlot(object = gex_4,  group.by = 'orig.ident',split.by ='orig.ident')

# df_metad_GEx <- gex_4@meta.data 
# df_metad_GEx$raw_cellBC_pdT <- row.names(df_metad_GEx)
# write.table(df_metad_GEx,"metadata_GEx_cellBC_SP_T7loxShuffle_20231121.txt",
#             quote=FALSE,row.names=FALSE,sep="\t")

# metadata_GEx_cellBC_SP_T7loxShuffle_20231121.txt

