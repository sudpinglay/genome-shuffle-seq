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

date_str <- "20241005"

# read raw data
reps <- c("Lane1_mm10","Lane1_hg38")

df_cell_identities <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_GEx/list_cellBC_hg38_mm10_barnyard_assignment_20240912.txt",
                                 header=TRUE)

df_cells_rearrangements <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/lane1_rearranged_cells.csv",sep=",",
                                      header=FALSE)

df_cells_rearrangements_mm10 <- df_cells_rearrangements %>% filter(V1 %in% (df_cell_identities %>% filter(species_id=="mm10") %>% pull(cellBC)))
df_cells_rearrangements_hg38 <- df_cells_rearrangements %>% filter(V1 %in% (df_cell_identities %>% filter(species_id=="hg38") %>% pull(cellBC)))


rep_oi <- reps[2]
# # # PRE FILTERING
for (rep_oi in reps){
  
  print(rep_oi)
  path_10x_mtx <- sprintf("%s/%s",wd_oi,rep_oi)
  gex_mat <- Read10X(data.dir =path_10x_mtx)
  
  
  # create separate Seurat objects
  gex <- CreateSeuratObject(counts = gex_mat, project = rep_oi, min.cells = 3, min.features = 50)
  
  if (str_detect(rep_oi,"mm10")){
    gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^mt-")
  } else {
    gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^MT-")
  }
  
  if (str_detect(rep_oi,"mm10")){
    mt_thresh <- c(1.3,6)
    UMI_thresh <- 1000
  } else {
    mt_thresh <- c(2,8)
    UMI_thresh <- 1500
  }
  
  
  df_metad_GEx <- gex@meta.data
  df_metad_GEx2 <- data.frame(cellBC=row.names(df_metad_GEx),
                             GEx_UMI=df_metad_GEx$nCount_RNA,
                             mito_frac=df_metad_GEx$percent_mt)
  
  if (str_detect(rep_oi,"mm10")){
    df_metad_GEx2_mm10 <- df_metad_GEx2
  } else {
    df_metad_GEx2_hg38 <- df_metad_GEx2
  }
  # 
  
  fig_QC <- ggplot()+
    stat_bin2d(aes(x=gex$nCount_RNA,y=gex$percent_mt,fill=log(..count..),color=log(..count..)),bins=100)+ 
    geom_point(data=df_metad_GEx2 %>% filter(cellBC %in% (df_cell_identities %>% filter(species_id=="mm10_hg38_doublet") %>% pull(cellBC))),
               aes(x=GEx_UMI,y=mito_frac),color="red",size=0.5,shape=4)+
    # geom_point(data=df_metad_GEx2 %>% filter(cellBC %in% df_cells_rearrangements_mm10$V1),
    #            aes(x=GEx_UMI,y=mito_frac),color="orange",size=0.5)+
    geom_point(data=df_metad_GEx2 %>% filter(cellBC %in% df_cells_rearrangements_hg38$V1),
               aes(x=GEx_UMI,y=mito_frac),color="orange",size=0.5)+
    scale_x_log10(limits=c(100,50000), 
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=expansion(c(0,0)))+
    scale_y_log10(limits=c(0.1,100),
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=expansion(c(0,0)))+
    coord_fixed()+
    annotation_logticks()+
    geom_vline(xintercept=UMI_thresh, color="red", linetype="dashed")+
    geom_hline(yintercept=mt_thresh, color="red", linetype="dashed")+
    theme_cowplot()+
    labs(title="Lane1, human",x="Total transcriptome UMI",y="Mitochondrial fraction (%)")+
    theme(panel.grid.major = element_line(linewidth=0.1,color="grey"))
  
  
  dev.set(4)
  fig_QC
  
  fig_QC_name <- sprintf("mt_vs_RNA_UMI_gex_w_rearranged_doublets_%s_%s.pdf",rep_oi,date_str)
  fig_QC_name <- sprintf("mt_vs_RNA_UMI_gex_w_rearranged_%s_%s.pdf",rep_oi,date_str)
  
  pdf(fig_QC_name,width=6,height=6)
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

# 
# # full data combination for a v2 barnyard
# 
# mt_thresh <- c(1.3,6)
# UMI_thresh <- 1000
# df_metad_GEx2_mm10_2 <- df_metad_GEx2_mm10 %>% transform(category=ifelse(GEx_UMI>UMI_thresh & mito_frac>mt_thresh[1] & mito_frac<mt_thresh[2],"mm10","mm10_non_cell"))
# mt_thresh <- c(2,8)
# UMI_thresh <- 1500
# df_metad_GEx2_hg38_2 <- df_metad_GEx2_hg38 %>% transform(category=ifelse(GEx_UMI>UMI_thresh & mito_frac>mt_thresh[1] & mito_frac<mt_thresh[2],"hg38","hg38_non_cell"))




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
reps <- c("Lane1_mm10","Lane1_hg38")
date_str <- "20240912"
# # # READ ALL W/ DOUBLET SCORES
for (rep_oi in reps){
  
  print(rep_oi)
  path_10x_mtx <- sprintf("%s/%s_filtered_10x_mtx_%s",wd_oi,rep_oi,date_str)
  gex_mat <- Read10X(data.dir =path_10x_mtx)
  
  # create separate Seurat objects
  gex <- CreateSeuratObject(counts = gex_mat, project = rep_oi, min.cells = 3, min.features = 50)
  
  
  if (str_detect(rep_oi,"mm10")){
    gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^mt-")
  } else {
    gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^MT-")
  }
  

  file_doublet <- sprintf("%s/%s_filtered_10x_mtx_%s/scrublet_score_scRNA%s_filtered_10x_mtx_%s.txt",
                          wd_oi,rep_oi,date_str, rep_oi,date_str)
  df_doublet <- read.table(file_doublet,sep="\t",header=TRUE)
  
  rownames(df_doublet) <- df_doublet$cell_bc
  gex$doublet_score <- df_doublet[Cells(gex),2]
  
  assign(sprintf("gex_%s",rep_oi),gex)

}


# # # # threshold on scrublet score
scrublet_thresh <- 0.4
gex_Lane1_mm10_filtered <- subset(gex_Lane1_mm10, subset = doublet_score<scrublet_thresh) 
dim(gex_Lane1_mm10)
dim(gex_Lane1_mm10_filtered)
cells_mm10_pre_scrubled <- colnames(gex_Lane1_mm10)
cells_mm10_post_scrubled <- colnames(gex_Lane1_mm10_filtered)
cells_mm10_scrublet_filtered <- setdiff(cells_mm10_pre_scrubled,cells_mm10_post_scrubled)


gex_Lane1_hg38_filtered <- subset(gex_Lane1_hg38, subset = doublet_score<scrublet_thresh) 
dim(gex_Lane1_hg38)
dim(gex_Lane1_hg38_filtered)
cells_hg83_pre_scrubled <- colnames(gex_Lane1_hg38)
cells_hg38_post_scrubled <- colnames(gex_Lane1_hg38_filtered)
cells_hg38_scrublet_filtered <- setdiff(cells_hg83_pre_scrubled,cells_hg38_post_scrubled)


df_metad_GEx2_mm10 <- gex_Lane1_mm10_filtered@meta.data
df_metad_GEx2_mm10 <- df_metad_GEx2_mm10 %>% transform(cellBC=rownames(df_metad_GEx2_mm10))
df_metad_GEx2_hg38 <- gex_Lane1_hg38_filtered@meta.data
df_metad_GEx2_hg38 <- df_metad_GEx2_hg38 %>% transform(cellBC=rownames(df_metad_GEx2_hg38))

# still need to load all the data to join the under threshold ones
df_raw_GEx_counts <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_GEx/GEx_cell_metadata_Lane1_combined_hg38_mm10_20241005.txt",
                                header=TRUE)

table(df_raw_GEx_counts$species_id)
# need to filter scrublet cells
all_cellBC <- c(df_metad_GEx2_mm10$cellBC,df_metad_GEx2_hg38$cellBC) %>% unique() 

# remove scrublet filtered in EITHER species:
all_cellBC2 <- all_cellBC[!(all_cellBC %in% c(cells_mm10_scrublet_filtered,cells_hg38_scrublet_filtered))]

df_metad_GEx2_mm10_hg38 <- data.frame(cellBC=all_cellBC2) %>% left_join(df_raw_GEx_counts)

# final call: 
df_metad_GEx2_mm10_hg38_2 <- df_metad_GEx2_mm10_hg38 %>%
  transform(species_id2=case_when(species_id=="mm10" & GEx_UMI_mm10/GEx_UMI_hg38<1/0.7 ~ "mm10_hg38_doublet",
                                  species_id=="hg38" & GEx_UMI_hg38/GEx_UMI_mm10<3 ~ "mm10_hg38_doublet",
                                  .default=species_id))

table(df_metad_GEx2_mm10_hg38_2$species_id2)

gex_Lane1_hg38_filtered_final <- subset(gex_Lane1_hg38_filtered,cells=df_metad_GEx2_mm10_hg38_2 %>% filter(species_id2=="hg38") %>% pull(cellBC))
gex_Lane1_mm10_filtered_final <- subset(gex_Lane1_mm10_filtered,cells=df_metad_GEx2_mm10_hg38_2 %>% filter(species_id2=="mm10") %>% pull(cellBC))


saveRDS(gex_Lane1_hg38_filtered_final,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_GEx/seurat_obj_GEx_filtered_SP_shuffleBC_revision_Lane1_final_K562_raw_20241014.RDS")
saveRDS(gex_Lane1_mm10_filtered_final,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_GEx/seurat_obj_GEx_filtered_SP_shuffleBC_revision_Lane1_final_mESC_raw_20241014.RDS")




table(df_metad_GEx2_mm10_hg38$species_id)
#   left_join(df_metad_GEx2_mm10 %>%  select(cellBC,GEx_UMI_mm10=nCount_RNA,mito_frac_mm10=percent_mt) %>% transform(cat_mm10="mm10")) %>%
#   left_join(df_metad_GEx2_hg38 %>%  select(cellBC,GEx_UMI_hg38=nCount_RNA,mito_frac_hg38=percent_mt) %>% transform(cat_hg38="hg38"))

# 4 categories: not a cell, a doublet, a mouse cell, a human cell. 

# df_metad_GEx2_mm10_hg38_w_cell_id <- df_metad_GEx2_mm10_hg38 %>%
#   transform(species_id=case_when(cat_mm10=="mm10" & cat_hg38=="hg38_non_cell" ~ "mm10",
#                                  cat_mm10=="mm10" & cat_hg38=="hg38" ~ "mm10_hg38_doublet",
#                                  cat_mm10=="mm10_non_cell" & cat_hg38=="hg38" ~ "hg38",
#                                  cat_mm10=="mm10_non_cell" & cat_hg38=="hg38_non_cell" ~ "non_cell"))


which(!(df_GEx_mm10_vs_hg38_2$cellBC[df_GEx_mm10_vs_hg38_2$species_id=="hg38"] %in% df_metad_GEx2_mm10_hg38_w_cell_id$cellBC[df_metad_GEx2_mm10_hg38_w_cell_id$species_id=="hg38"]))
df_GEx_mm10_vs_hg38_2$cellBC[df_GEx_mm10_vs_hg38_2$species_id=="hg38"]

cellBC_oi <- "AAACGCTAGGGTGAAA-1"
df_GEx_mm10_vs_hg38_2 %>% filter(cellBC==cellBC_oi)
df_metad_GEx2_mm10_hg38_w_cell_id %>% filter(cellBC==cellBC_oi)
df_raw_GEx_counts  %>% filter(cellBC==cellBC_oi)
df_metad_GEx2_mm10 %>% filter(cellBC==cellBC_oi)
table(df_metad_GEx2_mm10_hg38_w_cell_id$species_id)


dev.set(4)
ggplot(df_metad_GEx2_mm10_hg38_2) + 
  geom_point(aes(x=GEx_UMI_mm10,y=GEx_UMI_hg38,color=species_id2))+
  geom_abline(intercept=log10(3),linetype="dashed")+
  geom_abline(intercept=log10(0.7),linetype="dashed")+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  coord_fixed()+theme_cowplot()+
  annotation_logticks()+
  theme(panel.grid.major = element_line(linewidth=0.1,color="grey"))

ggplot(df_metad_GEx2_mm10_hg38_2 %>% filter(species_id!="non_cell")) + 
  geom_point(aes(x=GEx_UMI_mm10,y=GEx_UMI_hg38,color=species_id2),size=0.5)+
  geom_abline(slope=3,linetype="dashed")+
  geom_abline(slope=0.7,linetype="dashed")+
  scale_x_continuous(limits=c(0,25000),expand=expansion(c(0,0)))+
  scale_y_continuous(limits=c(0,45000),expand=expansion(c(0,0)))+
  
  # scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  # scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  coord_fixed()+theme_cowplot()+
  # annotation_logticks()+
  theme(panel.grid.major = element_line(linewidth=0.1,color="grey"))

df_metad_GEx2_mm10_hg38_w_cell_id %>% filter()

ggplot(df_metad_GEx2_mm10_hg38_w_cell_id %>% filter(species_id=="mm10")) +
  stat_bin(aes(x=GEx_UMI_mm10/GEx_UMI_hg38),geom="step")+
  geom_vline(xintercept=1/0.6,color='red',linetype='dotted')+
  scale_x_log10()

ggplot(df_metad_GEx2_mm10_hg38_w_cell_id %>% filter(species_id=="mm10")) +
  geom_step(aes(x=GEx_UMI_mm10/GEx_UMI_hg38),stat="ecdf")+
  geom_vline(xintercept=1/0.6,color='red',linetype='dotted')+
  scale_x_log10()

df_metad_GEx2_mm10_hg38_w_cell_id %>% filter(species_id=="mm10") %>% dim()
df_metad_GEx2_mm10_hg38_w_cell_id %>% filter(species_id=="mm10") %>% filter(GEx_UMI_mm10/GEx_UMI_hg38<1/0.6) %>% dim()

df_metad_GEx2_mm10_hg38_w_cell_id %>% filter(species_id=="hg38") %>% dim()
df_metad_GEx2_mm10_hg38_w_cell_id %>% filter(species_id=="hg38") %>% filter(GEx_UMI_hg38/GEx_UMI_mm10<3.3) %>% dim()

#gex_2 <- subset(gex, subset = (percent_mt < mt_thresh[2]) & (nCount_RNA >  UMI_thresh ) )
df_species_final <- df_metad_GEx2_mm10_hg38_2 %>% select(-cat_hg38,-cat_mm10,-species_id) %>% rename(species_id2="species_id")
write.table(df_species_final,"FINAL_GEx_cell_metadata_Lane1_combined_hg38_mm10_20241006.txt",
            sep="\t",row.names=FALSE,quote=FALSE)












union_cell_BC <- union(rownames(gex_Lane1_hg38_filtered@meta.data),rownames(gex_Lane1_mm10_filtered@meta.data))

mean(union_cell_BC %in% c(cells_mm10_scrublet_filtered,cells_hg38_scrublet_filtered))

df_GEx_mm10_vs_hg38 <- data.frame(cellBC=union_cell_BC)
df_mm10 <- gex_Lane1_mm10_filtered@meta.data %>% transform(cellBC=rownames(gex_Lane1_mm10_filtered@meta.data)) %>% select(cellBC,GEx_mm10_UMI=nCount_RNA)

df_mm10 %>% filter(cellBC==cellBC_oi)
df_hg38 <- gex_Lane1_hg38_filtered@meta.data %>% transform(cellBC=rownames(gex_Lane1_hg38_filtered@meta.data)) %>% select(cellBC,GEx_hg38_UMI=nCount_RNA)

df_GEx_mm10_vs_hg38 <- df_GEx_mm10_vs_hg38 %>% left_join(df_mm10) %>% left_join(df_hg38)
df_GEx_mm10_vs_hg38[is.na(df_GEx_mm10_vs_hg38)] <- 500


df_GEx_mm10_vs_hg38_2 <- df_GEx_mm10_vs_hg38 %>% transform(species_id=ifelse(GEx_mm10_UMI>500 & GEx_hg38_UMI>500,"mm10_hg38_doublet",ifelse(GEx_mm10_UMI>500,"mm10","hg38")))
table(df_GEx_mm10_vs_hg38_2$species_id)

df_GEx_mm10_vs_hg38_2[df_GEx_mm10_vs_hg38_2==500] <- 0



# write.table(df_GEx_mm10_vs_hg38_2,"list_cellBC_hg38_mm10_barnyard_assignment_20240912.txt",
#             sep="\t",quote=FALSE,row.names=FALSE)

dev.set(4)
ggplot(df_GEx_mm10_vs_hg38) + stat_bin2d(aes(x=GEx_mm10_UMI,y=GEx_hg38_UMI))+
  scale_x_log10()+scale_y_log10()+coord_fixed()


# # # # MERGE ALL IN ONE SEURAT OBJECT
# 
# Idents(gex_sample1_2) <- "Lane1"
# Idents(gex_sample2_2) <- "Lane2"
# 
# saveRDS(gex_sample2_2,"seurat_obj_GEx_filtered_SP_shuffleBC_revision_Lane2_raw_20240906.RDS")
# 
# gex_all <- merge(gex_sample1_2, y = c(gex_sample2_2),
#                        add.cell.ids = reps)
# 
# df_metad_GEx <- gex_all@meta.data 
# df_metad_GEx$raw_cellBC_pdT <- row.names(df_metad_GEx)
# 
# setwd(wd_oi)
# saveRDS(gex_all,"seurat_obj_GEx_filtered_SP_shuffleBC_revision_20240906.RDS")
# 
# write.table(df_metad_GEx,"metadata_GEx_cellBC__SP_shuffleBC_revision_20240906.txt",
#             quote=FALSE,row.names=FALSE,sep="\t")
# 
# table(df_metad_GEx$orig.ident)
# 
# 
# # saveRDS(gex_all,"gex_all_reps_no_dim_red_20220525.RDS")
# 
# 
# 
# # normalization/dimensional reduction
# gex_3 <- gex_all
# gex_3 <- NormalizeData(gex_3, normalization.method = "LogNormalize", scale.factor = 10000)
# gex_3 <- FindVariableFeatures(gex_3, selection.method = "vst", nfeatures = 1000, verbose = TRUE)
# all.genes <- rownames(gex_3)
# gex_3 <- ScaleData(gex_3, features = all.genes)
# gex_3 <- RunPCA(gex_3, features = VariableFeatures(object = gex_3), verbose = FALSE, npcs = 100)
# 
# 
# gex_4 <- gex_3
# top_pc <- 50 # before was 40
# gex_4 <- FindNeighbors(gex_4, dims = 1:top_pc)
# gex_4 <- FindClusters(gex_4, resolution = 0.2) # before was 0.5
# gex_4 <- RunUMAP(gex_4, dims = 1:top_pc, n.neighbors = 50, seed.use = 42)
# 
# plt0 <- DimPlot(gex_4, reduction = "umap", label=TRUE)
# plt1 <- FeaturePlot(object = gex_4,  features = 'nCount_RNA',  pt.size = 0.1,  max.cutoff = 'q99',  ncol = 1,  order = TRUE)
# plt2 <- DimPlot(object = gex_4,  group.by = 'orig.ident')
# plt0+plt1+plt2
# 
# dev.set(4)
# DimPlot(object = gex_4,  group.by = 'orig.ident',split.by ='orig.ident')

#






