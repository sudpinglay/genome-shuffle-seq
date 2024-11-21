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


setwd("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/gex_prepro_v2")

wd_oi <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/gex_prepro_v2"
date_str <- "20231220"
# read raw data
reps <- c("Cre","Parental")



# # # PRE FILTERING
for (rep_oi in reps){
  
  print(rep_oi)
  path_10x_mtx <- sprintf("%s/%s",wd_oi,rep_oi)
  gex_mat <- Read10X(data.dir =path_10x_mtx)
  
  
  # create separate Seurat objects
  gex <- CreateSeuratObject(counts = gex_mat, project = rep_oi, min.cells = 3, min.features = 50)
  gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^mt-")
  
  mt_thresh <- c(1,12)
  UMI_thresh <- 1000
  fig_QC <- ggplot()+stat_bin2d(aes(x=gex$nCount_RNA,y=gex$percent_mt,fill=log(..count..),color=log(..count..)),bins=100)+
    scale_x_log10(limits=c(30,50000))+
    scale_y_log10(limits=c(0.05,100))+
    annotation_logticks()+
    geom_vline(xintercept=UMI_thresh, color="red", linetype="dashed")+
    geom_hline(yintercept=mt_thresh, color="red", linetype="dashed")
  
  
  fig_QC_name <- sprintf("mt_vs_RNA_UMI_gex_%s_%s.pdf",rep_oi,date_str)
  
  pdf(fig_QC_name,width=4,height=4)
  print(fig_QC)
  dev.off()
  
  
  # subsetting on cells
  gex_2 <- subset(gex, subset = (percent_mt < mt_thresh[2]) & (percent_mt > mt_thresh[1]) & (nCount_RNA >  UMI_thresh ) )
  
  print(dim(gex))
  print(dim(gex_2))
  
  write10xCounts(sprintf("%s/%s_filtered_10x_mtx_%s/",wd_oi,rep_oi,date_str),
                 gex_2@assays$RNA@counts,
                 barcodes = colnames(gex_2),
                 gene.id = rownames(gex_2))
  
}

#gex_2 <- subset(gex, subset = (percent_mt < mt_thresh[2]) & (nCount_RNA >  UMI_thresh ) )




rep_oi <- "Cre"

print(rep_oi)
path_10x_mtx <- sprintf("%s/%s",wd_oi,rep_oi)
gex_mat <- Read10X(data.dir =path_10x_mtx)


# create separate Seurat objects
gex <- CreateSeuratObject(counts = gex_mat, project = rep_oi, min.cells = 3, min.features = 50)
gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^mt-")

setwd("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/gex_prepro_v2/")

df_raw_metad_Cre <- gex@meta.data %>% transform(cBC_pdT=row.names(gex@meta.data))
df_raw_metad_Parental <- gex@meta.data %>% transform(cBC_pdT=row.names(gex@meta.data))

df_raw_gex <- rbind(df_raw_metad_Cre %>% transform(origin="Cre"),
                    df_raw_metad_Parental %>% transform(origin="Parental"))

library(scales)
mt_thresh <- c(1,12)
UMI_thresh <- 1000
fig_QC <- ggplot(df_raw_gex)+stat_bin2d(aes(x=nCount_RNA,y=percent_mt,fill=log(..count..),color=log(..count..)),bins=100)+
  scale_x_log10(limits=c(30,50000),breaks=c(10,100,1000,10000),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(limits=c(0.2,70),breaks=c(1,10),labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks()+
  facet_wrap(~orig.ident)+
  geom_vline(xintercept=UMI_thresh, color="red", linetype="dashed")+
  geom_hline(yintercept=mt_thresh, color="red", linetype="dashed")+
  theme_cowplot()+
  labs(x="Transcriptome UMI/cell",y="Mitochondrial RNA fraction (%)")+
  theme(strip.background=element_blank())

dev.new()
pdf("QC_plot_raw_GEx_metadata_Cre_vs_Parental_20240114.pdf",width=5,height=3)
print(fig_QC)
dev.off()


write.table(gex@meta.data %>% transform(cBC_pdT=row.names(gex@meta.data)),"raw_GEx_data_sample1_Cre_20240114.txt",
            row.names=FALSE, quote=FALSE, sep="\t")

mt_thresh <- c(1,12)
UMI_thresh <- 750
fig_QC <- ggplot()+stat_bin2d(aes(x=gex$nCount_RNA,y=gex$percent_mt,fill=log(..count..),color=log(..count..)),bins=100)+
  scale_x_log10(limits=c(30,50000))+
  scale_y_log10(limits=c(0.05,100))+
  annotation_logticks()+
  geom_vline(xintercept=UMI_thresh, color="red", linetype="dashed")+
  geom_hline(yintercept=mt_thresh, color="red", linetype="dashed")









# read raw data
reps <- c("Cre","Parental")

# # # READ ALL W/ DOUBLET SCORES
for (rep_oi in reps){
  
  print(rep_oi)
  path_10x_mtx <- sprintf("%s/%s_filtered_10x_mtx_%s",wd_oi,rep_oi,date_str)
  gex_mat <- Read10X(data.dir =path_10x_mtx)
  
  # create separate Seurat objects
  gex <- CreateSeuratObject(counts = gex_mat, project = rep_oi, min.cells = 3, min.features = 50)
  gex[["percent_mt"]] <- PercentageFeatureSet(gex, pattern = "^mt-")
  

  file_doublet <- sprintf("%s/%s_filtered_10x_mtx_%s/scrublet_score_scRNA%s_filtered_10x_mtx_%s.txt",
                          wd_oi,rep_oi,date_str, rep_oi,date_str)
  df_doublet <- read.table(file_doublet,sep="\t",header=TRUE)
  
  rownames(df_doublet) <- df_doublet$cell_bc
  gex$doublet_score <- df_doublet[Cells(gex),2]
  
  assign(sprintf("gex_%s",rep_oi),gex)

}


# # # # threshold on scrublet score
scrublet_thresh <- 0.4
gex_sample1_2 <- subset(gex_Cre, subset = doublet_score<scrublet_thresh) 
gex_sample2_2 <- subset(gex_Parental, subset = doublet_score<scrublet_thresh) 



# # # MERGE ALL IN ONE SEURAT OBJECT

Idents(gex_sample1_2) <- "Cre"
Idents(gex_sample2_2) <- "Parental"



gex_all <- merge(gex_sample1_2, y = c(gex_sample2_2),
                       add.cell.ids = reps)

df_metad_GEx <- gex_all@meta.data 
df_metad_GEx$raw_cellBC_pdT <- row.names(df_metad_GEx)
write.table(df_metad_GEx,"metadata_GEx_cellBC_SP_T7loxShuffle_combined_20231220.txt",
            quote=FALSE,row.names=FALSE,sep="\t")

table(df_metad_GEx$orig.ident)


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




rownames(gex_4) %>% str_detect("lenti") %>% which
rownames(gex_4)[16718]
# cluster representation:
table_freq <- table(gex_4$seurat_clusters,gex_4$orig.ident)

df_exp_lenti <- gex_4@meta.data
df_exp_lenti$cellBC <- rownames(df_exp_lenti)

lenti_dCas9_exp <- gex_4@assays$RNA@counts[16718,]

df_exp_lenti$dCas9_exp <- lenti_dCas9_exp



dev.set(4)
ggplot(df_exp_lenti) + stat_bin_2d(aes(x=nCount_RNA,y=dCas9_exp))+
  scale_x_log10()+facet_wrap(~seurat_clusters)


ggplot(df_exp_lenti) + geom_violin(aes(x=seurat_clusters,y=(dCas9_exp+0.03)/(nCount_RNA/1000)))

ggplot(df_exp_lenti) + geom_boxplot(aes(x=seurat_clusters,y=dCas9_exp))

ggplot(df_exp_lenti) + geom_boxplot(aes(x=seurat_clusters,y=dCas9_exp+0.1))+scale_y_log10()


df_exp_lenti %>% group_by(seurat_clusters) %>% summarize(frac_dCas9=mean(dCas9_exp>0),
                                                         mean_dCas9=mean(dCas9_exp),
                                                         mean_norm_dCas9=mean(dCas9_exp/(nCount_RNA/10000)),
                                                         mean_GEx_UMI=mean(nCount_RNA),
                                                         median_GEx_UMI=median(nCount_RNA))

ggplot(df_exp_lenti) + geom_violin(aes(x=seurat_clusters,y=(dCas9_exp+0.03)/(nCount_RNA/1000)))+
  scale_y_log10()
# 
# library(ggbeeswarm)
# ggplot(df_exp_lenti) + geom_beeswarm(aes(x=seurat_clusters,y=(dCas9_exp+0.05)/(nCount_RNA/1000)))+
#   scale_y_log10()






raw_cellBC_pdT <- str_split(Cells(gex_4),"_") %>% lapply("[[",2) %>% substr(1,16)
df_metadata_cells <- data.frame(cellBC=Cells(gex_4),
                                raw_cellBC_pdT,
                                rep_id=gex_4$orig.ident,
                                gex_UMI=gex_4$nCount_RNA,
                                gex_percent_mt=round(gex_4$percent_mt,2),
                                cluster_id=gex_4$seurat_clusters)

write.table(df_metadata_cells,sprintf("metadata_filtered_cells_gex_scv2_%s.txt",date_str),
            sep="\t",row.names=FALSE,quote=FALSE)


                         





 DimPlot(gex_4, cells=Cells(gex_4)[gex_4$percent_mt<1], reduction = "umap", label=TRUE)+
  DimPlot(gex_4, cells=Cells(gex_4)[gex_4$percent_mt>1], reduction = "umap", label=TRUE)

gex_4$low_mt_flag <- gex_4$percent_mt<1








FeaturePlot(object = gex_4,  features = 'percent_mt',  pt.size = 0.1,  min.cutoff = 'q10',  ncol = 1,  order = TRUE)
FeaturePlot(object = gex_4,  features = 'low_mt_flag',  pt.size = 0.1,  ncol = 1,  order = TRUE)


FeaturePlot(object = gex_4,  cells=Cells(gex_4)[gex_4$percent_mt<1],  pt.size = 0.1,  ncol = 1,  order = TRUE)


FeaturePlot(object = gex_4,  features = c("Lama1","Sox17","Tubb2b","Col5a1","Col1a1","Sox2","Gata4"),  pt.size = 0.1,  max.cutoff = 'q95',  order = TRUE)

# highlight cells with low mt content? 

plt0+plt1

table(gex_scMPRA_4$seurat_clusters)


marker1 <- FindMarkers(object=gex_scMPRA_4, 
            ident.1=1,
            test.use = "negbinom",
            logfc.threshold = 1)
marker1 %>% arrange(desc(avg_log2FC)) %>% head(n=20)
marker1 %>% arrange(desc(avg_log2FC)) %>% tail(n=20)


marker2 <- FindMarkers(object=gex_scMPRA_4, 
                       ident.1=2,
                       test.use = "negbinom",
                       logfc.threshold = 1)
marker2 %>% arrange(desc(avg_log2FC)) %>% head(n=20)

marker0 <- FindMarkers(object=gex_scMPRA_4, 
                       ident.1=0,
                       test.use = "negbinom",
                       logfc.threshold = 1)
marker0 %>% arrange(desc(avg_log2FC)) %>% head(n=20)


# tentative
cell_line <- c("K562","HEK293","HepG2")

gex_scMPRA_4$cell_line <- NA
gex_scMPRA_4$cell_line[gex_scMPRA_4$seurat_clusters==0] <- "K562"
gex_scMPRA_4$cell_line[gex_scMPRA_4$seurat_clusters==1] <- "HEK293"
gex_scMPRA_4$cell_line[gex_scMPRA_4$seurat_clusters==2] <- "HepG2"


df_summary <- data.frame(cBC=Cells(gex_scMPRA_4),
                         gex_UMI=gex_scMPRA_4$nCount_RNA,
                         cell_line=gex_scMPRA_4$cell_line)

write.table(df_summary,"repA_gex_scMPRA_valid_cBC_20220202.txt", sep="\t", row.names=FALSE, quote=FALSE)
DimPlot(seur_obj, reduction = "umap", label=TRUE, split.by = 'orig.ident')


