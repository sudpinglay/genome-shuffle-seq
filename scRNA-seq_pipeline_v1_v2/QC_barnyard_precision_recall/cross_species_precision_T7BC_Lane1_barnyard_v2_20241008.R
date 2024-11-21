
library(tidyverse)
library(stringr)
library(Matrix)
library(Biostrings)


# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# set of T7 BCs with all UMIs (from species singleton cells and with reads/umi>=15)
# # # # # # # # # # # # # # # # # # # # # # # # # # # 

df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling_Lane1_round1/valid_T7BC_Lane1_RERUN_updated_cells_valid_species_singleton_20241007.txt.gz",
                            header=TRUE)


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

df_T7BC_valid %>% arrange(cBC_CS,CS1_BC)


df_T7BC_valid2 <- df_T7BC_valid %>% select(-BC,-BC_first,-BC_second) %>% transform(BC_pair=paste0(CS1_BC,"__",CS2_BC)) %>% select(-CS1_BC,-CS2_BC)

df_T7BC_valid2 %>% arrange(cBC_CS,BC_pair)


CS1_BC <- df_T7BC_valid2 %>% pull(BC_pair) %>% str_split("__") %>% lapply("[[",1) %>% unlist() 
df_T7BC_valid3 <- df_T7BC_valid2
df_T7BC_valid3$CS1_BC <- CS1_BC
# df_T7BC_valid3 %>% arrange(cBC_CS,CS1_BC)

df_T7BC_valid3 <- df_T7BC_valid3 %>% group_by(cBC_pdT,CS1_BC) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC)) 


# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# SPECIES SPECIFIC BC FROM BULK AMPLICON SEQ
# # # # # # # # # # # # # # # # # # # # # # # # # # # 

# mappable ones: "/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_shuffleBC/20240912_K562_parental_pos_df.tsv"
file_hg38_bulk <- "/Users/jbl/Downloads/df_K562_1063_amplicon_seq_allreps.tsv"

df_shuffleBC_bulk_hg38 <- read.table(file_hg38_bulk,
                                     header=TRUE) %>% filter(averaged_norm_readcount>20)
# df_shuffleBC_bulk_hg38_2 <- df_shuffleBC_bulk_hg38 %>% transform(BC=ifelse(strand_y_1=="bottom-CS1",
#                                                                            paste0(real_barcode2,"_CS1","__",real_barcode1,"_CS2"),
#                                                                            paste0(real_barcode2,"_CS1","__",real_barcode1,"_CS2")))
df_shuffleBC_bulk_hg38_2 <- df_shuffleBC_bulk_hg38 %>% transform(CS1_BC=paste0(real_barcode2,"_CS1")) # BC=paste0(real_barcode2,"_CS1","__",real_barcode1,"_CS2"))

df_shuffleBC_bulk_hg38_3 <- df_shuffleBC_bulk_hg38_2 %>% transform(in_sc= CS1_BC %in% df_T7BC_valid3$CS1_BC)

ggplot(df_shuffleBC_bulk_hg38_3) + stat_bin(aes(x=averaged_norm_readcount,color=in_sc),position="identity",geom="step")+scale_x_log10()


# mappable ones: "/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_shuffleBC/20240912_mESC_parental_pos_df.tsv"
file_mm10_bulk <- "/Users/jbl/Downloads/df_mESC_1063_amplicon_seq_allreps.tsv"
df_shuffleBC_bulk_mm10 <- read.table(file_mm10_bulk,
                                     header=TRUE) %>% filter(averaged_norm_readcount>100)
df_shuffleBC_bulk_mm10_2 <- df_shuffleBC_bulk_mm10  %>% transform(CS1_BC=paste0(real_barcode2,"_CS1"))
# df_shuffleBC_bulk_mm10_2 <- df_shuffleBC_bulk_mm10 %>% transform(BC=ifelse(strand_y_1=="bottom-CS1",
#                                                                            paste0(real_barcode2,"_CS1","__",real_barcode1,"_CS2"),
#                                                                            paste0(real_barcode2,"_CS1","__",real_barcode1,"_CS2")))

df_shuffleBC_bulk_mm10_3 <- df_shuffleBC_bulk_mm10_2 %>% transform(in_sc= CS1_BC  %in% df_T7BC_valid3$CS1_BC)

mean(df_shuffleBC_bulk_mm10_3$in_sc)
mean(df_shuffleBC_bulk_hg38_3$in_sc)

shared_BC <- intersect(df_shuffleBC_bulk_mm10_3$CS1_BC,df_shuffleBC_bulk_hg38_3$CS1_BC)

all_bulk_BC <- union(df_shuffleBC_bulk_mm10_3$CS1_BC,df_shuffleBC_bulk_hg38_3$CS1_BC)

mean(df_T7BC_valid3$CS1_BC %in% all_bulk_BC)
sum(df_T7BC_valid3$UMIs_BC[df_T7BC_valid3$CS1_BC %in% all_bulk_BC])/sum(df_T7BC_valid3$UMIs_BC)
# mean(df_T7BC_CS2$BC1 %in% df_shuffleBC_bulk_hg38$real_barcode1)
# # mean(df_T7BC_CS2$BC1 %in% (df_shuffleBC_bulk_hg38$real_barcode1 %>% DNAStringSet() %>% reverseComplement() %>% as.character()))
# mean(df_T7BC_CS2$BC1 %in% df_shuffleBC_bulk_hg38$real_barcode2)
# mean(df_T7BC_CS2$BC2 %in% df_shuffleBC_bulk_hg38$real_barcode1)
# mean(df_T7BC_CS2$BC2 %in% df_shuffleBC_bulk_hg38$real_barcode2)
# 
# 
# 
# ggplot(df_T7BC_CS2) + stat_bin2d(aes(x=reads_BC,y=UMIs_BC))+
#   scale_x_log10()+scale_y_log10()+geom_abline(intercept=-log10(10))+coord_fixed()

# # # # # # # # # # # # 
# SPECIES ASSIGNMENT 

# more restrictive set of cells without likely doublet FINAL CORRECTION
file_valid_cellBC_v2 <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_GEx/FINAL_GEx_cell_metadata_Lane1_combined_hg38_mm10_20241006.txt"
df_valid_cellBC_all2 <- read.table(file_valid_cellBC_v2,sep="\t",header=TRUE) %>% rename(cellBC="cBC_pdT")

valid_cells_mm10 <- df_valid_cellBC_all2 %>% filter(species_id=="mm10") %>% pull(cBC_pdT)
valid_cells_hg38 <- df_valid_cellBC_all2 %>% filter(species_id=="hg38") %>% pull(cBC_pdT)

length(valid_cells_mm10)
length(valid_cells_hg38)

# # # # # # # # # # # # 

# restrict to species-singleton assigned cells, non shared BC, and BC from the bulk set
df_T7BC_CS_valid_mm10 <- df_T7BC_valid3 %>% filter(cBC_pdT %in% valid_cells_mm10) %>% filter(!(CS1_BC %in% shared_BC) & (CS1_BC %in% all_bulk_BC))
df_T7BC_CS_valid_hg38 <- df_T7BC_valid3 %>% filter(cBC_pdT %in% valid_cells_hg38) %>% filter(!(CS1_BC %in% shared_BC)& (CS1_BC %in% all_bulk_BC))

df_sum_T7 <- (df_T7BC_valid3 %>% filter( (cBC_pdT %in% valid_cells_mm10) | (cBC_pdT %in% valid_cells_hg38))) %>%
  group_by(CS1_BC) %>% summarize(sum_UMI=sum(UMIs_BC)) %>% transform(in_bulk_set=(CS1_BC %in% all_bulk_BC))
# df_sum_T7_mm10 <- df_T7BC_CS_valid_mm10 %>% group_by(BC_pair) %>% summarize(sum_UMI=sum(UMIs_BC))
# df_sum_T7_hg38 <- df_T7BC_CS_valid_hg38 %>% group_by(BC_pair) %>% summarize(sum_UMI=sum(UMIs_BC))
# 
hi_BC <- df_sum_T7  %>% filter(sum_UMI>10) %>% pull(CS1_BC)

ggplot(df_sum_T7) + stat_bin(aes(x=sum_UMI,color=in_bulk_set),position="identity",geom="step")+scale_x_log10()+scale_y_log10()

ggplot(df_sum_T7) 

# hi_BC_mm10 <- df_sum_T7_mm10 %>% filter(sum_UMI>10) %>% pull(BC_pair)
# hi_BC_hg38 <- df_sum_T7_hg38 %>% filter(sum_UMI>10) %>% pull(BC_pair)


# proportion of incorrect detection as a function of UMI threshold
UMI_threshs <- seq(1,10)
df_prop_detections_mESC <- data.frame()
df_prop_detections_K562 <- data.frame()

n_tot_hg38_umis <- df_T7BC_CS_valid_hg38 %>% filter(UMIs_BC>=1) %>% pull(UMIs_BC) %>% sum()
n_tot_mm10_umis <- df_T7BC_CS_valid_mm10 %>% filter(UMIs_BC>=1) %>% pull(UMIs_BC) %>% sum()


df_filtered_cells_all <- data.frame()

for (UMI_thresh_oi in UMI_threshs){
  # df_filtered_mm10 <- df_T7BC_CS_valid_mm10 %>% filter(UMIs_BC>=UMI_thresh_oi) %>% filter(BC_pair %in% union(hi_BC_mm10,hi_BC_hg38))
  # df_filtered_hg38 <- df_T7BC_CS_valid_hg38 %>% filter(UMIs_BC>=UMI_thresh_oi) %>% filter(BC_pair %in% union(hi_BC_mm10,hi_BC_hg38))
  # 
  df_filtered_mm10 <- df_T7BC_CS_valid_mm10 %>% filter(UMIs_BC>=UMI_thresh_oi) #%>% filter(BC_pair %in% hi_BC)
  df_filtered_hg38 <- df_T7BC_CS_valid_hg38 %>% filter(UMIs_BC>=UMI_thresh_oi) #%>% filter(BC_pair %in% hi_BC)
  
  df_filtered_mm10 <- df_filtered_mm10 %>% transform(in_hg38_list=CS1_BC %in% df_shuffleBC_bulk_hg38_3$CS1_BC,
                                                     # frac_umi_in_hg38_list=sum(UMIs_BC[BC_pair %in% df_shuffleBC_bulk_hg38_3$BC])/sum(UMIs_BC),
                                                     in_mm10_list=CS1_BC %in% df_shuffleBC_bulk_mm10_3$CS1_BC)
                                                     # frac_umi_in_mm10_list=sum(UMIs_BC[BC_pair %in% df_shuffleBC_bulk_mm10_3$BC])/sum(UMIs_BC))
  
  df_filtered_hg38 <- df_filtered_hg38 %>% transform(in_hg38_list=CS1_BC %in% df_shuffleBC_bulk_hg38_3$CS1_BC,
                                                     # frac_umi_in_hg38_list=sum(UMIs_BC[BC_pair %in% df_shuffleBC_bulk_hg38_3$BC])/sum(UMIs_BC),
                                                     in_mm10_list=CS1_BC %in% df_shuffleBC_bulk_mm10_3$CS1_BC)
                                                     # frac_umi_in_mm10_list=sum(UMIs_BC[BC_pair %in% df_shuffleBC_bulk_mm10_3$BC])/sum(UMIs_BC))
  
  
  df_filtered_mm10_cells <- df_filtered_mm10 %>% group_by(cBC_pdT) %>% summarize(mm10_UMI=sum(UMIs_BC[in_mm10_list]),
                                                                                 hg38_UMI=sum(UMIs_BC[in_hg38_list]))
  
  df_filtered_hg38_cells <- df_filtered_hg38 %>% group_by(cBC_pdT) %>% summarize(mm10_UMI=sum(UMIs_BC[in_mm10_list]),
                                                                                 hg38_UMI=sum(UMIs_BC[in_hg38_list]))
  
  df_filtered_cells <- rbind(df_filtered_mm10_cells %>% transform(species="mouse"),
                             df_filtered_hg38_cells %>% transform(species="human"))
  
  df_filtered_cells_all <- rbind(df_filtered_cells_all,df_filtered_cells %>% transform(UMI_thresh=UMI_thresh_oi))
  
  
  df_prop_detections_mESC <- rbind(df_prop_detections_mESC,data.frame(UMI_thresh=UMI_thresh_oi,
                                                                      n_BCs=length(df_filtered_mm10$BC_pair),
                                                                      n_umis=sum(df_filtered_mm10$UMIs_BC),
                                                                      prop_n_mm10=mean(df_filtered_mm10$in_mm10_list),
                                                                      prop_umi_mm10=sum(df_filtered_mm10$UMIs_BC[df_filtered_mm10$in_mm10_list])/sum(df_filtered_mm10$UMIs_BC),
                                                                      prop_n_hg38=mean(df_filtered_mm10$in_hg38_list),
                                                                      prop_umi_hg38=sum(df_filtered_mm10$UMIs_BC[df_filtered_mm10$in_hg38_list])/sum(df_filtered_mm10$UMIs_BC)))
  
  df_prop_detections_K562 <- rbind(df_prop_detections_K562,data.frame(UMI_thresh=UMI_thresh_oi,
                                                                      n_BCs=length(df_filtered_hg38$BC_pair),
                                                                      n_umis=sum(df_filtered_hg38$UMIs_BC),
                                                                      prop_n_mm10=mean(df_filtered_hg38$in_mm10_list),
                                                                      prop_umi_mm10=sum(df_filtered_hg38$UMIs_BC[df_filtered_hg38$in_mm10_list])/sum(df_filtered_hg38$UMIs_BC),
                                                                      prop_n_hg38=mean(df_filtered_hg38$in_hg38_list),
                                                                      prop_umi_hg38=sum(df_filtered_hg38$UMIs_BC[df_filtered_hg38$in_hg38_list])/sum(df_filtered_hg38$UMIs_BC)))
}


color_mESC <- "#E1812B"
color_K562 <- "#3175A1"


p_count <- 1
dev.set(4)


ggplot(df_filtered_cells_all) + stat_bin2d(aes(x=mm10_UMI+p_count,y=hg38_UMI+p_count),bins=20)+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits=c(0.3,300))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                limits=c(0.3,300))+
  annotation_logticks()+mytheme()+coord_fixed()+
  scale_color_manual(values=c(color_K562,color_mESC))+
  facet_grid(cols=vars(UMI_thresh),rows=vars(species))
  
plt_T7_barnyard <- ggplot(df_filtered_cells_all %>% filter(UMI_thresh<=6)) + 
  # geom_point(aes(x=mm10_UMI+p_count,y=hg38_UMI+p_count,color=species),size=0.3)+
  geom_jitter(aes(x=mm10_UMI+p_count,y=hg38_UMI+p_count,color=species),size=0.5,alpha=0.05,width=0.03,height=0.03,shape=16)+
  
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits=c(0.3,300))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                limits=c(0.3,300))+
  annotation_logticks()+mytheme()+coord_fixed()+
  scale_color_manual(values=c(color_K562,color_mESC))+
  # facet_wrap(~UMI_thresh,nrow=1)+
  # facet_grid(cols=vars(UMI_thresh),rows=vars(species))
  facet_grid(cols=vars(UMI_thresh))+
  labs(x="Sum UMIs to mESC barcodes (+1)",
       y="Sum UMIs to\nK562 barcodes (+1)")

dev.set(4)
plt_T7_barnyard

pdf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/T7_barnyard_all_UMI_thresh_20240511.pdf",width=8,height=2)
pdf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/T7_barnyard_6lower_UMI_thresh_20240511.pdf",width=6.5,height=2.5)

print(plt_T7_barnyard)
dev.off()




mytheme <- function() {
  theme_cowplot() +
    theme(
      # Text elements
      # text = element_text(family = "Arial", color = "black"),
      # plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      # plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      # panel.grid.major.x=element_line(linewidth=0.1,color="grey"),
      panel.grid.major=element_line(linewidth=0.1,color="grey"),
      
      strip.background = element_blank(),
      strip.text = element_text(hjust=0,size=10),
      
      # Legend
      legend.position = "none",
  
    )
}




dev.set(4)
plt_mESC <- ggplot(df_prop_detections_mESC) + 
  geom_point(aes(x=UMI_thresh,y=n_umis/n_tot_mm10_umis),color="black")+
  geom_line(aes(x=UMI_thresh,y=n_umis/n_tot_mm10_umis),color="black")+
  geom_point(aes(x=UMI_thresh,y=prop_umi_mm10),color=color_mESC)+
  geom_line(aes(x=UMI_thresh,y=prop_umi_mm10),color=color_mESC)+
  geom_point(aes(x=UMI_thresh,y=prop_umi_hg38),color=color_K562)+
  geom_line(aes(x=UMI_thresh,y=prop_umi_hg38),color=color_K562)+

  
  scale_y_continuous(limits=c(0,1.03),expand=expansion(c(0,0)))+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  
  labs(x="UMI threshold (inclusive) for detection",y="Proportion of UMIs detected",title="mESC")+
  mytheme()+
  theme(panel.grid.major = element_line(color="grey",linewidth=0.1))

plt_mESC_umi <-ggplot(df_prop_detections_mESC,aes(x=UMI_thresh,y=n_umis)) + 
  geom_point()+
  geom_line()+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits=c(1e-2,1)*1e6)+
  annotation_logticks(sides="l")+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  mytheme()+
  labs(x="UMI threshold (inclusive) for detection",y="Number of UMIs detected")+
  theme(panel.grid.major = element_line(color="grey",linewidth=0.1),)

# plt_mESC_umi

plt_K562 <- ggplot(df_prop_detections_K562) + 
  geom_point(aes(x=UMI_thresh,y=n_umis/n_tot_hg38_umis),color="black")+
  geom_line(aes(x=UMI_thresh,y=n_umis/n_tot_hg38_umis),color="black")+
  geom_point(aes(x=UMI_thresh,y=prop_umi_mm10),color=color_mESC)+
  geom_line(aes(x=UMI_thresh,y=prop_umi_mm10),color=color_mESC)+
  geom_point(aes(x=UMI_thresh,y=prop_umi_hg38),color=color_K562)+
  geom_line(aes(x=UMI_thresh,y=prop_umi_hg38),color=color_K562)+
  
  scale_y_continuous(limits=c(0,1.03),expand=expansion(c(0,0)))+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  mytheme()+
  labs(x="UMI threshold (inclusive) for detection",y="Proportion of UMIs detected",title="K562")+
  theme(panel.grid.major = element_line(color="grey",linewidth=0.1),
        axis.title=element_blank())

plt_K562_umi <-ggplot(df_prop_detections_K562,aes(x=UMI_thresh,y=n_umis)) + 
  geom_point()+
  geom_line()+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits=c(1e-2,1)*1e6)+
  annotation_logticks(sides="l")+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  mytheme()+
  labs(x="UMI threshold (inclusive) for detection",y="Number of UMIs detected")+
  theme(panel.grid.major = element_line(color="grey",linewidth=0.1),
        axis.title=element_blank())

plt_all <- plt_mESC+plt_K562+plt_mESC_umi+plt_K562_umi
plt_all

pdf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/cross_detection_analysis_Lane1_K562_mESC_v2_20241008.pdf",width=6,height=5)
print(plt_all)
dev.off()


mm10_UMI_thresh <- 20
hg38_UMI_thresh <- 20
ggplot(df_sum_T7_mm10) + stat_bin(aes(x=sum_UMI),geom="step")+scale_x_log10()+scale_y_log10()+geom_vline(aes(xintercept=mm10_UMI_thresh),color="red",linewidth=0.2)





# mean(df_T7BC_CS2$)

table(df_shuffleBC_bulk_hg38$strand_y_1)



ggplot(df_sum_T7_hg38) + stat_bin(aes(x=sum_UMI),geom="step")+scale_x_log10()+scale_y_log10()+geom_vline(aes(xintercept=hg38_UMI_thresh),color="red",linewidth=0.2)

hi_count_BC_mm10 <- df_sum_T7_mm10 %>% filter(sum_UMI>mm10_UMI_thresh) %>% pull(BC)
hi_count_BC_hg38 <- df_sum_T7_hg38 %>% filter(sum_UMI>hg38_UMI_thresh) %>% pull(BC)

length(intersect(hi_count_BC_mm10,hi_count_BC_hg38))/length(union(hi_count_BC_mm10,hi_count_BC_hg38))



df_valid_cellBC_all <-  df_valid_cellBC_all %>% 
  transform(raw_cellBC_pdT=str_split(raw_cellBC_pdT,"_") %>% lapply("[[",2) %>% unlist())

# restrict to sample of interest
df_valid_cellBC_all2 <- df_valid_cellBC_all %>% filter(orig.ident==sample_oi)

                     
df_T7BC_CS3 <- df_T7BC_CS2 %>% transform(valid_cBC=(cBC_pdT %in% df_valid_cellBC_all2$raw_cellBC_pdT))

df_sum_T7 <- df_T7BC_CS3 %>% group_by(BC) %>% summarize(sum_UMI=sum(UMIs_BC)) 

hi_BC <- df_sum_T7 %>% filter(sum_UMI>10) %>% pull(BC)
list_all_cBC <- df_T7BC_CS3 %>% pull(cBC_pdT) %>% unique()
empty_cBC <- list_all_cBC[!(list_all_cBC %in% df_valid_cellBC_all2$raw_cellBC_pdT)]

df_empty_cBC <- data.frame(cBC_pdT=empty_cBC) %>% left_join(df_T7BC_CS3 %>% select(cBC_pdT,UMIs_BC,BC))
ggplot(df_sum_T7) + stat_bin(aes(x=sum_UMI),geom="step")+scale_x_log10()+scale_y_log10()

df_T7BC_CS3_empty <- df_T7BC_CS3 %>% filter(cBC_pdT %in% empty_cBC) 
df_BC_av_empty <- df_T7BC_CS3_empty %>% group_by(BC) %>% summarize(mean_UMI_per_empty=sum(UMIs_BC)/length(empty_cBC)) %>% arrange(desc(mean_UMI_per_empty))

df_tot_shuffle_BC <- df_T7BC_CS3 %>% group_by(raw_cellBC_pdT=cBC_pdT,valid_cBC) %>% summarize(sum_shuffleBC_umi=sum(UMIs_BC))

df_valid_cellBC_all2 <- df_valid_cellBC_all %>% left_join(df_tot_shuffle_BC)
ggplot(df_valid_cellBC_all2) + stat_bin2d(aes(x=nCount_RNA,y=sum_shuffleBC_umi))+
  scale_x_log10()+scale_y_log10()



dev.set(4)
ggplot(df_tot_shuffle_BC) + stat_bin(aes(x=sum_shuffleBC_umi,color=valid_cBC),position="identity",geom="step")+
  scale_x_log10()+scale_y_log10()

df_T7BC_CS3 %>% filter(valid_cBC) %>% dim()

# loose filter on counts
df_T7BC_CS_valid <- df_T7BC_CS3 %>% filter(valid_cBC) %>% select(-valid_cBC) %>% filter(reads_BC/UMIs_BC>=9)
  # filter(UMIs_BC>=2 | reads_BC/UMIs_BC>=8)


write.table(df_T7BC_CS_valid,
            sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling/valid_T7BC_%s_20240907.txt",sample_oi),
            sep="\t",row.names=FALSE,quote=FALSE)

dev.set(4)
ggplot(df_T7BC_CS_valid) + geom_bin_2d(aes(x=reads_BC,y=UMIs_BC))+
  scale_x_log10()+scale_y_log10()+coord_fixed()+labs(title=sample_oi)

