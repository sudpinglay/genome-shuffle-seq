
# # # # # # # # # # # # # # # # # # # # 
# get the raw counts with umi = 1 too
# # # # # # # # # # # # # # # # # # # # 

# 10x T7 data from 'valid' cellBC: 
df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/valid_T7BC_Cre_v3_CS_20231220.txt.gz",
                            header=TRUE)
# # parental cells: 
# df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/valid_T7BC_Parental_v3_CS_20231220.txt.gz",
#                             header=TRUE)

# # paired BCs
# BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split("_CS")
# df_T7BC_valid$BC_pair <- BC_all %>% lapply("[[",1) %>% unlist()
# df_T7BC_valid$capture_direction <- paste0("CS",BC_all %>% lapply("[[",2) %>% unlist())
# df_T7BC_valid <- df_T7BC_valid %>% select(-BC)
# df_T7BC_valid2_1umi <- df_T7BC_valid %>% group_by(cBC_pdT,BC_pair) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC))
# df_T7BC_valid2 <- df_T7BC_valid2_1umi %>% filter(UMIs_BC>=2)




# considering BC separately:
BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split("_")
df_T7BC_valid$BC1 <- BC_all %>% lapply("[[",1) %>% unlist()
df_T7BC_valid$BC2 <- BC_all %>% lapply("[[",2) %>% unlist()
df_T7BC_valid$capture_direction <- BC_all %>% lapply("[[",3) %>% unlist()
df_T7BC_valid <- df_T7BC_valid %>% select(-BC)
df_T7BC_valid2_1umi <- df_T7BC_valid %>% group_by(cBC_pdT,BC1,BC2) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC)) %>% data.frame()
df_T7BC_valid2 <- df_T7BC_valid2_1umi %>% filter(UMIs_BC>=2) %>% data.frame()













# # # # # # # # # # # # # # # # # # # # 
# generating different assignments with different thresholds and repeating the PR analysis 
# # # # # # # # # # # # # # # # # # # # 


# df_mapping_clonotype_final_round1_3 <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/FINAL_ASSIGNMENTS_round2_Cre_cells_to_clonotypes_denovo3umi_round1_round2_bulk_v2_separate_BC_20240112.txt",
#                                                   header=TRUE)


df_mapping_to_clonotype <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/mapping_Cre_T7BC_CS_denovo3umi_round2_w_bulk_v2_separate_BC_20240112.txt",
                                      header=TRUE)

thresh_capt_low <- 0.1
thresh_capt_hi <- 0.3
# thresh_precision <- 0.75
# 
thresh_precision <- 0.5

# # best precision (purity)
df_mapping_to_clonotype_best <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_10xBC1_in_clone,n=1,with_ties=FALSE)
df_mapping_to_clonotype_best2 <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_10xBC1_in_clone,n=2,with_ties=FALSE)

df_mapping_to_clonotype_best2_only <- anti_join(df_mapping_to_clonotype_best2,df_mapping_to_clonotype_best)
colnames(df_mapping_to_clonotype_best2_only)[c(2,3,4,5,6)] <- c("consensus_clone_id_2ndBest","frac_clone_in_10x1_2ndBest","frac_clone_in_10x2_2ndBest","frac_10xBC1_in_clone_2ndBest","frac_10xBC2_in_clone_2ndBest")
df_mapping_to_clonotype_best_w_2nd_best <- df_mapping_to_clonotype_best %>% left_join(df_mapping_to_clonotype_best2_only)
colnames(df_mapping_to_clonotype_best_w_2nd_best)[c(2,3,4,5,6)] <- c("clone_id","recall_BC1","recall_BC2","precision_BC1","precision_BC2")
colnames(df_mapping_to_clonotype_best_w_2nd_best)[c(7,8,9,10,11)] <- c("clone_id_2nd","recall_BC1_2nd","recall_BC2_2nd","precision_BC1_2nd","precision_BC2_2nd")

# categorization

df_cat_no_capt <-  df_mapping_to_clonotype_best %>% filter(!(frac_clone_in_10x1>0 | frac_10xBC1_in_clone>0)) %>% 
  transform(assignment_category="no_capture")

df_cat_lo_capt <- df_mapping_to_clonotype_best %>% filter((frac_clone_in_10x1>0 | frac_10xBC1_in_clone>0) & frac_clone_in_10x1<thresh_capt_low) %>% 
  transform(assignment_category="low_capture")

df_cat_mid_capt_hi_purity <- df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x1>=thresh_capt_low & frac_clone_in_10x1<thresh_capt_hi & frac_10xBC1_in_clone>=thresh_precision)%>% 
  transform(assignment_category="mid_capture_high_purity")

df_cat_hi_capt_hi_purity <- df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x1>=thresh_capt_hi & frac_10xBC1_in_clone>=thresh_precision) %>% 
  transform(assignment_category="high_capture_high_purity")

df_cat_lo_purity <- df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x1>=thresh_capt_low & frac_10xBC1_in_clone<thresh_precision) %>% 
  transform(assignment_category="low_purity")

df_mapping_to_clonotype_w_category <- rbind(df_cat_no_capt,
                                            df_cat_lo_capt,
                                            df_cat_mid_capt_hi_purity,
                                            df_cat_hi_capt_hi_purity,
                                            df_cat_lo_purity)

# table(df_mapping_to_clonotype_w_category$assignment_category)


# refinement based on second best match
df_mapping_clonotype_final_round1 <- df_mapping_to_clonotype_best_w_2nd_best %>% left_join(df_mapping_to_clonotype_w_category %>% select(cellBC,assignment_category))

df_mapping_clonotype_final_round1_2 <- df_mapping_clonotype_final_round1 #%>% transform(ifelse(assignment_category=="high_capture_high_purity"))
df_mapping_clonotype_final_round1_2$assignment_category[df_mapping_clonotype_final_round1_2$assignment_category=="high_capture_high_purity" &
                                                          df_mapping_clonotype_final_round1_2$recall_BC1_2nd>=0.1] <- "high_capture_high_purity_too_high_recall_2nd"
df_mapping_clonotype_final_round1_2$assignment_category[df_mapping_clonotype_final_round1_2$assignment_category=="mid_capture_high_purity" &
                                                          df_mapping_clonotype_final_round1_2$recall_BC1_2nd>=0.1] <- "mid_capture_high_purity_too_high_recall_2nd"

df_mapping_clonotype_final_round1_3 <- df_mapping_clonotype_final_round1_2 %>% 
  transform(valid_assignment=ifelse(assignment_category %in% c("high_capture_high_purity","mid_capture_high_purity"),TRUE,FALSE))



df_BC_clonotypes <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_round2/FINAL_CLONOTYPES_v2_denovo3umi_round1_round2_w_bulk_joined_BC_pairs_20240106.txt",
                               header=TRUE)
df_BC_clonotypes <- df_BC_clonotypes %>% transform(
  BC_pair=str_replace(BC_pair,"-","_"),
  BC1=str_split(BC_pair,"-") %>% lapply("[[",1) %>% unlist(),
  BC2=str_split(BC_pair,"-") %>% lapply("[[",2) %>% unlist())



# precision recall analysis on final round1 set. 
umi_ct_thresh <- seq(1,10)
df_precision_recall_high_conf_clones <- data.frame()

df_valid_assignment <- df_mapping_clonotype_final_round1_3 %>% filter(valid_assignment)
for (id in seq(dim(df_valid_assignment)[1])){
  
  if ((id %% 100)==0){
    print(id)
  }
  
  cellBC_oi <- df_valid_assignment$cellBC[id]
  clone_oi <- df_valid_assignment$clone_id[id]
  clone_BC1_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(BC1)
  clone_BC2_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(BC2)
  
  df_clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi)
  
  df_all_capture_confident_clones <- df_T7BC_valid2_1umi %>% filter(cBC_pdT %in% cellBC_oi) %>% 
    transform(in_best_clone_BC1= BC1 %in% clone_BC1_oi,
              in_best_clone_BC2= BC2 %in% clone_BC2_oi)
  # tic()
  precision_oi_BC1 <- c()
  precision_oi_BC2 <- c()
  recall_oi_BC1 <- c()
  recall_oi_BC2 <- c()
  for (thres in umi_ct_thresh){
    precision_oi_BC1[thres] <- mean(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres) %>% pull(in_best_clone_BC1) %>% mean())
    precision_oi_BC2[thres] <- mean(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres) %>% pull(in_best_clone_BC2) %>% mean())
    
    recall_oi_BC1[thres] <- mean(!is.na(df_clone_BCs_oi %>% 
                                      left_join(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres),by = join_by(BC1)) %>% pull(UMIs_BC)))
    recall_oi_BC2[thres] <- mean(!is.na(df_clone_BCs_oi %>% 
                                          left_join(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres),by = join_by(BC2)) %>% pull(UMIs_BC)))
  }
  
  df_PR_oi <- data.frame(cellBC=cellBC_oi,
                         high_conf_clone=clone_oi,
                         umi_thresh=umi_ct_thresh,
                         precision_BC1=precision_oi_BC1,
                         precision_BC2=precision_oi_BC2,
                         recall_BC1=recall_oi_BC1,
                         recall_BC2=recall_oi_BC2)
  
  
  df_precision_recall_high_conf_clones <- rbind(df_precision_recall_high_conf_clones,df_PR_oi)

}


# 
# 
# write.table(df_precision_recall_high_conf_clones,
#             "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/table_PR_high_confidence_round2_assignments_v2_precthresh0p75_sep_BC_20240112.txt",
#             row.names=FALSE, quote = FALSE)


write.table(df_precision_recall_high_conf_clones,
            "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/table_PR_high_confidence_round2_assignments_v2_precthresh0p5_sep_BC_20240112.txt",
            row.names=FALSE, quote = FALSE)

# df_precision_recall_high_conf_clones <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/table_PR_high_confidence_round2_assignments_v2_precthresh0p75_20240107.txt",
#            header=TRUE)
# 
# write.table(df_precision_recall_high_conf_clones,
#             "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/table_PR_high_confidence_round2_assignments_v2_precthresh0p5_20240107.txt",
#             row.names=FALSE, quote = FALSE)


df_av_PR_all <- df_precision_recall_high_conf_clones %>% group_by(umi_thresh) %>% 
  summarize(median_prec=median(precision_BC1,na.rm=TRUE),
            mean_prec=mean(precision_BC1,na.rm=TRUE),
            median_recall=median(recall_BC1,na.rm=TRUE),
            mean_recall=mean(recall_BC1,na.rm=TRUE))

df_av_PR_all

dev.set(4)
plt_PR <- ggplot(df_av_PR_all) + 
  geom_line(aes(x=umi_thresh,y=mean_prec),color="firebrick",linewidth=0.3)+
  geom_line(aes(x=umi_thresh,y=mean_recall),color="dodgerblue",linewidth=0.3)+
  geom_point(aes(x=umi_thresh,y=mean_prec),color="firebrick",size=0.3)+
  geom_point(aes(x=umi_thresh,y=mean_recall),color="dodgerblue",size=0.3)+
  # coord_fixed()+
  theme_cowplot()+
  # scale_x_continuous(limits=c(0,1.01),expand=expansion(0,0.05))+
  scale_y_continuous(limits=c(0,1.0),expand=expansion(0,0))+
  # facet_wrap(~high_conf_clone,nrow=5)+ 
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=9),
    panel.grid.major=element_line(size=0.1,color="grey"),
    axis.ticks = element_line(size=0.2),
    axis.text=element_text(size=9),
    axis.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.title=element_text(size=9))+
  labs(x="UMI threshold",y="Precision (red) or Recall (blue)")

dev.set(4)

dev.new()
plt_PR

pdf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/PR_assigned_clone_vs_UMI_threshold_sep_BC_20240113.pdf",width=3.5,height=3)
print(plt_PR)
dev.off()
