

library(data.table)


# load consensus clonotype BC table:
df_BC_clonotypes <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/denovo3umi_w_bulk_joined_round1_clonotypes_BC_pairs_20231220.txt",
                               header=TRUE)

# df_BC_clonotypes %>% group_by(consensus_clone_id) %>% summarize(MOI=length(real_barcode1)) %>% pull(MOI) %>% quantile()
df_BC_clonotypes <- df_BC_clonotypes %>% transform(
  BC_pair=str_replace(BC_pair,"-","_"),
  real_barcode1=str_split(BC_pair,"-") %>% lapply("[[",1) %>% unlist(),
  real_barcode2=str_split(BC_pair,"-") %>% lapply("[[",2) %>% unlist())

list_clones <- df_BC_clonotypes %>% pull(consensus_clone_id) %>% unique()

# # convert to long version: 
# df_BC_clonotypes_lg <- df_BC_clonotypes 
# colnames(df_BC_clonotypes_lg) <- c("consensus_clone_id","CS1","CS2")
# df_BC_clonotypes_lg <- df_BC_clonotypes_lg %>% pivot_longer(cols=c(CS1,CS2))
# list_clones <- df_BC_clonotypes_lg %>% pull(consensus_clone_id) %>% unique()


# 10x T7 data from 'valid' cellBC: 
df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/valid_T7BC_Cre_v3_CS_20231220.txt",
                            header=TRUE)


# considering BC separately:
# BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split("_")
# df_T7BC_valid$BC1 <- BC_all %>% lapply("[[",1) %>% unlist()
# df_T7BC_valid$BC2 <- BC_all %>% lapply("[[",2) %>% unlist()
# df_T7BC_valid$capture_direction <- BC_all %>% lapply("[[",3) %>% unlist()
# df_T7BC_valid <- df_T7BC_valid %>% select(-BC)
# df_T7BC_valid2 <- df_T7BC_valid %>% group_by(cBC_pdT,BC1,BC2) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC))


BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split("_CS")
df_T7BC_valid$BC_pair <- BC_all %>% lapply("[[",1) %>% unlist()
df_T7BC_valid$capture_direction <- paste0("CS",BC_all %>% lapply("[[",2) %>% unlist())
df_T7BC_valid <- df_T7BC_valid %>% select(-BC)


df_T7BC_valid2_1umi <- df_T7BC_valid %>% group_by(cBC_pdT,BC_pair) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC))
df_T7BC_valid2 <- df_T7BC_valid2_1umi %>% filter(UMIs_BC>=2)


# df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/valid_T7BC_Cre_v2_pdT_20231210.txt"
#                             header=TRUE)




cellBC <- df_T7BC_valid2 %>% pull(cBC_pdT) %>% unique()


list_df_mapping_to_clonotype <- list()

counter <- 0
# tic()
for (cellBC_oi in cellBC){
  
  if ((counter%%100)==0){
    print(counter)
    # if (counter==1000){
    #   break
    # }
  }
  
  # # keeping BC separate
  # T7BC1_oi <- df_T7BC_valid %>% filter(cBC_pdT %in% cellBC_oi) %>% pull(BC1)
  # T7BC2_oi <- df_T7BC_valid %>% filter(cBC_pdT %in% cellBC_oi) %>% pull(BC2)
  # df_frac_from_clones_in_10xBC1 <- df_BC_clonotypes %>% group_by(consensus_clone_id) %>% 
  #   summarize(frac_clone_in_10x1=mean(real_barcode1 %in% T7BC1_oi)) %>% data.frame()
  # df_frac_from_clones_in_10xBC2 <- df_BC_clonotypes %>% group_by(consensus_clone_id) %>% 
  #   summarize(frac_clone_in_10x2=mean(real_barcode2 %in% T7BC2_oi)) %>% data.frame()
  #
  # df_frac_from_10xBC1_in_clones <- data.frame()
  # df_frac_from_10xBC2_in_clones <- data.frame()
  # 
  # for (clone_oi in list_clones){
  #   clone_BC1s <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(real_barcode1)
  #   clone_BC2s <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(real_barcode2)
  #   
  #   df_frac_from_10xBC1_in_clones <- rbind(df_frac_from_10xBC1_in_clones,
  #                                         data.frame(consensus_clone_id=clone_oi,
  #                                                    frac_10xBC1_in_clone=mean(T7BC1_oi %in% clone_BC1s)))
  #   
  #   df_frac_from_10xBC2_in_clones <- rbind(df_frac_from_10xBC2_in_clones,
  #                                          data.frame(consensus_clone_id=clone_oi,
  #                                                     frac_10xBC2_in_clone=mean(T7BC2_oi %in% clone_BC2s)))
  #   
  # }
  # 
  # # combine all information at this stage: 
  # 
  # df_mapping_to_clonotype_oi <- data.frame(cellBC=cellBC_oi,
  #                                       df_frac_from_clones_in_10xBC1) %>% 
  #   left_join(df_frac_from_clones_in_10xBC2,by="consensus_clone_id") %>% 
  #   left_join(df_frac_from_10xBC1_in_clones,by="consensus_clone_id") %>% 
  #   left_join(df_frac_from_10xBC2_in_clones,by="consensus_clone_id")

  
  T7BC_oi <- df_T7BC_valid2 %>% filter(cBC_pdT %in% cellBC_oi) %>% pull(BC_pair)

  df_frac_from_clones_in_10x <- df_BC_clonotypes %>% group_by(consensus_clone_id) %>%
    summarize(frac_clone_in_10x=mean(BC_pair %in% T7BC_oi)) %>% data.frame()
  

  df_frac_from_10x_in_clones <- data.frame()

  for (clone_oi in list_clones){
    clone_BCs <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(BC_pair)

    df_frac_from_10x_in_clones <- rbind(df_frac_from_10x_in_clones,
                                          data.frame(consensus_clone_id=clone_oi,
                                                     frac_10x_in_clone=mean(T7BC_oi %in% clone_BCs)))
  }

  # combine all information at this stage:

  df_mapping_to_clonotype_oi <- data.frame(cellBC=cellBC_oi,
                                           df_frac_from_clones_in_10x) %>%
    left_join(df_frac_from_10x_in_clones,by="consensus_clone_id")
  
  
  list_df_mapping_to_clonotype[[cellBC_oi]] <- df_mapping_to_clonotype_oi

  
  counter <- counter+1
}


library(data.table)
df_mapping_to_clonotype <- rbindlist(list_df_mapping_to_clonotype)

write.table(df_mapping_to_clonotype,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_v2/mapping_Cre_T7BC_CS_denovo3umi_w_bulk_20231220.txt",
            quote=FALSE,row.names=FALSE,sep="\t")

# both
# df_mapping_to_clonotype_best1 <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_clone_in_10x1,n=1,with_ties=FALSE)
# df_mapping_to_clonotype_best2 <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_clone_in_10x2,n=1,with_ties=FALSE)

# dev.new()
# dev.set(4)
# pltBC1 <- ggplot(df_mapping_to_clonotype_best1 %>% filter(frac_clone_in_10x1>0 | frac_10xBC1_in_clone>0)) +
#   geom_bin_2d(aes(x=frac_clone_in_10x1,y=frac_10xBC1_in_clone))+xlim(c(NA,1))+coord_fixed()
# 
# 
# pltBC2 <- ggplot(df_mapping_to_clonotype_best2 %>% filter(frac_clone_in_10x2>0 | frac_10xBC2_in_clone>0)) +
#   geom_bin_2d(aes(x=frac_clone_in_10x2,y=frac_10xBC2_in_clone))+xlim(c(NA,1))+coord_fixed()
# 
# pltBC1+pltBC2


# best recall (capture)
df_mapping_to_clonotype_best <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_clone_in_10x,n=1,with_ties=FALSE)
df_mapping_to_clonotype_best2 <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_clone_in_10x,n=2,with_ties=FALSE)

# best precision (purity)
df_mapping_to_clonotype_best <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_10x_in_clone,n=1,with_ties=FALSE)
df_mapping_to_clonotype_best2 <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_10x_in_clone,n=2,with_ties=FALSE)


df_mapping_to_clonotype_best2_only <- anti_join(df_mapping_to_clonotype_best2,df_mapping_to_clonotype_best)
colnames(df_mapping_to_clonotype_best2_only)[c(2,3,4)] <- c("consensus_clone_id_2ndBest","frac_clone_in_10x_2ndBest","frac_10x_in_clone_2ndBest")
df_mapping_to_clonotype_best_w_2nd_best <- df_mapping_to_clonotype_best %>% left_join(df_mapping_to_clonotype_best2_only)
colnames(df_mapping_to_clonotype_best_w_2nd_best)[c(2,3,4)] <- c("clone_id","recall","precision")
colnames(df_mapping_to_clonotype_best_w_2nd_best)[c(5,6,7)] <- c("clone_id_2","recall_2","precision_2")


# df_mapping_to_clonotype_best <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_clone_in_10x,n=1,with_ties=FALSE)

# df_mapping_to_clonotype_best_precision <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_10x_in_clone,n=1,with_ties=FALSE)


df_mapping_to_clonotype_best %>% filter(!(frac_clone_in_10x>0 | frac_10x_in_clone>0)) %>% dim()
df_mapping_to_clonotype_best %>% filter(!(frac_clone_in_10x>0 | frac_10x_in_clone>0)) %>% dim()



thresh_capt_low <- 0.1
thresh_capt_hi <- 0.3
thresh_precision <- 0.75

dev.set(4)
ggplot(df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x>0 | frac_10x_in_clone>0)) +
  geom_bin_2d(aes(x=frac_clone_in_10x,y=frac_10x_in_clone))+xlim(c(NA,1))+coord_fixed()+
  geom_vline(xintercept=thresh_capt_low,color="orange",linetype="dashed")+
  geom_vline(xintercept=thresh_capt_hi,color="red",linetype="dashed")+
  geom_hline(yintercept=thresh_precision,color="green",linetype="dashed")+
  scale_x_continuous(expand=expansion(0,0))+
  scale_y_continuous(expand=expansion(0,0))+
  labs(x="Recall (fraction barcodes from best\n clone found in 10x data)",
       y="Precision (fraction barcodes from\n 10x from best clone)")+theme_cowplot()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=9),
    panel.grid.major=element_line(size=0.1,color="grey"),
    axis.ticks = element_line(size=0.2),
    axis.text=element_text(size=9),
    axis.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.title=element_text(size=9))


# categories
df_cat_no_capt <-  df_mapping_to_clonotype_best %>% filter(!(frac_clone_in_10x>0 | frac_10x_in_clone>0)) %>% 
  transform(assignment_category="no_capture")

df_cat_lo_capt <- df_mapping_to_clonotype_best %>% filter((frac_clone_in_10x>0 | frac_10x_in_clone>0) & frac_clone_in_10x<thresh_capt_low) %>% 
  transform(assignment_category="low_capture")

# df_mapping_to_clonotype_best %>% filter((frac_clone_in_10x>0 | frac_10x_in_clone>0) & frac_clone_in_10x<thresh_capt_low) %>% dim
# df_mapping_to_clonotype_best %>% filter((frac_clone_in_10x>0 | frac_10x_in_clone>0) & frac_clone_in_10x<thresh_capt_low & frac_10x_in_clone>=thresh_precision) %>% dim
# df_mapping_to_clonotype_best %>% filter((frac_clone_in_10x>0 | frac_10x_in_clone>0) & frac_clone_in_10x<thresh_capt_low & frac_10x_in_clone<thresh_precision) %>% dim


df_cat_mid_capt_hi_purity <- df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x>=thresh_capt_low & frac_clone_in_10x<thresh_capt_hi & frac_10x_in_clone>=thresh_precision)%>% 
  transform(assignment_category="mid_capture_high_purity")

df_cat_hi_capt_hi_purity <- df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x>=thresh_capt_hi & frac_10x_in_clone>=thresh_precision) %>% 
  transform(assignment_category="high_capture_high_purity")

df_cat_lo_purity <- df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x>=thresh_capt_low & frac_10x_in_clone<thresh_precision) %>% 
  transform(assignment_category="low_purity")


df_mapping_to_clonotype_w_category <- rbind(df_cat_no_capt,
                                            df_cat_lo_capt,
                                            df_cat_mid_capt_hi_purity,
                                            df_cat_hi_capt_hi_purity,
                                            df_cat_lo_purity)
table(df_mapping_to_clonotype_w_category$assignment_category)



df_mapping_clonotype_final_round1 <- df_mapping_to_clonotype_best_w_2nd_best %>% left_join(df_mapping_to_clonotype_w_category %>% select(cellBC,assignment_category))

df_mapping_clonotype_final_round1_recall <- df_mapping_clonotype_final_round1
df_mapping_clonotype_final_round1_precision <- df_mapping_clonotype_final_round1


df_comp_selection <- df_mapping_clonotype_final_round1_precision %>% select(cellBC,clone_id_PRE=clone_id,assignment_category_PRE=assignment_category) %>% left_join(
  df_mapping_clonotype_final_round1_recall %>% select(cellBC,clone_id_REC=clone_id,assignment_category_REC=assignment_category)
)

df_comp_selection %>% group_by(assignment_category_PRE) %>% summarize(frac_congruent=mean(clone_id_PRE==clone_id_REC))
df_comp_selection %>% group_by(assignment_category_REC) %>% summarize(frac_congruent=mean(clone_id_PRE==clone_id_REC))



df_mapping_clonotype_final_round1 %>% filter(recall_2>=thresh_capt_hi & precision_2>=thresh_precision)%>% data.frame()
# df_mapping_clonotype_final_round1 %>% filter(recall_2>=thresh_capt_hi & precision_2>=thresh_precision)%>% data.frame()


df_mapping_clonotype_final_round1_2 <- df_mapping_clonotype_final_round1_precision #%>% transform(ifelse(assignment_category=="high_capture_high_purity"))
df_mapping_clonotype_final_round1_2$assignment_category[df_mapping_clonotype_final_round1_2$assignment_category=="high_capture_high_purity" &
                                                          df_mapping_clonotype_final_round1_2$recall_2>=0.1] <- "high_capture_high_purity_too_high_recall_2nd"
df_mapping_clonotype_final_round1_2$assignment_category[df_mapping_clonotype_final_round1_2$assignment_category=="mid_capture_high_purity" &
                                                          df_mapping_clonotype_final_round1_2$recall_2>=0.1] <- "mid_capture_high_purity_too_high_recall_2nd"

table(df_mapping_clonotype_final_round1_2$assignment_category)

df_mapping_clonotype_final_round1_3 <- df_mapping_clonotype_final_round1_2 %>% 
  transform(valid_assignment=ifelse(assignment_category %in% c("high_capture_high_purity","mid_capture_high_purity"),TRUE,FALSE))

table(df_mapping_clonotype_final_round1_3$valid_assignment)


write.table(df_mapping_clonotype_final_round1_3,
            "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_v2/round1_assignment_cells_to_denovo3umi_bulk_clonotypes_20231221.txt",
            row.names=FALSE,sep="\t",quote=FALSE)

dev.set(4)
ggplot(df_mapping_clonotype_final_round1_3 %>% filter(valid_assignment)) + 
  geom_bin_2d(aes(x=recall,y=precision))+
  facet_wrap(~clone_id)

df_summary_per_clone <- df_mapping_clonotype_final_round1_3 %>% filter(valid_assignment) %>% group_by(clone_id) %>% 
  summarize(
    mean_recall=mean(recall),
    mean_precision=mean(precision),
    n_cells=length(cellBC))

plt1 <- ggplot(df_summary_per_clone)+stat_bin(aes(x=mean_recall),geom="step")
plt2 <- ggplot(df_summary_per_clone)+stat_bin(aes(x=mean_precision),geom="step")
plt3 <- ggplot(df_summary_per_clone)+stat_bin(aes(x=n_cells),geom="step")
plt4 <- ggplot(df_summary_per_clone) + geom_point(aes(y=mean_precision,x=mean_recall,size=n_cells),shape=1)
plt1+plt2+plt3+plt4



# 
# dev.set(4)
# plt1 <- ggplot(df_mapping_clonotype_final_round1) + 
#   stat_bin2d(aes(x=precision,y=precision_2))+
#   facet_wrap(~assignment_category)
# plt2 <- ggplot(df_mapping_clonotype_final_round1) + 
#   stat_bin2d(aes(x=recall,y=recall_2))+
#   facet_wrap(~assignment_category)
# plt3 <- ggplot(df_mapping_clonotype_final_round1) + 
#   stat_bin2d(aes(y=precision_2,x=recall_2))+
#   facet_wrap(~assignment_category)
# plt1+plt2+plt3+plot_layout(ncol=1)


# df_sum_T7_umi <- df_T7BC_valid2 %>% group_by(cBC_pdT) %>% summarize(total_T7_UMI=sum(UMIs_BC))
# colnames(df_sum_T7_umi)[1] <- "cellBC" 
# 
# df_mapping_to_clonotype_w_category2 <- df_mapping_to_clonotype_w_category %>% left_join(df_sum_T7_umi)
# 
# dev.set(4)
# ggplot(df_mapping_to_clonotype_w_category2) + stat_bin(aes(x=total_T7_UMI),geom="step")+
#   scale_x_log10()+facet_wrap(~assignment_category)
# ggplot(df_mapping_to_clonotype_w_category2) + 
#   geom_step(aes(x=total_T7_UMI,color=assignment_category),stat="ecdf")+
#   scale_x_log10(limits=c(1,1000))

  #%>% filter(frac_clone_in_10x>0 | frac_10x_in_clone>0)

# ggplot(df_mapping_to_clonotype_best_precision %>% filter(frac_clone_in_10x>0 | frac_10x_in_clone>0)) +
#   geom_bin_2d(aes(x=frac_clone_in_10x,y=frac_10x_in_clone,fill=sqrt(..count..)))+xlim(c(NA,1))+coord_fixed()


# # fetch all the data from the high confidence >=2 umi/cell 
# df_mapping_to_clonotype_best_valid <- df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x>0.25 & frac_10x_in_clone>0.75)
# write.table(df_mapping_to_clonotype_best_valid,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/high_confidence_assignment_to_bulk_derived_clonotypes_20231214.txt",
#             sep="\t",quote=FALSE,row.names=FALSE)


# df_mapping_to_clonotype_best_valid_loose <- df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x>0.125 & frac_10x_in_clone>0.25)
# df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x > 0.25 & frac_10x_in_clone>0.75) %>% dim()
# df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x > 0.05 & frac_clone_in_10x <= 0.25 & frac_10x_in_clone>0.75) %>% dim()
# df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x > 0.25 & frac_10x_in_clone<=0.75) %>% dim()
# df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x > 0.05 & frac_clone_in_10x <= 0.25 &
#                                           frac_10x_in_clone<=0.75) %>% dim()




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
  clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(BC_pair)
  
  df_clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi)

  df_all_capture_confident_clones <- df_T7BC_valid2_1umi %>% filter(cBC_pdT %in% cellBC_oi) %>% 
    transform(in_best_clone=BC_pair %in% clone_BCs_oi)
  tic()
  precision_oi <- c()
  recall_oi <- c()
  for (thres in umi_ct_thresh){
      precision_oi[thres] <- mean(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres) %>% pull() %>% mean())
      recall_oi[thres] <- mean(!is.na(df_clone_BCs_oi %>% 
      left_join(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres),by = join_by(BC_pair)) %>% pull(UMIs_BC)))
  }
  
  df_PR_oi <- data.frame(cellBC=cellBC_oi,
                         high_conf_clone=clone_oi,
                         umi_thresh=umi_ct_thresh,
                         precision=precision_oi,
                         recall=recall_oi)

  
  df_precision_recall_high_conf_clones <- rbind(df_precision_recall_high_conf_clones,df_PR_oi)

}




write.table(df_precision_recall_high_conf_clones,
            "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_v2/table_PR_high_confidence_round1_assignments_20231221.txt",
            row.names=FALSE, quote = FALSE)


df_av_PR_all <- df_precision_recall_high_conf_clones %>% group_by(umi_thresh) %>% 
  summarize(median_prec=median(precision,na.rm=TRUE),
            mean_prec=mean(precision,na.rm=TRUE),
            median_recall=median(recall,na.rm=TRUE),
            mean_recall=mean(recall,na.rm=TRUE))



dev.set(5)
ggplot(df_av_PR_all) + 
  geom_line(aes(x=umi_thresh,y=mean_prec),color="firebrick",linewidth=0.3)+
  geom_line(aes(x=umi_thresh,y=mean_recall),color="dodgerblue",linewidth=0.3)+
  geom_point(aes(x=umi_thresh,y=mean_prec),color="firebrick",size=0.3)+
  geom_point(aes(x=umi_thresh,y=mean_recall),color="dodgerblue",size=0.3)+
  # coord_fixed()+
  theme_cowplot()+
  # scale_x_continuous(limits=c(0,1.01),expand=expansion(0,0.05))+
  scale_y_continuous(limits=c(0,1.01),expand=expansion(0,0))+
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


# ggplot()+geom_point(aes(x=precision_oi,y=recall_oi))+
#   coord_fixed()+
#   xlim(c(0,1.05))+ylim(c(0,1.05))+theme_cowplot()+
#   # facet_wrap(~umi_thresh,nrow=2)+
#   theme(
#     strip.background=element_blank(),
#     strip.text=element_text(size=9),
#     panel.grid.major=element_line(size=0.1,color="grey"),
#     axis.ticks = element_line(size=0.2),
#     axis.text=element_text(size=9),
#     axis.title=element_text(size=10),
#     legend.text=element_text(size=7),
#     legend.title=element_text(size=9))


# ggplot(df_precision_recall_high_conf_clones)+geom_bin_2d(aes(x=precision,y=recall))+
#   coord_fixed()+
#   xlim(c(0,1.05))+ylim(c(0,1.05))+theme_cowplot()+
#   facet_wrap(~umi_thresh,nrow=2)+
#   theme(
#     strip.background=element_blank(),
#     strip.text=element_text(size=9),
#     panel.grid.major=element_line(size=0.1,color="grey"),
#     axis.ticks = element_line(size=0.2),
#     axis.text=element_text(size=9),
#     axis.title=element_text(size=10),
#     legend.text=element_text(size=7),
#     legend.title=element_text(size=9))


















umi_ct_thresh <- seq(1,10)

df_precision_recall_mid_conf_clones <- data.frame()
for (id in seq(dim(df_mapping_to_clonotype_best_valid_loose)[1])){
  
  if ((id %% 100)==0){
    print(id)
  }
  
  cellBC_oi <- df_mapping_to_clonotype_best_valid_loose$cellBC[id]
  clone_oi <- df_mapping_to_clonotype_best_valid_loose$consensus_clone_id[id]
  clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(BC_pair)
  
  df_clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi)
  
  df_all_capture_confident_clones <- df_T7BC_valid2_1umi %>% filter(cBC_pdT %in% cellBC_oi) %>% 
    transform(in_best_clone=BC_pair %in% clone_BCs_oi)
  tic()
  precision_oi <- c()
  recall_oi <- c()
  for (thres in umi_ct_thresh){
    precision_oi[thres] <- mean(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres) %>% pull() %>% mean())
    recall_oi[thres] <- mean(!is.na(df_clone_BCs_oi %>% 
                                      left_join(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres),by = join_by(BC_pair)) %>% pull(UMIs_BC)))
  }
  
  df_PR_oi <- data.frame(cellBC=cellBC_oi,
                         high_conf_clone=clone_oi,
                         umi_thresh=umi_ct_thresh,
                         precision=precision_oi,
                         recall=recall_oi)
  
  
  df_precision_recall_mid_conf_clones <- rbind(df_precision_recall_mid_conf_clones,df_PR_oi)
  
}












hi_conf_cBC <- df_mapping_to_clonotype_best_valid %>% pull(cellBC)
mid_conf_cBC <- df_mapping_to_clonotype_best_valid_loose %>% pull(cellBC)
mid_conf_cBC <- mid_conf_cBC[!(mid_conf_cBC %in% hi_conf_cBC)]


df_intermediate <- df_mapping_to_clonotype_best_valid_loose %>% filter(cellBC %in% mid_conf_cBC)

# id <- 1

# cellBC_oi <- df_intermediate$cellBC[id]
# clone_oi <- df_intermediate$consensus_clone_id[id]
# clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(BC_pair)
# df_clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi)
# df_all_capture_confident_clones <- df_T7BC_valid2_1umi %>% filter(cBC_pdT %in% cellBC_oi) %>% 
#   transform(in_best_clone=BC_pair %in% clone_BCs_oi)


for (id in seq(50)){ #dim(df_mapping_to_clonotype_best_valid)[1])
  
  cellBC_oi <- df_intermediate$cellBC[id]
  clone_oi <- df_intermediate$consensus_clone_id[id]
  clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(BC_pair)
  
  df_example <- df_T7BC_valid2_1umi %>% filter(cBC_pdT %in% cellBC_oi) %>% 
    transform(in_best_clone=BC_pair %in% clone_BCs_oi)
  
  plt <- ggplot(df_example) + geom_point(aes(x=reads_BC,y=UMIs_BC,color=in_best_clone))+
    scale_x_log10(limits=c(1,3000))+scale_y_log10(limits=c(0.5,40))+coord_fixed()+
    labs(title=sprintf("cell: %s\nclone: %s, capt=%.2f, purity=%.2f",
                       cellBC_oi,clone_oi,
                       df_intermediate$frac_clone_in_10x[id],
                       df_intermediate$frac_10x_in_clone[id]),
         size=8)+
    theme(title=element_text(size=7),
          legend.text=element_text(size=6))
  
  # pdf(sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/example_mapped_clones_CS/cell_%d_20231210.pdf",id),
  #     width=4,height=3)
  pdf(sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/example_intermediate_quality_clones/cell_%d_20231215.pdf",id),
      width=4,height=3)
  print(plt)
  dev.off()
}








df_av_PR_all <- df_precision_recall_high_conf_clones %>% group_by(umi_thresh) %>% 
  summarize(median_prec=median(precision,na.rm=TRUE),
            mean_prec=mean(precision,na.rm=TRUE),
            median_recall=median(recall,na.rm=TRUE),
            mean_recall=mean(recall,na.rm=TRUE))

df_av_PR_all <- df_precision_recall_mid_conf_clones %>% filter(cellBC %in% mid_conf_cBC) %>% group_by(umi_thresh) %>% 
  summarize(median_prec=median(precision,na.rm=TRUE),
            mean_prec=mean(precision,na.rm=TRUE),
            median_recall=median(recall,na.rm=TRUE),
            mean_recall=mean(recall,na.rm=TRUE))



dev.set(5)
ggplot(df_av_PR_all) + 
  geom_line(aes(x=umi_thresh,y=mean_prec),color="firebrick",linewidth=0.3)+
  geom_line(aes(x=umi_thresh,y=mean_recall),color="dodgerblue",linewidth=0.3)+
  geom_point(aes(x=umi_thresh,y=mean_prec),color="firebrick",size=0.3)+
  geom_point(aes(x=umi_thresh,y=mean_recall),color="dodgerblue",size=0.3)+
  # coord_fixed()+
  theme_cowplot()+
  # scale_x_continuous(limits=c(0,1.01),expand=expansion(0,0.05))+
  scale_y_continuous(limits=c(0,1.01),expand=expansion(0,0))+
  # facet_wrap(~high_conf_clone,nrow=5)+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=9),
    panel.grid.major=element_line(size=0.1,color="grey"),
    axis.ticks = element_line(size=0.2),
    axis.text=element_text(size=9),
    axis.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.title=element_text(size=9))


df_av_PR <- df_precision_recall_high_conf_clones %>% group_by(high_conf_clone,umi_thresh) %>% 
  summarize(av_prec=median(precision,na.rm=TRUE),
            av_recall=median(recall,na.rm=TRUE))


ggplot(df_av_PR) + geom_point(aes(x=av_prec,y=av_recall))+
  coord_fixed()+
  theme_cowplot()+
  scale_x_continuous(limits=c(0,1.01),expand=expansion(0,0.05))+
  scale_y_continuous(limits=c(0,1.01),expand=expansion(0,0))+
  facet_wrap(~high_conf_clone,nrow=5)+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=9),
    panel.grid.major=element_line(size=0.1,color="grey"),
    axis.ticks = element_line(size=0.2),
    axis.text=element_text(size=9),
    axis.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.title=element_text(size=9))

ggplot(df_av_PR) + 
  geom_line(aes(x=umi_thresh,y=av_prec),color="firebrick",linewidth=0.3)+
  geom_line(aes(x=umi_thresh,y=av_recall),color="dodgerblue",linewidth=0.3)+
  geom_point(aes(x=umi_thresh,y=av_prec),color="firebrick",size=0.3)+
  geom_point(aes(x=umi_thresh,y=av_recall),color="dodgerblue",size=0.3)+
  # coord_fixed()+
  theme_cowplot()+
  # scale_x_continuous(limits=c(0,1.01),expand=expansion(0,0.05))+
  scale_y_continuous(limits=c(0,1.01),expand=expansion(0,0))+
  facet_wrap(~high_conf_clone,nrow=5)+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=9),
    panel.grid.major=element_line(size=0.1,color="grey"),
    axis.ticks = element_line(size=0.2),
    axis.text=element_text(size=9),
    axis.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.title=element_text(size=9))

setwd("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/PR_consensus_clonotypes")

# for (cid in df_av_PR$high_conf_clone %>% unique){
  
  # df_av_PR_oi <- df_av_PR %>% 
  # 
# }


dev.set(5)
ggplot()+geom_point(aes(x=precision_oi,y=recall_oi))+coord_fixed()+
  xlim(c(0,1))+ylim(c(0,1))+theme_cowplot()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=9),
    panel.grid.major=element_line(size=0.1,color="grey"),
    axis.ticks = element_line(size=0.2),
    axis.text=element_text(size=9),
    axis.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.title=element_text(size=9))

df_clone_quant <- df_mapping_to_clonotype_best_valid %>% group_by(consensus_clone_id) %>% summarize(mean_capture=mean(frac_clone_in_10x),
                                                                                                    mean_purity=mean(frac_10x_in_clone)) %>% 
  left_join(df_BC_clonotypes %>% group_by(consensus_clone_id) %>% summarize(MOI=length(real_barcode1)))

df_clone_quant %>% arrange(desc(mean_capture))


dev.set(5)
ggplot(df_clone_quant) + geom_point(aes(x=mean_capture,y=mean_purity,color=MOI))+
  xlim(0,1)+ylim(0,1)+coord_fixed()

df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x<0.25 & frac_10x_in_clone<0.25) %>% dim()


df_mapping_to_clonotype_best_valid <- df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x>0.25 | frac_10x_in_clone>0.25)

for (id in seq(50)){ #dim(df_mapping_to_clonotype_best_valid)[1])
  
  cellBC_oi <- df_mapping_to_clonotype_best_valid$cellBC[id]
  clone_oi <- df_mapping_to_clonotype_best_valid$consensus_clone_id[id]
  clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(BC_pair)
  
  df_example <- df_T7BC_valid2 %>% filter(cBC_pdT %in% cellBC_oi) %>% 
    transform(in_best_clone=BC_pair %in% clone_BCs_oi)
  
  plt <- ggplot(df_example) + geom_point(aes(x=reads_BC,y=UMIs_BC,color=in_best_clone))+
    scale_x_log10(limits=c(1,3000))+scale_y_log10(limits=c(0.5,40))+coord_fixed()+
    labs(title=sprintf("cell: %s\nclone: %s, capt=%.2f, purity=%.2f",
                       cellBC_oi,clone_oi,
                       df_mapping_to_clonotype_best_valid$frac_clone_in_10x[id],
                       df_mapping_to_clonotype_best_valid$frac_10x_in_clone[id]),
         size=8)+
    theme(title=element_text(size=7),
          legend.text=element_text(size=6))
  
  # pdf(sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/example_mapped_clones_CS/cell_%d_20231210.pdf",id),
  #     width=4,height=3)
  pdf(sprintf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/example_mapped_clones_pdT/cell_%d_20231210.pdf",id),
      width=4,height=3)
  print(plt)
  dev.off()
}
plt



dev.set(4)









df_mapping_to_clonotype_best1 %>% filter(cellBC=="AAACCCAAGAAATGGG-1")















plt2 <- ggplot(df_mapping_to_clonotype_best ) + stat_bin(aes(x=frac_clone_in_10x))+
  scale_y_log10()
plt3 <- ggplot(df_mapping_to_clonotype_best ) + stat_bin(aes(y=frac_10x_in_clone))+
  scale_x_log10()


plt_design <- "B#
AC"
plt1+plt2+plt3+plot_layout(design=plt_design)

# %>% filter(frac_clone_in_10x>0 | frac_10x_in_clone>0)



df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x>0.1 & frac_10x_in_clone>0.75) %>% dim()

df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x<0.1 & frac_10x_in_clone<0.1) %>% dim()







df_mapping_to_clonotype_best_1umi <- df_mapping_to_clonotype_best
colnames(df_mapping_to_clonotype_best_1umi) <- c("cellBC","consensus_clone_id_1umi",
                                                 "frac_clone_in_10x_1umi","frac_10x_in_clone_1umi")

df_mapping_to_clonotype_best_2umi <- df_mapping_to_clonotype_best
colnames(df_mapping_to_clonotype_best_2umi) <- c("cellBC","consensus_clone_id_2umi",
                                                 "frac_clone_in_10x_2umi","frac_10x_in_clone_2umi")

df_compare_1_2_umi <- df_mapping_to_clonotype_best_2umi %>% left_join(df_mapping_to_clonotype_best_1umi)

mean(df_compare_1_2_umi$consensus_clone_id_2umi == df_compare_1_2_umi$consensus_clone_id_1umi,na.rm=TRUE)

plt1 <- ggplot(df_compare_1_2_umi %>% filter(consensus_clone_id_2umi==consensus_clone_id_1umi)) + geom_bin_2d(aes(x=frac_clone_in_10x_1umi,y=frac_clone_in_10x_2umi))+coord_fixed()+
  geom_abline()
plt2 <- ggplot(df_compare_1_2_umi) + geom_bin_2d(aes(x=frac_10x_in_clone_1umi,y=frac_10x_in_clone_2umi))+coord_fixed()+
  geom_abline()

plt1+plt2

length(df_mapping_to_clonotype_best %>% pull(cellBC) %>% unique())
# 
# table(df_mapping_to_clonotype_best %>% pull(cellBC)) %>% sort()

dir_clones <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq055SP_T7_shuffle_10X/cloneBCs"
setwd(dir_clones)
files_cloneBC <- dir(path=dir_clones,pattern=glob2rx("*clone*.tsv"))

df_mapping_to_clonotype_best %>% filter(cellBC=="ACGTTCCTCCTACGGG-1")




df_all_clones <- data.frame()
for (file_oi in files_cloneBC){
  
  df_clone <- read.table(file_oi,header=TRUE)
  
  clone_id <- str_replace(file_oi,".tsv","") %>% str_split("_") %>% lapply("[[",4) %>% unlist()
  
  df_clone2 <- df_clone %>% transform(clone_id=clone_id)
  
  df_all_clones <- rbind(df_all_clones,
                         df_clone2 %>% select(clone_id,real_barcode1,real_barcode2) )
  
}

head(df_T7BC_valid)


clones_id <- df_all_clones$clone_id %>% unique()

df_pairs_clones <- CombPairs(clones_id)

for (id in seq(dim(df_pairs_clones)[1])){
  
  
  if ((id %% 100)==0){
    print(id)
  }
  
  BC1_clone1 <- df_all_clones %>% filter(clone_id==df_pairs_clones$X1[id]) %>% pull(real_barcode1)
  BC2_clone1 <- df_all_clones %>% filter(clone_id==df_pairs_clones$X1[id]) %>% pull(real_barcode2)
  BC1_clone2 <- df_all_clones %>% filter(clone_id==df_pairs_clones$X2[id]) %>% pull(real_barcode1)
  BC2_clone2 <- df_all_clones %>% filter(clone_id==df_pairs_clones$X2[id]) %>% pull(real_barcode2)
  
  jaccard_pair_BC1 <- length(intersect(BC1_clone1,BC1_clone2))/length(union(BC1_clone1,BC1_clone2))
  jaccard_pair_BC2 <- length(intersect(BC2_clone1,BC2_clone2))/length(union(BC2_clone1,BC2_clone2))
  
  df_pairs_clones$jaccard_pair_BC1[id] <- jaccard_pair_BC1 
  df_pairs_clones$jaccard_pair_BC2[id] <- jaccard_pair_BC2 
  
}

dev.set(4)
ggplot(df_pairs_clones) + stat_bin(aes(x=jaccard_pair_BC1),color="firebrick",geom="step")+
  stat_bin(aes(x=jaccard_pair_BC2),color="blue",geom="step")+
  scale_y_log10()

df_pairs_clones %>% filter(jaccard_pair_BC1>0.125 | jaccard_pair_BC2>0.125)

# Convert the data frame to a matrix
jaccard_matrix_BC1 <- df_pairs_clones %>% select(-jaccard_pair_BC2) %>%
  spread(key = X1, value = jaccard_pair_BC1) %>% column_to_rownames(var = "X2") %>% as.matrix()

jac_mat_BC1_v2 <- matrix(data=NA,nrow=length(clones_id),ncol=length(clones_id))
rownames(jac_mat_BC1_v2) <- clones_id
colnames(jac_mat_BC1_v2) <- clones_id
for (id in seq(dim(df_pairs_clones)[1])){
  jac_mat_BC1_v2[df_pairs_clones$X1[id],df_pairs_clones$X2[id]] <- df_pairs_clones$jaccard_pair_BC1[id]
}


# col_names_mat <- colnames(jaccard_matrix_BC1)
# row_names_mat <- rownames(jaccard_matrix_BC1)
# 
# intersect(col_names_mat,row_names_mat)
# union(col_names_mat,row_names_mat)
# 
# sort.int(row_names_mat,index.return = TRUE)
# sort.int(col_names_mat,index.return = TRUE)


dist_thresh <- 0.125
conn_mat <- jac_mat_BC1_v2>=dist_thresh
g  <- graph.adjacency(conn_mat)
clu <- components(g)

graph <- graph_from_adjacency_matrix(conn_mat, mode = 'undirected')
dev.set(4)
plot(graph,vertex.size=2)





# update matrix to remove doublets: 
list_doublet_clones <- c("B3","C4","C7","C10","H4")
jac_mat_BC1_v2_no_doub <- jac_mat_BC1_v2[!(rownames(jac_mat_BC1_v2) %in% list_doublet_clones),
                                         !(colnames(jac_mat_BC1_v2) %in% list_doublet_clones)]


dist_thresh <- 0.125
conn_mat_no_doub <- jac_mat_BC1_v2_no_doub>=dist_thresh
g_no_doub  <- graph.adjacency(conn_mat_no_doub)
clu_no_doub <- components(g_no_doub)
graph_no_doub <- graph_from_adjacency_matrix(conn_mat_no_doub, mode = 'undirected')
dev.set(4)
plot(graph_no_doub,vertex.size=2)

df_pairs_clones_no_doub <- df_pairs_clones %>% filter(! ((X1 %in% list_doublet_clones) | (X2 %in% list_doublet_clones)))
ggplot(df_pairs_clones_no_doub) + stat_bin(aes(x=jaccard_pair_BC1),color="firebrick",geom="step")+
  stat_bin(aes(x=jaccard_pair_BC2),color="blue",geom="step")+
  scale_y_log10()


df_clone_membership_no_doub <- data.frame(clone_id=names(clu_no_doub$membership),
                                  conn_clu=clu_no_doub$membership)

df_consensus_clone_name <- df_clone_membership_no_doub %>% group_by(conn_clu) %>% summarize(consensus_clone_id=paste0(clone_id,collapse="_"))

df_clone_membership_no_doub2 <- df_clone_membership_no_doub %>% left_join(df_consensus_clone_name)
# go back to origina list of clones:

df_all_clones2 <- df_all_clones %>% filter(!(clone_id %in% list_doublet_clones)) %>% left_join(df_clone_membership_no_doub2)

df_consensus_clones_union <- df_all_clones2 %>% select(consensus_clone_id,real_barcode1,real_barcode2) %>% unique()
table(df_consensus_clones_union$consensus_clone_id)
table(df_all_clones2$clone_id)

write.table(df_consensus_clones_union,"table_consensus_clones_w_BCs_20231203.txt",
      sep="\t",row.names=FALSE,quote=FALSE)

# 
# 
# 
# clones_oi <-df_clone_membership %>% filter(conn_clu==2) %>% pull(clone_id)
# jaccard_matrix_BC1[clones_oi[1],clones_oi[2]]
# jaccard_matrix_BC1[clones_oi[2],clones_oi[1]]
# 
# df_pairs_clones %>% filter( (X1 %in% clones_oi) & (X2 %in% clones_oi))
# 
# 
# df_pairs_clones %>% filter(jaccard_pair_BC1>0.1 & jaccard_pair_BC1<0.8)
# df_pairs_clones %>% filter(jaccard_pair_BC1>0.8)
# 
# 
# # df_pairs_clones2 <- df_pairs_clones %>% rename(X1="clone_id") %>% left_join(df_clone_membership)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

# 
# 
# 
# 
# 
# 
# RC_BCs <- df_T7BC_valid$BC2 %>% DNAStringSet() %>% reverseComplement() %>% as.character()
# mean(df_T7BC_valid$BC2 %in% df_clone2$real_barcode1)
# mean(df_T7BC_valid$BC2 %in% df_clone2$real_barcode2)
# 
# mean(RC_BCs %in% df_clone2$real_barcode1)
# mean(RC_BCs %in% df_clone2$real_barcode2)
# 
# 
# # check unicity/overlap between the clones