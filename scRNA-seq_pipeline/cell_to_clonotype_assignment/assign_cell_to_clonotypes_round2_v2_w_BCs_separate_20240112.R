

library(data.table)
library(cowplot)

# load consensus clonotype BC table:

df_BC_clonotypes <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_round2/FINAL_CLONOTYPES_v2_denovo3umi_round1_round2_w_bulk_joined_BC_pairs_20240106.txt",
                               header=TRUE)

# df_BC_clonotypes %>% select(consensus_clone_id,modality_generation) %>% unique %>% pull(modality_generation) %>% table

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

# # Cre cells
# df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/valid_T7BC_Cre_v3_CS_20231220.txt.gz",
                            # header=TRUE)

# # parental cells: 
df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/clonotype_mapping_v2/valid_T7BC_Parental_v3_CS_20231220.txt.gz",
                            header=TRUE)


# considering BC separately:
BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split("_")
df_T7BC_valid$BC1 <- BC_all %>% lapply("[[",1) %>% unlist()
df_T7BC_valid$BC2 <- BC_all %>% lapply("[[",2) %>% unlist()
df_T7BC_valid$capture_direction <- BC_all %>% lapply("[[",3) %>% unlist()
df_T7BC_valid <- df_T7BC_valid %>% select(-BC)

df_T7BC_valid2_1umi <- df_T7BC_valid %>% group_by(cBC_pdT,BC1,BC2) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC))
df_T7BC_valid2 <- df_T7BC_valid2_1umi %>% filter(UMIs_BC>=2) %>% data.frame()


# # considering BC in pairs
# BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split("_CS")
# df_T7BC_valid$BC_pair <- BC_all %>% lapply("[[",1) %>% unlist()
# df_T7BC_valid$capture_direction <- paste0("CS",BC_all %>% lapply("[[",2) %>% unlist())
# df_T7BC_valid <- df_T7BC_valid %>% select(-BC)
# df_T7BC_valid2_1umi <- df_T7BC_valid %>% group_by(cBC_pdT,BC_pair) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC))
# df_T7BC_valid2 <- df_T7BC_valid2_1umi %>% filter(UMIs_BC>=2)



# length(df_T7BC_valid2_1umi %>% pull(cBC_pdT) %>% unique() )
# length(df_T7BC_valid2 %>% pull(cBC_pdT) %>% unique() )


# df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/valid_T7BC_Cre_v2_pdT_20231210.txt"
#                             header=TRUE)




cellBC <- df_T7BC_valid2 %>% pull(cBC_pdT) %>% unique()


list_df_mapping_to_clonotype <- list()

counter <- 0
tic()
for (cellBC_oi in cellBC){
  
  if ((counter%%100)==0){
    toc()
    print(counter)
    tic()
    
    # if (counter==1000){
    #   break
    # }
  }
  
  # keeping BC separate
  
  # tic()
  
  T7BC1_oi <- df_T7BC_valid2 %>% filter(cBC_pdT %in% cellBC_oi) %>% pull(BC1)
  T7BC2_oi <- df_T7BC_valid2 %>% filter(cBC_pdT %in% cellBC_oi) %>% pull(BC2)
  
  df_frac_from_clones_in_10xBC1 <- df_BC_clonotypes %>% group_by(consensus_clone_id) %>%
    summarize(frac_clone_in_10x1=mean(real_barcode1 %in% T7BC1_oi)) %>% data.frame()
  df_frac_from_clones_in_10xBC2 <- df_BC_clonotypes %>% group_by(consensus_clone_id) %>%
    summarize(frac_clone_in_10x2=mean(real_barcode2 %in% T7BC2_oi)) %>% data.frame()
  
  # toc()
  # 
  # tic()
  df_frac_from_10xBC1_in_clones <- data.frame()
  df_frac_from_10xBC2_in_clones <- data.frame()

  for (clone_oi in list_clones){
    clone_BC1s <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(real_barcode1)
    clone_BC2s <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(real_barcode2)

    df_frac_from_10xBC1_in_clones <- rbind(df_frac_from_10xBC1_in_clones,
                                          data.frame(consensus_clone_id=clone_oi,
                                                     frac_10xBC1_in_clone=mean(T7BC1_oi %in% clone_BC1s)))
    df_frac_from_10xBC2_in_clones <- rbind(df_frac_from_10xBC2_in_clones,
                                           data.frame(consensus_clone_id=clone_oi,
                                                      frac_10xBC2_in_clone=mean(T7BC2_oi %in% clone_BC2s)))
  }
  # toc()

  # combine all information at this stage:

  df_mapping_to_clonotype_oi <- data.frame(cellBC=cellBC_oi,
                                        df_frac_from_clones_in_10xBC1) %>%
    left_join(df_frac_from_clones_in_10xBC2,by="consensus_clone_id") %>%
    left_join(df_frac_from_10xBC1_in_clones,by="consensus_clone_id") %>%
    left_join(df_frac_from_10xBC2_in_clones,by="consensus_clone_id")

  
  
  
  # BC in pairs
  # T7BC_oi <- df_T7BC_valid2 %>% filter(cBC_pdT %in% cellBC_oi) %>% pull(BC_pair)
  # 
  # df_frac_from_clones_in_10x <- df_BC_clonotypes %>% group_by(consensus_clone_id) %>%
  #   summarize(frac_clone_in_10x=mean(BC_pair %in% T7BC_oi)) %>% data.frame()
  # 
  # 
  # df_frac_from_10x_in_clones <- data.frame()
  # 
  # for (clone_oi in list_clones){
  #   clone_BCs <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(BC_pair)
  # 
  #   df_frac_from_10x_in_clones <- rbind(df_frac_from_10x_in_clones,
  #                                         data.frame(consensus_clone_id=clone_oi,
  #                                                    frac_10x_in_clone=mean(T7BC_oi %in% clone_BCs)))
  # }

  # combine all information at this stage:
  # df_mapping_to_clonotype_oi <- data.frame(cellBC=cellBC_oi,
  #                                          df_frac_from_clones_in_10x) %>%
  #   left_join(df_frac_from_10x_in_clones,by="consensus_clone_id")
  
  
  list_df_mapping_to_clonotype[[cellBC_oi]] <- df_mapping_to_clonotype_oi

  
  counter <- counter+1
}




df_mapping_to_clonotype <- rbindlist(list_df_mapping_to_clonotype)


# write.table(df_mapping_to_clonotype,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/mapping_Cre_T7BC_CS_denovo3umi_round2_w_bulk_v2_separate_BC_20240112.txt",
#             quote=FALSE,row.names=FALSE,sep="\t")

write.table(df_mapping_to_clonotype,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/mapping_Parental_T7BC_CS_denovo3umi_round2_w_bulk_v2_separate_BC_20240112.txt",
            quote=FALSE,row.names=FALSE,sep="\t")


# df_mapping_to_clonotype <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/mapping_Cre_T7BC_CS_denovo3umi_round2_w_bulk_v2_20240106.txt",
#            header=TRUE)
# write.table(df_mapping_to_clonotype,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/mapping_Parental_T7BC_CS_denovo3umi_round2_w_bulk_v2_20240106.txt",
#             quote=FALSE,row.names=FALSE,sep="\t")

# df_mapping_to_clonotype <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/mapping_Parental_T7BC_CS_denovo3umi_round2_w_bulk_v2_20240106.txt",
#                                       header=TRUE)


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


# # best recall (capture)
# df_mapping_to_clonotype_best <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_clone_in_10x,n=1,with_ties=FALSE)
# df_mapping_to_clonotype_best2 <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_clone_in_10x,n=2,with_ties=FALSE)

# # best precision (purity)
df_mapping_to_clonotype_best <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_10xBC1_in_clone,n=1,with_ties=FALSE)
df_mapping_to_clonotype_best2 <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_10xBC1_in_clone,n=2,with_ties=FALSE)

# df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x1 != frac_clone_in_10x2)


df_mapping_to_clonotype_best2_only <- anti_join(df_mapping_to_clonotype_best2,df_mapping_to_clonotype_best)
colnames(df_mapping_to_clonotype_best2_only)[c(2,3,4,5,6)] <- c("consensus_clone_id_2ndBest","frac_clone_in_10x1_2ndBest","frac_clone_in_10x2_2ndBest","frac_10xBC1_in_clone_2ndBest","frac_10xBC2_in_clone_2ndBest")
df_mapping_to_clonotype_best_w_2nd_best <- df_mapping_to_clonotype_best %>% left_join(df_mapping_to_clonotype_best2_only)
colnames(df_mapping_to_clonotype_best_w_2nd_best)[c(2,3,4,5,6)] <- c("clone_id","recall_BC1","recall_BC2","precision_BC1","precision_BC2")
colnames(df_mapping_to_clonotype_best_w_2nd_best)[c(7,8,9,10,11)] <- c("clone_id_2nd","recall_BC1_2nd","recall_BC2_2nd","precision_BC1_2nd","precision_BC2_2nd")


# df_mapping_to_clonotype_best <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_clone_in_10x,n=1,with_ties=FALSE)
# df_mapping_to_clonotype_best_precision <- df_mapping_to_clonotype %>% group_by(cellBC) %>% slice_max(frac_10x_in_clone,n=1,with_ties=FALSE)

# df_mapping_to_clonotype_best %>% filter(!(frac_clone_in_10x>0 | frac_10x_in_clone>0)) %>% dim()
# df_mapping_to_clonotype_best %>% filter(!(frac_clone_in_10x>0 | frac_10x_in_clone>0)) %>% dim()


thresh_capt_low <- 0.1
thresh_capt_hi <- 0.3
thresh_precision <- 0.75

max_axis <- 1.05
shading_alpha <- 0.15

dev.set(4)
plt_assignment <- ggplot(df_mapping_to_clonotype_best %>% filter(frac_clone_in_10x1>0 | frac_10xBC1_in_clone>0)) +
  geom_bin_2d(aes(x=frac_clone_in_10x1,y=frac_10xBC1_in_clone))+xlim(c(NA,1))+coord_fixed()+
  annotate(geom="rect",xmin=0,xmax=thresh_capt_low,ymin=0,ymax=max_axis,fill="red",alpha=shading_alpha)+
  annotate(geom="rect",xmin=thresh_capt_low,xmax=max_axis,ymin=0,ymax=thresh_precision,fill="orange",alpha=shading_alpha)+
  annotate(geom="rect",xmin=thresh_capt_low,xmax=max_axis,ymin=thresh_precision,ymax=max_axis,fill="green",alpha=shading_alpha)+
  # geom_vline(xintercept=thresh_capt_low,color="orange",linetype="dashed")+
  # geom_vline(xintercept=thresh_capt_hi,color="red",linetype="dashed")+
  # geom_hline(yintercept=thresh_precision,color="green",linetype="dashed")+
  scale_x_continuous(expand=expansion(0,0),limits=c(-0.01,max_axis))+
  scale_y_continuous(expand=expansion(0,0),limits=c(-0.01,max_axis))+
  labs(x="Recall (fraction barcodes from best\n clone in 10x)",
       y="Precision (fraction barcodes in\n 10x from best clone)")+
  theme_cowplot()+
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=9),
    panel.grid.major=element_line(size=0.1,color="grey"),
    axis.ticks = element_line(size=0.2),
    axis.text=element_text(size=9),
    axis.title=element_text(size=10),
    legend.text=element_text(size=7),
    legend.title=element_text(size=9))
plt_assignment

pdf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/assignment_precision_recall_plane_final_Cre_20240107.pdf",width=3.5,height=3)
print(plt_assignment)
dev.off()


# categories
# df_mapping_to_clonotype_best %>% filter((frac_clone_in_10x>0 | frac_10x_in_clone>0) & frac_clone_in_10x<thresh_capt_low) %>% dim
# df_mapping_to_clonotype_best %>% filter((frac_clone_in_10x>0 | frac_10x_in_clone>0) & frac_clone_in_10x<thresh_capt_low & frac_10x_in_clone>=thresh_precision) %>% dim
# df_mapping_to_clonotype_best %>% filter((frac_clone_in_10x>0 | frac_10x_in_clone>0) & frac_clone_in_10x<thresh_capt_low & frac_10x_in_clone<thresh_precision) %>% dim

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
table(df_mapping_to_clonotype_w_category$assignment_category)



df_mapping_clonotype_final_round1 <- df_mapping_to_clonotype_best_w_2nd_best %>% left_join(df_mapping_to_clonotype_w_category %>% select(cellBC,assignment_category))

# df_mapping_clonotype_final_round1_recall <- df_mapping_clonotype_final_round1
# df_mapping_clonotype_final_round1_precision <- df_mapping_clonotype_final_round1

# df_comp_selection <- df_mapping_clonotype_final_round1_precision %>% select(cellBC,clone_id_PRE=clone_id,assignment_category_PRE=assignment_category) %>% left_join(
#   df_mapping_clonotype_final_round1_recall %>% select(cellBC,clone_id_REC=clone_id,assignment_category_REC=assignment_category)
# )
# df_comp_selection %>% group_by(assignment_category_PRE) %>% summarize(frac_congruent=mean(clone_id_PRE==clone_id_REC))
# df_comp_selection %>% group_by(assignment_category_REC) %>% summarize(frac_congruent=mean(clone_id_PRE==clone_id_REC))
# cells_conflict <- df_comp_selection %>% filter(assignment_category_PRE=="high_capture_high_purity") %>% filter(clone_id_PRE!=clone_id_REC) %>% pull(cellBC)
# cells_conflict <- df_comp_selection %>% filter(assignment_category_PRE=="mid_capture_high_purity") %>% filter(clone_id_PRE!=clone_id_REC) %>% pull(cellBC)


# df_mapping_clonotype_final_round1_precision %>% filter(cellBC %in% cells_conflict) %>% data.frame()
# df_mapping_clonotype_final_round1 %>% filter(recall_2>=thresh_capt_hi & precision_2>=thresh_precision)%>% data.frame()
# df_mapping_clonotype_final_round1 %>% filter(recall_2>=thresh_capt_hi & precision_2>=thresh_precision)%>% data.frame()


df_mapping_clonotype_final_round1_2 <- df_mapping_clonotype_final_round1 #%>% transform(ifelse(assignment_category=="high_capture_high_purity"))
df_mapping_clonotype_final_round1_2$assignment_category[df_mapping_clonotype_final_round1_2$assignment_category=="high_capture_high_purity" &
                                                          df_mapping_clonotype_final_round1_2$recall_BC1_2nd>=0.1] <- "high_capture_high_purity_too_high_recall_2nd"
df_mapping_clonotype_final_round1_2$assignment_category[df_mapping_clonotype_final_round1_2$assignment_category=="mid_capture_high_purity" &
                                                          df_mapping_clonotype_final_round1_2$recall_BC1_2nd>=0.1] <- "mid_capture_high_purity_too_high_recall_2nd"

table(df_mapping_clonotype_final_round1_2$assignment_category)

df_mapping_clonotype_final_round1_3 <- df_mapping_clonotype_final_round1_2 %>% 
  transform(valid_assignment=ifelse(assignment_category %in% c("high_capture_high_purity","mid_capture_high_purity"),TRUE,FALSE))

table(df_mapping_clonotype_final_round1_3$valid_assignment)


# write.table(df_mapping_clonotype_final_round1_3,
#             "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/FINAL_ASSIGNMENTS_round2_Parental_cells_to_clonotypes_denovo3umi_round1_round2_bulk_v2_20240107.txt",
#             row.names=FALSE,sep="\t",quote=FALSE)

# write.table(df_mapping_clonotype_final_round1_3,
#             "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/FINAL_ASSIGNMENTS_round2_Cre_cells_to_clonotypes_denovo3umi_round1_round2_bulk_v2_separate_BC_20240112.txt",
#             row.names=FALSE,sep="\t",quote=FALSE)

write.table(df_mapping_clonotype_final_round1_3,
            "/Users/jbl/Documents/UW/Data/sequencing_runs/seq057SP_reprep_T7_BCs_20231206/assign_cell_to_clonotype_round2/FINAL_ASSIGNMENTS_round2_Parental_cells_to_clonotypes_denovo3umi_round1_round2_bulk_v2_separate_BC_20240112.txt",
            row.names=FALSE,sep="\t",quote=FALSE)

dev.set(4)
# ggplot(df_mapping_clonotype_final_round1_3 %>% filter(valid_assignment)) + 
#   geom_bin_2d(aes(x=recall,y=precision))+
#   facet_wrap(~clone_id)

df_summary_per_clone <- df_mapping_clonotype_final_round1_3 %>% filter(valid_assignment) %>% group_by(clone_id) %>% 
  summarize(
    mean_recall=mean(recall_BC1),
    mean_precision=mean(precision_BC1),
    n_cells=length(cellBC))

df_summary_per_clone2 <- df_summary_per_clone %>% transform(clonotype_cat=ifelse(str_detect(clone_id,"b_"),"bulk",ifelse(str_detect(clone_id,"_r1"),"denovo_sc_round1","denovo_sc_round2")))

df_summary_per_clone2 %>% pull(clonotype_cat) %>% table()

plt1 <- ggplot(df_summary_per_clone2)+stat_bin(aes(x=mean_recall,color=clonotype_cat),geom="step",position="identity")+theme(legend.position="none")
plt2 <- ggplot(df_summary_per_clone2)+stat_bin(aes(x=mean_precision,color=clonotype_cat),geom="step",position="identity")+theme(legend.position="none")
plt3 <- ggplot(df_summary_per_clone2)+stat_bin(aes(x=n_cells,color=clonotype_cat),geom="step",position="identity")+theme(legend.position="none")
plt4 <- ggplot(df_summary_per_clone2) + geom_point(aes(y=mean_precision,x=mean_recall,size=n_cells,color=clonotype_cat),shape=1)
plt1+plt2+plt3+plt4
plt1


# 
# dev.set(4)
# plt1 <- ggplot(df_mapping_clonotype_final_round1) + 
#   stat_bin2d(aes(x=precision,y=precision_2))+
#   facet_wrap(~assignment_category)
# plt2 <- ggplot(df_mapping_clonotype_final_round1) + 
#   stat_bin2d(aes(x=recall,y=recall_2))+
#   facet_wrap(~assignment_category)

# ggplot(df_mapping_clonotype_final_round1_2) +
#   stat_bin2d(aes(y=precision_2,x=recall_2))+
#   facet_wrap(~assignment_category)+
#   theme_cowplot()+
#   theme(
#     strip.background=element_blank(),
#     strip.text=element_text(size=8),
#     panel.grid.major=element_line(size=0.1,color="grey"),
#     axis.ticks = element_line(size=0.2),
#     axis.text=element_text(size=9),
#     axis.title=element_text(size=10),
#     legend.text=element_text(size=7),
#     legend.title=element_text(size=9))

# plt1+plt2+plt3+plot_layout(ncol=1)


