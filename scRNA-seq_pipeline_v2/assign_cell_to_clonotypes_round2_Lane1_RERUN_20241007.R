

library(data.table)


# load consensus clonotype BC table:


df_BC_clonotypes <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling_Lane1_round2/FINAL_denovo_clonotypes_round2_Lane1_3umi_table_CS1_BC_RERUN_20241007.txt",
                               header=TRUE) %>% transform(CS1_BC=str_replace(CS1_BC,"-","_"))


df_BC_clonotypes %>% group_by(consensus_clone_id) %>% summarize(MOI=length(CS1_BC)) %>% arrange(MOI)


list_clones <- df_BC_clonotypes %>% pull(consensus_clone_id) %>% unique()

df_T7BC_valid <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling_Lane1_round1/valid_T7BC_Lane1_RERUN_updated_cells_valid_species_singleton_20241007.txt.gz",
                            header=TRUE)

BC_all <- df_T7BC_valid %>% pull(BC) %>% str_split(regex("_att.__"))

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
CS1_BC <- df_T7BC_valid2 %>% pull(BC_pair) %>% str_split("__") %>% lapply("[[",1) %>% unlist() 
df_T7BC_valid3 <- df_T7BC_valid2
df_T7BC_valid3$CS1_BC <- CS1_BC
df_T7BC_valid3 <- df_T7BC_valid3 %>% group_by(cBC_pdT,CS1_BC) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC))

# restrict to 2 umi or more to decrease cross-detection: 
df_T7BC_valid_final <- df_T7BC_valid3 %>% filter(UMIs_BC>1)



# 
# # exclude BC from large clone in a second iteration. 
# BC_large_clones <- df_BC_clonotypes %>% filter(str_detect(consensus_clone_id,"large")) %>% pull(CS1_BC)
# df_T7BC_valid_final <- df_T7BC_valid_final %>% filter(!(CS1_BC %in% BC_large_clones))

cellBC <- df_T7BC_valid_final %>% pull(cBC_pdT) %>% unique()


list_df_mapping_to_clonotype <- list()

counter <- 0
# tic()
for (cellBC_oi in cellBC){
  
  if ((counter%%100)==99){
    print(counter)
   
  }
  
  
  
  T7BC_oi <- df_T7BC_valid_final %>% filter(cBC_pdT %in% cellBC_oi) %>% pull(CS1_BC)

  df_frac_from_clones_in_10x <- df_BC_clonotypes %>% group_by(consensus_clone_id) %>%
    summarize(frac_clone_in_10x=mean(CS1_BC %in% T7BC_oi)) %>% data.frame()
  

  df_frac_from_10x_in_clones <- data.frame()

  for (clone_oi in list_clones){
    clone_BCs <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(CS1_BC)

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

write.table(df_mapping_to_clonotype,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/assignment_cell_to_clonotype_Lane1_round2/mapping_Lane1_round2_denovo3umi_RERUN_20241007.txt",
            quote=FALSE,row.names=FALSE,sep="\t")



# best precision (purity) --> take that one!
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

# df_mapping_clonotype_final_round1_recall <- df_mapping_clonotype_final_round1
df_mapping_clonotype_final_round1_precision <- df_mapping_clonotype_final_round1



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
            "/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/assignment_cell_to_clonotype_Lane1_round2/FINAL_round2_assignment_cells_to_denovo3umi_Lane1_RERUN_20241007.txt",
            row.names=FALSE,sep="\t",quote=FALSE)


df_mapping_clonotype_final_round1_3 <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/assignment_cell_to_clonotype_Lane1_round2/FINAL_round2_assignment_cells_to_denovo3umi_Lane1_RERUN_20241007.txt",
                                                  header=TRUE)
dev.set(4)
ggplot(df_mapping_clonotype_final_round1_3 %>% filter(valid_assignment)) + 
  geom_bin_2d(aes(x=recall,y=precision))+
  facet_wrap(~clone_id)

df_summary_per_clone <- df_mapping_clonotype_final_round1_3 %>% filter(valid_assignment) %>% group_by(clone_id) %>% 
  summarize(
    mean_recall=mean(recall),
    mean_precision=mean(precision),
    n_cells=length(cellBC))

write.table(df_summary_per_clone,"/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/assignment_cell_to_clonotype_Lane1_round2/summary_stats_round2_assignment_Lane1_RERUN_20241007.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

 plt1 <- ggplot(df_summary_per_clone)+stat_bin(aes(x=mean_recall),geom="step")
plt2 <- ggplot(df_summary_per_clone)+stat_bin(aes(x=mean_precision),geom="step")
plt3 <- ggplot(df_summary_per_clone)+stat_bin(aes(x=n_cells),geom="step")+scale_x_log10()
plt4 <- ggplot(df_summary_per_clone) + geom_point(aes(y=mean_precision,x=mean_recall,size=n_cells,color=n_cells<10),shape=1)
dev.new()
dev.set(4)
plt1+plt2+plt3+plt4

df_summary_per_clone %>% arrange(mean_precision)


dev.set(4)
plt1 <- ggplot(df_mapping_clonotype_final_round1) +
  stat_bin2d(aes(x=precision,y=precision_2))+
  facet_wrap(~assignment_category)
plt2 <- ggplot(df_mapping_clonotype_final_round1) +
  stat_bin2d(aes(x=recall,y=recall_2))+
  facet_wrap(~assignment_category)
plt3 <- ggplot(df_mapping_clonotype_final_round1) +
  stat_bin2d(aes(y=precision_2,x=recall_2))+
  facet_wrap(~assignment_category)
plt1+plt2+plt3+plot_layout(ncol=1)
