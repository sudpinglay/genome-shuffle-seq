

# getting the starting data 
df_BC_clonotypes <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/clonotype_calling_Lane1_round2/FINAL_denovo_clonotypes_round2_Lane1_3umi_table_CS1_BC_RERUN_20241007.txt",
                               header=TRUE) %>% transform(CS1_BC=str_replace(CS1_BC,"-","_"))

df_mapping_clonotype_final_round1_3 <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/assignment_cell_to_clonotype_Lane1_round2/FINAL_round2_assignment_cells_to_denovo3umi_Lane1_RERUN_20241007.txt",
                                                  header=TRUE)

# only retain 'major' clonotypes with >10 cells assigned
df_summary_assignment <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/assignment_cell_to_clonotype_Lane1_round2/summary_stats_round2_assignment_Lane1_RERUN_20241007.txt",
                                    header=TRUE)

major_clonotypes <- df_summary_assignment %>% filter(n_cells>10) %>% pull(clone_id)




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
df_T7BC_valid3_1umi <- df_T7BC_valid3 %>% group_by(cBC_pdT,CS1_BC) %>% summarize(reads_BC=sum(reads_BC),UMIs_BC=sum(UMIs_BC))


# precision recall analysis on final round1 set. 
umi_ct_thresh <- seq(1,10)
df_precision_recall_high_conf_clones <- data.frame()

df_valid_assignment <- df_mapping_clonotype_final_round1_3 %>% filter(valid_assignment & clone_id %in% major_clonotypes)
for (id in seq(dim(df_valid_assignment)[1])){
  
  if ((id %% 100)==0){
    print(id)
  }
  
  cellBC_oi <- df_valid_assignment$cellBC[id]
  clone_oi <- df_valid_assignment$clone_id[id]
  clone_BCs_oi <- df_BC_clonotypes %>% filter(consensus_clone_id==clone_oi) %>% pull(CS1_BC)
  
  df_clone_BCs_oi <- df_BC_clonotypes %>% transform(CS1_BC=str_replace(CS1_BC,"-","_")) %>% filter(consensus_clone_id==clone_oi)
  
  df_all_capture_confident_clones <- df_T7BC_valid3_1umi %>% filter(cBC_pdT %in% cellBC_oi) %>% 
    transform(in_best_clone= (CS1_BC %in% clone_BCs_oi))

  precision_oi <- c()
  recall_oi <- c()
  for (thres in umi_ct_thresh){
    precision_oi[thres] <- mean(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres) %>% pull() %>% mean())
    recall_oi[thres] <- mean(!is.na(df_clone_BCs_oi %>% 
                                      left_join(df_all_capture_confident_clones %>% filter(UMIs_BC>=thres),by = join_by(CS1_BC)) %>% pull(UMIs_BC)))
  }
  
  df_PR_oi <- data.frame(cellBC=cellBC_oi,
                         high_conf_clone=clone_oi,
                         umi_thresh=umi_ct_thresh,
                         precision=precision_oi,
                         recall=recall_oi)
  
  
  df_precision_recall_high_conf_clones <- rbind(df_precision_recall_high_conf_clones,df_PR_oi)
  
}




write.table(df_precision_recall_high_conf_clones,
            "/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/assignment_cell_to_clonotype_Lane1_round2/table_PR_high_confidence_round2_Lane1_RERUN_assignments_20241008.txt",
            row.names=FALSE, quote = FALSE)


df_av_PR_all <- df_precision_recall_high_conf_clones %>% group_by(umi_thresh) %>% 
  summarize(median_prec=median(precision,na.rm=TRUE),
            mean_prec=mean(precision,na.rm=TRUE),
            median_recall=median(recall,na.rm=TRUE),
            mean_recall=mean(recall,na.rm=TRUE))



dev.set(4)
plt_PR <- ggplot(df_av_PR_all) + 
  geom_line(aes(x=umi_thresh,y=mean_prec),color="firebrick",linewidth=0.5)+
  geom_line(aes(x=umi_thresh,y=mean_recall),color="dodgerblue",linewidth=0.5)+
  geom_point(aes(x=umi_thresh,y=mean_prec),color="firebrick",size=1)+
  geom_point(aes(x=umi_thresh,y=mean_recall),color="dodgerblue",size=1)+
  # coord_fixed()+
  theme_cowplot()+
  # scale_x_continuous(limits=c(0,1.01),expand=expansion(0,0.05))+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
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

pdf("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/Precision_Recall_Lane1_clonotypes_assignment_20241008.pdf",width=3,height=2.5)
print(plt_PR)
dev.off()
