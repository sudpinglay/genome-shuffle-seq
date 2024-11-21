
library(tidyverse)
library(stringr)
library(Matrix)
library(Biostrings)


# read T7BC data

sample_oi <- "Lane2"
df_T7BC_CS <- read.table("/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_shuffleBC/Lane2_FB_CS_S2_w_47bpR1R2transfer_get_R2_2BC_umi_cleaned_20240904.txt.gz",
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

# starting
# CGTAAGTAGGATAACC CS-cellBC
# CGTAAGTTCGATAACC JBs converted
# CGTAAGTCCGATAACC Suds converted



# get set of "bona fide cells" from the GEx criteria 
file_valid_cellBC <- "/Users/jbl/Documents/UW/Data/sequencing_runs/seq075_SP_10x_shuffle_BC_revision_20240831/prepro_GEx/metadata_GEx_cellBC_SP_shuffleBC_revision_20240906.txt"
df_valid_cellBC_all <- read.table(file_valid_cellBC,sep="\t",header=TRUE) #%>% transform(raw_cellBC_pdT=paste0(raw_cellBC_pdT,))

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

ggplot(df_T7BC_CS3) + stat_bin2d(aes(x=reads_BC,y=UMIs_BC))+
  scale_x_log10()+scale_y_log10()+geom_abline(intercept=-log10(8))+coord_fixed()

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

