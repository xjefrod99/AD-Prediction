install.packages('tidyverse')

install.packages("stringr")
library("stringr")
library(tidyverse)

setwd("~/Desktop/AI_project")

donor_info <- read.csv('Donor_info.csv')
donor_info <- subset(donor_info, select = - c(2,3,4,6,8,11,15:17) ) 
gene_rows <- read.csv('rows-genes.csv')

#only grabbing genes having to do mitochondria & hippocampus
#genes_of_interes <- gene_rows %>%
#                  filter(str_detect(gene_name, 'mitochondria|hippocampus|mitochondria|factor|nucleus|protein|micro'))

#genes_of_interes <- gene_rows %>%
#                  filter(str_detect(gene_name, 'mitochondria|hippocampus|mitochondria|factor|actin'))



genes_of_interes <- gene_rows %>%
                  filter(str_detect(gene_name, 
                  'mitochondria|hippocampus|tay|factor|cis|MAP kinase|
                  CLU|VASP|phosphoprotein|tyrosine phosphatase|PTPRH'))



list_genes <- genes_of_interes$gene_id

#in fpkm, rows are genes IDs, cols 

#contains info regarding where the samples came from, like donor name and brain area
columns <- read.csv('columns-samples.csv')

#cols_of_interes <- columns %>%
  #filter(str_detect(structure_name, 'hippocampus|hippocampal|temporal'))


cols_of_interes <- columns 
data <- read.csv('fpkm_table_normalized.csv')
#rows are genes
#cols are samples
data_cols <- ( colnames(data))

data_cols <- str_replace_all(pattern="X", replacement="", data_cols)
data<- set_names(data, data_cols)


#data_of_interest <- inner_join(data,genes_of_interes, by = c('rnaseq_profile_id'),copy = FALSE , suffix = c('.x', '.y'))

data_of_interest<- subset(data, rnaseq_profile_id %in% list_genes)
col.names <- colnames(data_of_interest) #sample names 
row.names <-data_of_interest$rnaseq_profile_id #genes 
row.names = append(row.names, "samples", 0)
#rows are genes
#cols are samples

#we want genes as columns and cols are samples bc we're trying to predict sample label given all the genes
final_df <- as.data.frame(t(data_of_interest))

final_df <- cbind(newColName = rownames(final_df), final_df)

colnames(final_df)<- row.names


final_df <- final_df[-1,] 
#we have 382 obsevations, 610 RNA  features , now we append label data


final_df$samples <- as.integer(final_df$samples)

cols_of_interes$rnaseq_profile_id <-as.integer(cols_of_interes$rnaseq_profile_id)

final_df <- left_join(final_df, cols_of_interes, by = c('samples' = 'rnaseq_profile_id'), copy = FALSE , suffix = c('.x', '.y'))

final_df <- left_join(final_df, donor_info, by = c('donor_id') )


a <- final_df %>%
  mutate(Dementia = case_when( ('Dementia'%in% nincds_arda_diagnosis & !('No Dementia' %in% nincds_arda_diagnosis) )  ~ 1,
                              TRUE ~ 0) ) %>%
  mutate(Alzheimer = case_when( (nincds_arda_diagnosis == "Probable Alzheimer'S Disease") | (nincds_arda_diagnosis == "Possible Alzheimer'S Disease") ~ 1,
                               TRUE ~0)) %>%
  mutate(APOE4 = case_when('Y' %in% apo_e4_allele ~ 1,
                               TRUE ~0)) %>%
  mutate(TBI = case_when('Y' %in% ever_tbi_w_loc ~ 1,
                           TRUE ~0)) 
  
a <- subset(a, select = -c(3817:3823))

write.csv(a, 'samples_org_dataset_3kfeatures.csv')


b = read.csv('PCA_important.txt')


b$gene_id <-as.integer(b$gene_id)

pca_genes <-  left_join(b, genes_of_interes, by = c('gene_id') )

pca_gene_names <- pca_genes$gene_name

write.csv(pca_genes, 'pca_genes_in10pcs.csv')

write.csv(pca_gene_names, 'pca_10importantgenes_in10pcs.csv')

