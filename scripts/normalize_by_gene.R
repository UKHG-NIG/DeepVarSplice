library(stringr)
library(plyr)
library(biomaRt)
require(dplyr)
library(tidyverse)
library(tibble)
library(ape)

intron_data <- read.table('intron_inf.txt',sep = '\t',fileEncoding='UTF-8', col.names = c('transcript_id',	'Start'	,'End'))
data_mutation_all <- read.table('shared_only_mutations_gene.maf',sep = '\t',fileEncoding='UTF-8', col.names = c('Hugo_Symbol','Chromosome','Start_Position'))

gff_general <- read.gff("Homo_sapiens.GRCh38.105_intron_added.gff3")

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

GeneID2 <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id"),
                 filters = "hgnc_symbol", values = data_mutation_all$Hugo_Symbol,
                 mart = mart)
                 
transcript_ids_intron <- str_split_fixed(gff_general$attributes[gff_general$type=='intron'], "[;:=]", 7)[,3]
transcript_ids_exon <- str_split_fixed(gff_general$attributes[gff_general$type=='exon'], "[;:=]", 6)[,3]
transcript_ids_cds <- str_split_fixed(gff_general$attributes[gff_general$type=='CDS'], "[;:=]", 7)[,6]
transcript_ids_3primeUTR<- str_split_fixed(gff_general$attributes[gff_general$type=='three_prime_UTR'], "[;:=]", 3)[,3]
transcript_ids_5primeUTR <- str_split_fixed(gff_general$attributes[gff_general$type=='five_prime_UTR'], "[;:=]", 3)[,3]

intron_start <- gff_general$start[gff_general$type=='intron']
intron_end <- gff_general$end[gff_general$type=='intron']
intron_lengths <- intron_end-intron_start+1
intron_df <- data.frame(transcript_ids_intron,intron_start,intron_end)
intron_totals <- aggregate(intron_df['intron_lengths'], by=intron_df['transcript_ids_intron'], sum)
number_of_introns <- ddply(intron_data,.(transcript_ids_intron),nrow)
mean_intron <- aggregate(intron_df['intron_lengths'], by=intron_df['transcript_ids_intron'], mean)
intron_df <- merge(intron_totals, mean_intron, by.x="transcript_ids_intron", by.y = "transcript_ids_intron", all=T)
colnames(intron_df)[2] <- 'total_intron_length'
colnames(intron_df)[3] <- 'mean_intron_length'
intron_df2 <- merge(intron_df, number_of_introns, by.x="transcript_ids_intron", by.y = "transcript_ids_intron", all=T)
colnames(intron_df2)[4] <- 'total_intron_number'
merged_data_introns <- merge(intron_df2,GeneID2, by.x="transcript_ids_intron", by.y = "ensembl_transcript_id", all=T)
merged_data_introns<-merged_data_introns[complete.cases(merged_data_introns), ]
merged_data_introns <- merged_data_introns[,2:5]

rec_introns <- aggregate(merged_data_introns['total_intron_number'], by=merged_data_introns['hgnc_symbol'], mean)
rec_introns_mean <- aggregate(merged_data_introns['total_intron_length'], by=merged_data_introns['hgnc_symbol'], mean)
final_data <- merge(rec_introns, intron_df, by.x="hgnc_symbol", by.y = "hgnc_symbol", all=T)

final_data_all <- merge(data_mutation, final_data, by.x="Hugo_Symbol", by.y = "hgnc_symbol", all=T)
final_data_all<-final_data_all[complete.cases(final_data_all), ]
final_data_bar <- merge(GeneID2, intron_df, by.x="ensembl_transcript_id", by.y = "transcript_ids_intron", all=F)


GeneID_for_gene_length <- getBM(attributes = c("hgnc_symbol","start_position","end_position"),
                 filters = "hgnc_symbol", values = final_data_all$Hugo_Symbol,
                 mart = mart) #data$gene was the previous one filter was ensembl_transcript_id

GeneID_for_gene_length$gene_length <- GeneID_for_gene_length$end_position-GeneID_for_gene_length$start_position
final_data_all2 <- merge(GeneID_for_gene_length[,c(1,4)], final_data_all, by.x="hgnc_symbol", by.y = "Hugo_Symbol", all=F)
final_data_all2$normalized_intron_length <- final_data_all2$gene_length/final_data_all2$total_intron_length
final_data_all2$normalized_deep_muts <- final_data_all2$Number_Deep/final_data_all2$normalized_intron_length  
final_data_all2$normalized_ss_muts <- final_data_all2$Number_SS/final_data_all2$total_intron_number  
final_data_all2$normalized_prox_muts <- final_data_all2$Number_Prox/final_data_all2$total_intron_number

res_data <- final_data_all2[!duplicated(final_data_all2$hgnc_symbol),]

mutation_num <- ddply(data_mutation,.(Hugo_Symbol),nrow)

