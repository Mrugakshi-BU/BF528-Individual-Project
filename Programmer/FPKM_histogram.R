#import libraries
library('tidyverse')
library(ggplot2)

#data processing
genes_fpkm <- read.table('/projectnb/bf528/students/mruga77/P0_1_cufflinks/genes.fpkm_tracking')
colnames(genes_fpkm) <- genes_fpkm[1,] #column names are in the first row
genes_fpkm <- genes_fpkm[-1,] #get rid of that first row
genes_fpkm$log_FPKM <- log10(as.numeric(genes_fpkm$FPKM)) #convert to numeric

#creating histogram
filtered_genes_fpkm <- genes_fpkm[which(genes_fpkm$log_FPKM > 0), ] #removing the genes that have an FPKM value of zero
hist(filtered_genes_fpkm$log_FPKM, xlab = 'log10(FPKM)', main = 'Histogram of log adjusted FPKM values', breaks = 22)

dim(genes_fpkm)
dim(filtered_genes_fpkm)





