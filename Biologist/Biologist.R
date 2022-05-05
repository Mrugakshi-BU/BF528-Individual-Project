##Author: Mrugakshi Chidrawar
##BF528 Individual Project - Biologist role

#import libraries
library(pheatmap)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

########### 7.1 ########### 

#reading the tables
P0_1 <- read.table("/projectnb/bf528/students/mruga77/P0_1_cufflinks/genes.fpkm_tracking",header = T)
P0_2 <- read.table("/project/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking",header = T)
AD_1 <- read.table("/project/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking",header = T)
AD_2 <- read.table("/project/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking",header = T)
P4_1 <- read.table("/project/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking",header = T)
P4_2 <- read.table("/project/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking",header = T)
P7_1 <- read.table("/project/bf528/project_2/data/samples/P7_1/genes.fpkm_tracking",header = T)
P7_2 <- read.table("/project/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking",header = T)

# Lists of important genes
sarcomere_genes <- c('Pdlim5', 'Pygm', 'Myoz2', 'Des', 'Csrp3', 'Tcap', 'Cryab')
mitochondria_genes <- c("Mpc1","Prdx3","Acat1","Echs1","Slc25a11","Phyh")
cell_cycle_genes <- c("Cdc7","E2f8","Cdk7","Cdc26","Cdc6","Cdc27","E2f1","Cdc45","Rad51","Aurkb","Cdc23")

#Plot FPKM values for Sacromere genes
#subset sarcomere genes between all the samples
sarcomere_genes_fpkm <- c(P0_1[P0_1$gene_short_name %in% sarcomere_genes,"FPKM"],
                P0_2[P0_2$gene_short_name %in% sarcomere_genes,"FPKM"],
                P4_1[P4_1$gene_short_name %in% sarcomere_genes,"FPKM"],
                P4_2[P4_2$gene_short_name %in% sarcomere_genes,"FPKM"],
                P7_1[P7_1$gene_short_name %in% sarcomere_genes,"FPKM"],
                P7_2[P7_2$gene_short_name %in% sarcomere_genes,"FPKM"],
                AD_1[AD_1$gene_short_name %in% sarcomere_genes,"FPKM"],
                AD_2[AD_2$gene_short_name %in% sarcomere_genes,"FPKM"])

#creating dataframe for the fpkm values of sarcomere genes
sarcomere_df <- data.frame("FPKM" = sarcomere_genes_fpkm, "Sample" = rep(c("P0","P4","P7","Ad"), each = 14),
                           "Gene" = rep(P0_1[P0_1$gene_short_name %in% sarcomere_genes, "gene_short_name"], 8))
sarcomere_df$Sample <- factor(sarcomere_df$Sample, levels = c("P0","P4","P7","Ad"))

#creating plot data for fpkm values of sacromere genes
sarcomere_genes_plot_data <- sarcomere_df %>% 
  group_by(Sample, Gene) %>% 
  summarise("Mean" = mean(FPKM), "SD" = sd(FPKM))
  
#create plot
ggplot(sarcomere_genes_plot_data, aes(x = Sample, y = Mean, group = Gene)) + 
  geom_line(aes(color = Gene)) + 
  geom_point(aes(color = Gene)) + 
  labs(title="Sarcomere", x ="Sample", y = "FPKM") + 
  theme_bw()

#Plot FPKM values for mitochodria genes
##subset mitochondria genes between all the samples
mitochondria_genes_fpkm <- c(P0_1[P0_1$gene_short_name %in% mitochondria_genes,"FPKM"],
                 P0_2[P0_2$gene_short_name %in% mitochondria_genes,"FPKM"],
                 P4_1[P4_1$gene_short_name %in% mitochondria_genes,"FPKM"],
                 P4_2[P4_2$gene_short_name %in% mitochondria_genes,"FPKM"],
                 P7_1[P7_1$gene_short_name %in% mitochondria_genes,"FPKM"],
                 P7_2[P7_2$gene_short_name %in% mitochondria_genes,"FPKM"],
                 AD_1[AD_1$gene_short_name %in% mitochondria_genes,"FPKM"],
                 AD_2[AD_2$gene_short_name %in% mitochondria_genes,"FPKM"])

##converting the subset to table
mitochondria_df <- data.frame("FPKM" = mitochondria_genes_fpkm, "Sample" = rep(c("P0","P4","P7","Ad"), each = 10),
                         "Gene" = rep(P0_1[P0_1$gene_short_name %in% mitochondria_genes,"gene_short_name"], 8))
mitochondria_df$Sample <- factor(mitochondria_df$Sample, levels = c("P0","P4","P7","Ad"))

#creating plot data for fpkm values of mitochodria genes
mitochondria_genes_plot_data <- mitochondria_df %>% 
  group_by(Sample, Gene) %>% 
  summarise("Mean" = mean(FPKM), "SD" = sd(FPKM))

#create plot
ggplot(mitochondria_genes_plot_data, aes(x = Sample, y = Mean, group = Gene)) + 
  geom_line(aes(color = Gene)) + 
  geom_point(aes(color = Gene)) + 
  labs(title="Mitochondria", x ="Sample", y = "FPKM") + 
  theme_bw()

#Plot FPKM values for cell cycle genes
#subset cell cycle genes between all the samples
cell_cycle_genes_fpkm <- c(P0_1[P0_1$gene_short_name %in% cell_cycle_genes,"FPKM"],
                 P0_2[P0_2$gene_short_name %in% cell_cycle_genes,"FPKM"],
                 P4_1[P4_1$gene_short_name %in% cell_cycle_genes,"FPKM"],
                 P4_2[P4_2$gene_short_name %in% cell_cycle_genes,"FPKM"],
                 P7_1[P7_1$gene_short_name %in% cell_cycle_genes,"FPKM"],
                 P7_2[P7_2$gene_short_name %in% cell_cycle_genes,"FPKM"],
                 AD_1[AD_1$gene_short_name %in% cell_cycle_genes,"FPKM"],
                 AD_2[AD_2$gene_short_name %in% cell_cycle_genes,"FPKM"])

#converting the subset to table
cell_cycle_df <- data.frame("FPKM" = cell_cycle_genes_fpkm, "Sample" = rep(c("P0","P4","P7","Ad"), each = 22),
                        "Gene" = rep(P0_1[P0_1$gene_short_name %in% cell_cycle_genes,"gene_short_name"], 8))
cell_cycle_df$Sample <- factor(cell_cycle_df$Sample, levels = c("P0","P4","P7","Ad"))

#creating plot data for fpkm values of cell cycle genes
cell_cycle_genes_plot_data <- cell_cycle_df %>% 
  group_by(Sample, Gene) %>% 
  summarise("Mean" = mean(FPKM), "SD" = sd(FPKM))

#create plot
ggplot(cell_cycle_genes_plot_data, aes(x = Sample, y = Mean, group = Gene)) + 
  geom_line(aes(color = Gene)) + 
  geom_point(aes(color = Gene)) + 
  labs(title="Cell Cycle", x ="Sample", y = "FPKM") + 
  theme_bw()


########### 7.3 ########### 

fpkm_matrix <- read.table("/project/bf528/project_2/data/fpkm_matrix.csv", header = T)
cuffdiff_gene_exp <- read.table("/projectnb/bf528/students/mruga77/cuffdiff_out/gene_exp.diff", header = T)

#Filter at most top 1k genes found to be differentially expressed between P0 and Ad 
cuffdiff_gene_exp <- arrange(cuffdiff_gene_exp, q_value)
cuffdiff_gene_exp_1000 <- head(cuffdiff_gene_exp, 1000)

top_id <- P0_1[P0_1$gene_short_name %in% cuffdiff_gene_exp_1000$gene,"tracking_id"] %>% unique()

df_join <- P0_1[, c("tracking_id", "FPKM")] %>%
  inner_join(fpkm_matrix, by = "tracking_id")

colnames(df_join)[2] <- "P0_1_FPKM"
top_1000 <- df_join[df_join$tracking_id %in% top_id,c(1:5)]
top_1000 <- top_1000[rowSums(top_1000[,-1]) != 0, ]

#creating heatmap
pheatmap(as.matrix(top_1000[,-1]), scale = "row",show_rownames = F)


