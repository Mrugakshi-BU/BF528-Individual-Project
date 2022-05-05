library(readxl)

enriched_terms <- function(list, sheetname) {
  
  list <- tolower(list)
  list_ref <- read_excel("/projectnb/bf528/users/dreadlocks/project-2-dreadlocks/outputs/NIHMS647083-supplement-2.xlsx", sheet = sheetname, skip = 1) %>% 
    na.omit(list_ref)
  list_ref <- list_ref[(list_ref$Category=="GOTERM_CC_FAT" |
                        list_ref$Category=="GOTERM_BP_FAT" |
                        list_ref$Category=="GOTERM_MF_FAT" |
                        grepl("Annotation Cluster", list_ref$Category)),]
  list_ref <- list_ref %>% distinct()
  list_ref <- separate(data = list_ref, col = Term, into = c("GO", "Term"), sep = "~")
  return(list_ref[list_ref$Term %in% list, 'Term'])
}

list <- c('Ion binding',
          'Mitochondrion',
          'Organic acid metabolic process',
          'Phosphorous metabolic process',
          'Identical protein binding')

enriched_terms(list, "1D_COMMON_UP")

list <- c('Ion binding',
          'Regulation of macromolecule biosynthetic process',
          'Regulation of cellular component organization',
          'Chromosome',
          'Chromosome organization')

enriched_terms(list, '1E_COMMON_DWN')