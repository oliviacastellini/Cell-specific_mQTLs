library(dplyr)
library(data.table)

dat <- lapply(1:7, \(x) {fread(paste0("/user/home/uy23281/scratch/rmd.results/model4/bcell/res.cis.", x, ".gz")) }) %>% 
  bind_rows() 
dat <- lapply(1:7, \(x) {fread(paste0("/home/uy23281/model3/EPIC_mQTL/results/05/neu/Neu_estINT.res.cis.", x, ".RData")) }) %>% 
  bind_rows() 

#dat <- lapply(1:7, \(x) {fread(paste0("/home/uy23281/model4/neu/EPIC_mQTL/results/05/res.cis.", x, ".gz")) }) %>%  bind_rows()

dat$FDR <- p.adjust(dat$`p-value`, "fdr")

dat_bcell <- dat %>% 
  filter(FDR < 0.05) %>%
  arrange(`p-value`) %>% 
  group_by(gene) %>%
  slice_head(n=1) %>%
  mutate(celltype="bcell")

# These pairs are significant in at least one cell type - 
# we want their effects estimates from all cell types and all models 
snp_cpg_m3 <- unique(c(paste(dat_neut$SNP, dat_neut$gene),
                    paste(dat_tcell$SNP,dat_tcell$gene),
                    paste(dat_bcell$SNP,dat_bcell$gene),
                    paste(dat_mono$SNP,dat_mono$gene)))

celltypes <- c("neu", "mono", "tcell", "bcell")
snp_cpg<-snp_cpg$x

results_list <- lapply(celltypes, function(x) {
  chunk_list <- lapply(1:7, function(y) {
    file_path <- paste0("/user/home/uy23281/scratch/rmd.results/model1/", x, "/res.cis.", y, ".gz")
    #file_path <- paste0("/home/uy23281/model3/", x, "/EPIC_mQTL/results/05/res.cis.", y, ".gz")
    data <- fread(file_path)
    data_filtered <- data %>%
      filter(paste0(SNP,"_",gene) %in% snp_cpg)
    
    data_filtered$celltype <- x
    data_filtered
  })
  bind_rows(chunk_list)
})

m1_sig <- bind_rows(results_list)
  
save(m1_sig, file="/user/home/uy23281/scratch/rmd.results/model1/model1_extracted.RData")
  
#Extract 20,000 null SNPs for mashr 
dat <- lapply(1:7, \(x) {fread(paste0("/user/home/uy23281/scratch/rmd.results/model4/bcell/res.cis.", x, ".gz")) }) %>% 
  bind_rows() 

dat$FDR <- p.adjust(dat$`p-value`, "fdr")

dat_bcell <- dat %>%
  filter(FDR > 0.05) %>%
  arrange(`p-value`) %>%
  group_by(gene) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  top_n(20000, wt = -`p-value`) %>%
  mutate(celltype = "bcell")



              