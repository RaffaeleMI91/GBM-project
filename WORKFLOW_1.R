# PRE-PROCESSED SANGER-BROAD GEX 
library(pbapply)
library(dplyr)
library(tidyr)
library(fgsea)
library(doParallel)
library(pheatmap)
library(ggrepel)
library(purrr)
library(stringr)
library(rstatix)
library(coin)

Sanger_expr <- read.csv("~/Downloads/rnaseq_all_20220624/rnaseq_tpm_20220624.csv", stringsAsFactors = FALSE, header = TRUE)
tpm <- Sanger_expr %>% dplyr::filter(model_id != "SIDG06476" & model_id != "SIDG03156")
genes <- tpm$X[5:length(tpm$X)]
tpm <- tpm[,(tpm[2,] == "Sanger & Broad Cell Lines RNASeq")]
colnames(tpm) <- tpm[1,]
tpm <- tpm[-c(1:4),]
rownames(tpm) <- genes
dim(tpm)
rm(Sanger_expr)

SD_GEX <- apply(tpm, 1, sd)
threshold <- quantile(SD_GEX, 0.35, na.rm=T)
HVG <- names(SD_GEX[SD_GEX >= threshold])
tpm <- tpm[which(rownames(tpm) %in% HVG),]
genes <- rownames(tpm) 
tpm <- apply(tpm,2,as.numeric)
rownames(tpm) <- genes

#tpm <- tpm[rowSums(tpm >= 1) >= (ncol(tpm) * 0.5), ]
#genes <- rownames(tpm)
#tpm <- apply(tpm,2,as.numeric)
#rownames(tpm) <- genes

#load("~/Downloads/GDSC_MATRIX.RData")
load("~/Downloads/GEX_drugRES_correlation_Ranks.RData")

#dim(GDSC_MATRIX)
#dim(tpm)

#cor_matrix <- pbapply(GDSC_MATRIX, 1, function(x) {
#  apply(tpm, 1, function(y) {
#    cor(x, y, use = "complete.obs")
#  })
#})

cor_matrix <- GEX_drugRES_correlation_Ranks
rm(GEX_drugRES_correlation_Ranks)

#save(cor_matrix, file="/Volumes/iorio/Raffaele/Project_ADD3/review/data/IC50_EXPR_CORR_MAT.RData")
#load("/Volumes/iorio/Raffaele/Project_ADD3/review/data/IC50_EXPR_CORR_MAT.RData")
drug_names <- str_replace(colnames(cor_matrix), "\\(gdsc1\\)","")
drug_names <- str_replace(drug_names, "\\(gdsc2\\)","")
drug_names <- str_replace(drug_names, " [^ ]*$","")
colnames(cor_matrix) <- drug_names
# load DEGs

file_path <- "/Volumes/iorio/Raffaele/Project_ADD3/TOT_DEGs_CELL_TYPE_NEFTEL_bis.csv"

deg_df <- read.csv(file=file_path, stringsAsFactors = FALSE, header=TRUE) %>%
  filter(pvals_adj < 0.05 & (logfoldchanges > 5 | logfoldchanges < -5)) %>%
  mutate(
    class = ifelse(logfoldchanges > 5, "Upregulated", "Downregulated")
  ) %>%
  group_by(sample, cell_type, class) %>%
  arrange(sample, cell_type, class, desc(abs(logfoldchanges))) %>% 
  slice_head(n = 300)

LIST_BIS <- lapply(unique(deg_df$sample), function(s){
  deg_by_ct <- deg_df %>% filter(sample==s)
  ct_annotation <- unique(deg_by_ct$cell_type)
  deg_list_by_ct <- lapply(ct_annotation, function(ct){
    degs <- deg_by_ct %>% filter(cell_type==ct) %>% dplyr::select(gene) %>% pull()
    class <- deg_by_ct %>% filter(cell_type==ct) %>% dplyr::select(class) %>% pull()
    my_gmt <- NULL
    my_gmt <- rbind(my_gmt, data.frame(
      CT=c(rep(ct, length(degs))),
      NAME=degs,
      CLASS=class,
      stringsAsFactors = FALSE))
    return(my_gmt)
  })
  return(deg_list_by_ct)
})
  
names(LIST_BIS) <- unique(deg_df$sample)

MY_GMT_LIST <- lapply(LIST_BIS, function(x) {
  my_gmt <- as.data.frame(do.call(rbind, x))
  return(my_gmt)
}
)

lapply(names(MY_GMT_LIST), function(sample){
  gmt_data <- MY_GMT_LIST[[sample]] %>%
    group_by(CT, CLASS) %>%
    summarize(
      genes = paste(NAME, collapse = "\t"),  # Tab-separated genes
      .groups = "drop")
  gmt_lines <- apply(gmt_data, 1, function(row) {
    set_name <- paste(row["CT"], row["CLASS"], sep = "_")  # Combine CT and CLASS as set name
    paste(set_name, "NA", row["genes"], sep = "\t")
  })
  writeLines(gmt_lines, paste0("/Users/raffaele.iannuzzi/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/",sample,"_sclevel_degs.gmt"))
})

file_list <- list.files("/Users/raffaele.iannuzzi/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/")
idx <- grep("_sclevel_degs.gmt",file_list)
gmt_list <- file_list[idx]

myGO_LIST <- lapply(gmt_list, function(x) {
  GO_file <- paste0("/Users/raffaele.iannuzzi/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/",x)
  myGO <- fgsea::gmtPathways(GO_file)
  return(myGO)
})

names(myGO_LIST) <- c("C3N-00662", "C3N-01814", "C3N-02190", "C3N-02783","C3N-03186", "C3N-03188")

sorted_list <- lapply(colnames(cor_matrix),function(i){
  ref <- cor_matrix[,i]
  idx <- order(ref,decreasing = T)
  vec <- ref[idx]
  names(vec) <- rownames(cor_matrix)[idx]
  return(vec)
})

names(sorted_list) <- colnames(cor_matrix)

sorted_list <- lapply(sorted_list, function(x) {x[is.finite(x)]})

num_cores <- detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)

run_fgsea <- function(drug_name, rank, query_sets, minSize=2, maxSize=300, nperm=1000) {
  result <- fgsea(
    pathways=query_sets, 
    stats=rank, 
    minSize=minSize, 
    maxSize=maxSize,
    nperm=nperm)
  result$condition <- drug_name 
  return(result)
}

fgsea_results <- foreach(i=names(sorted_list), .combine = rbind, .packages = "fgsea") %dopar% {
  ranked_list <- sorted_list[[i]]
  res_list <- lapply(
    myGO_LIST, function(x){
      run_fgsea(
        drug_name=i, 
        rank=ranked_list, 
        query_sets=x)})
  return(res_list)
}

stopCluster(cl)

df_fgsea <- do.call(
  rbind.data.frame, 
  lapply(1:ncol(fgsea_results), function(col){
    df <- do.call(rbind, fgsea_results[,col]) %>%
    as.data.frame() %>% 
    #filter(padj < 0.05) %>%
    mutate(NES=ifelse(padj > 0.05, 0, NES)) %>% # same as in DREEP
    mutate(TCGA_sample=colnames(fgsea_results)[col]) %>%
    mutate(id=paste0(pathway,"_",TCGA_sample)) %>%
    separate(pathway, into=c("ctype","class"),sep = "_")
    return(df)
  }
))

NES_matrix <- df_fgsea %>%
  dplyr::select(id,condition,NES) %>% 
  pivot_wider(names_from=id, values_from=NES, values_fn = median)

DRUGS <- NES_matrix$condition
CL_ID <- colnames(NES_matrix)[-1]
NES_matrix$condition <- NULL
NES_matrix <- t(NES_matrix) 
colnames(NES_matrix) <- DRUGS

write.csv(x = NES_matrix, file="/Volumes/iorio/Raffaele/Project_ADD3/review/NES_MAT_DRESP_FROMR.csv")

library(heatmaply)

plot(hclust(dist(NES_matrix)))

mat <- as.matrix(dist(NES_matrix))
scaled_mat <- (mat - min(mat)) / (max(mat) - min(mat))
colnames(mat)

hc <- hclust(dist(NES_matrix), method = "ward.D2")

# Convert to dendrogram and plot with ggplot2
dendro <- ggdendrogram(hc, theme_dendro = TRUE)
dendro
ggsave("~/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/NES_DRESP_DENDRO.pdf",dendro, width = 7, height = 6)

cell_type_ann <- data.frame(idx=rownames(NES_matrix)) %>% separate(idx, into = c("cell type", "class", "sample"), sep = "_")
rownames(cell_type_ann) <- rownames(NES_matrix)
cell_type_ann$sample <- NULL

set.seed(666)

pheatmap(NES_matrix,
         cluster_cols =T,
         cluster_rows = T,
         border_color = "black",
         annotation_row = cell_type_ann,
         show_rownames = F, 
         show_colnames = F)#,
         #filename="~/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/pheatmap_neftel_sclevel_degs_NES.pdf")

##

df_fgsea_up <- df_fgsea %>% filter(class=="Upregulated")

plt_sens_up <- df_fgsea_up %>% filter(padj < 0.05) %>%
  filter(ctype %in% "OPC-like") %>%
  mutate("type"=ifelse(NES > 0, "resistant","sensitive")) %>%
  arrange(NES) %>%
  ggplot(., aes(x=1:length(NES),y=-(NES))) + 
  #geom_text_repel(aes(label=ifelse(NES > 1.5, condition, "")), force=10, max.overlaps =100, size=2) +
  geom_text_repel(aes(label=ifelse(NES < -1.5, condition, "")), force=100, max.overlaps =Inf, size=2) +
  geom_point(aes(color=factor(type))) +
  geom_hline(yintercept=0, color="grey") + 
  #facet_wrap(~ , scales = "free_y")+
  labs(title="", x ="rank", y = "sensitivity score") +
  ylim(-3, 3) +
  theme_classic()

plot(plt_sens_up)

ggsave("/Volumes/iorio/Raffaele/Project_ADD3/review/sensitivity_OPClike.pdf",plt_sens_up, height = 5)

test <- df_fgsea_up %>% filter(padj < 0.05) %>%
  filter(ctype=="OPC-like") %>%
  mutate("type"=ifelse(NES > 0, "resistant","sensitive")) %>%
  filter(type=="sensitive") %>%
  arrange(NES) %>%
  dplyr::select(condition) %>% pull() %>% unique()


NES_matrix_up <- df_fgsea_up %>% 
  dplyr::select(id,condition,NES) %>% 
  pivot_wider(names_from=id, values_from=NES, values_fn = median)

DRUGS <- NES_matrix_up$condition
CL_ID <- colnames(NES_matrix_up)[-1]
NES_matrix_up$condition <- NULL
NES_matrix_up <- t(NES_matrix_up) 
colnames(NES_matrix_up) <- DRUGS

cell_type_ann_up <- data.frame(idx=rownames(NES_matrix_up)) %>% separate(idx, into = c("cell type", "class", "sample"), sep = "_")
rownames(cell_type_ann_up) <- rownames(NES_matrix_up)
cell_type_ann_up$class <- NULL
cell_type_ann_up$sample <- NULL

dim(scale(NES_matrix_up))

set.seed(666)
pheatmap(NES_matrix_up,
         cluster_cols =T,
         cluster_rows = T,
         annotation_row = cell_type_ann_up,
         show_rownames = F, 
         show_colnames = F)
         #filename="~/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/pheatmap_neftel_sclevel_up_NES.pdf",
         #height = 6,
         #width = 5)

NES_matrix_up %>% 
  as.data.frame() %>% 
  mutate("id"=rownames(NES_matrix_up)) %>%
  separate(id, into = c("cell_type", "class", "sample"), "_") %>%
  filter(cell_type=="OPC-like") %>%
  
  
dim(NES_matrix_up)
data_scaled <- scale(NES_matrix_down)
sum(is.na(data_scaled))
library(FactoMineR)
res.pca <- PCA(NES_matrix_up, scale.unit = T, graph=F)
res.pca <- as.data.frame(res.pca$ind$coord)
res.pca$Idx <- rownames(res.pca)
res.pca <- res.pca %>% separate(Idx, into = c("ctype", "class","sample"), sep = "_")
library(scales)
hex <- hue_pal()(8)

pca_nes <- ggplot(res.pca, aes(x=Dim.1,y=Dim.2, color=factor(ctype))) +
  geom_point(size=5) +
  scale_color_manual(
    values = c("AC-like"=hex[6], 
               "NPC1-like"=hex[5],
               "NPC2-like"=hex[2],
               "OPC-like"=hex[4],
               "MES1-like"=hex[3],
               "MES2-like"=hex[8],
               "G1/S"=hex[1],
               "G2/M"=hex[7])) +
  labs(y= "PC2", x = "PC1") +
  theme_classic()

plot(pca_nes)
ggsave("~/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/pca_nes_upregulated.pdf", pca_nes,width = 7, height = 5)

##
df_fgsea_down <- df_fgsea %>% filter(class=="Downregulated")

NES_matrix_down <- df_fgsea_down %>% 
  dplyr::select(id,condition,NES) %>% 
  pivot_wider(names_from=id, values_from=NES, values_fn = median)

DRUGS <- NES_matrix_down$condition
CL_ID <- colnames(NES_matrix_down)[-1]
NES_matrix_down$condition <- NULL
NES_matrix_down <- t(NES_matrix_down) 
colnames(NES_matrix_down) <- DRUGS

cell_type_ann_down <- data.frame(idx=rownames(NES_matrix_down)) %>% separate(idx, into = c("cell type", "class", "code"), sep = "_")
rownames(cell_type_ann_down) <- rownames(NES_matrix_down)
cell_type_ann_down$class <- NULL
cell_type_ann_down$code <- NULL

set.seed(666)
pheatmap(NES_matrix_down,
         cluster_cols =T,
         cluster_rows = T,
         annotation_row = cell_type_ann_down,
         show_rownames = F, 
         show_colnames = F,
         filename="Desktop/pheatmap_baseline_sclevel_down_NES.pdf",
         height = 6,
         width = 5)


dim(NES_matrix_down)

res.pca <- PCA(NES_matrix_down, scale.unit = F, graph=F)
res.pca <- as.data.frame(res.pca$ind$coord)
res.pca$Idx <- rownames(res.pca)
res.pca <- res.pca %>% separate(Idx, into = c("ctype", "class","sample"), sep = "_")
hex <- hue_pal()(7)

pca_degs <- ggplot(res.pca, aes(x=Dim.1,y=Dim.2, color=factor(ctype))) +
  geom_point(size=5) +
  scale_color_manual(
    values = c("Astrocytes"=hex[6], 
               "Ependymal cells"=hex[1],
               "Neurons"=hex[7],
               "Oligodendrocyte progenitor cells"=hex[3],
               "Oligodendrocytes"=hex[5],
               "Retinal ganglion cells"=hex[4],
               "Serotonergic neurons"=hex[2])) +
  labs(y= "PC2", x = "PC1") +
  theme_classic()

plot(pca_degs)


plt <- lapply(names(myGO_LIST), function(tcga_code){
  df_fgsea_down %>%
    filter(TCGA_sample==tcga_code) %>%
    filter(padj < 0.05) %>%
    filter(ctype %in% "OPC-like") %>%
    mutate("type"=ifelse(NES > 0, "resistant","sensitive")) %>%
    arrange(NES) %>%
    ggplot(., aes(x=1:length(NES),y=NES)) + 
    geom_text_repel(aes(label=ifelse(NES > 1.5, condition, "")), force=10, max.overlaps =30, size=2) +
    geom_text_repel(aes(label=ifelse(NES < -1.5, condition, "")), force=10, max.overlaps =30, size=2) +
    geom_point(aes(color=factor(type))) +
    geom_hline(yintercept=0, color="grey") + 
    #facet_wrap(~ ctype, scales = "free_y")+
    labs(title=paste0(tcga_code," OPC-like Downregulated DEGs"), x ="rank", y = "sensitivity score",color="predicted") +
    theme_classic()
})

names(plt) <- names(myGO_LIST)
plt$`C3N-03186`

write.table(file = "~/Downloads/c3n-01814_drugs_up.txt", plt$`C3N-01814`$data$condition, row.names = F)

# downstream drugsea 

df_to_dsea <- do.call(
  rbind.data.frame, 
  lapply(1:ncol(fgsea_results), function(col){
    df <- do.call(rbind, fgsea_results[,col]) %>%
      as.data.frame() %>%
      filter(padj < 0.05) %>%
      mutate(TCGA_sample=colnames(fgsea_results)[col]) %>%
      mutate(id=paste0(pathway,"_",TCGA_sample)) %>%
      group_by(id) %>%
      arrange(desc(NES), .by.group=TRUE) %>%
      summarise(ranked_drugs=list(setNames(NES,condition))) %>%
      ungroup() %>%
      as.data.frame()
  })
)

df_to_dsea_list <- list()

for (i in 1:nrow(df_to_dsea)) {
  vec <- as.vector(unlist(df_to_dsea[i,2]))
  names(vec) <- names(unlist(df_to_dsea[i,2]))
  df_to_dsea_list[[i]] <- vec
}

# on single targets (slide 30)
single_targets <- unique(screened_compounds$TARGET)
single_targets <- unlist(str_split(single_targets,","))
idx <- grep('[:lower:]', single_targets)
single_targets <- single_targets[-idx]
single_targets <- single_targets[-c(301,302,303)]
idx <- grep(" ",single_targets)
single_targets <- single_targets[-idx]
na.omit(single_targets)
single_targets[443] <- "PARP5b"


library(tibble)
drug_tab <- do.call(rbind, drug_name)
drug_tab <- as.data.frame(drug_tab) %>% rownames_to_column()
colnames(drug_tab) <- c("TARGET","NAME")

drug_tab <- drug_tab %>% separate(.,col="TARGET", into =c("idx","GENE_TARGET"), sep=":")
drug_tab$idx <- NULL

gmt_target_gene <- drug_tab
rWikiPathways::writeGMT(gmt_target_gene, "~/OneDrive - Htechnopole/OD_Old/Documents/Project-Regorafenib/untitled folder/drug_target_gene.gmt")
GO_file <- "~/OneDrive - Htechnopole/OD_Old/Documents/Project-Regorafenib/untitled folder/drug_target_gene.gmt"
myGO_target_gene <- fgsea::gmtPathways(GO_file)

run_fgsea_for_all <- function(gene_set_name, stats) {
  res <- fgsea(
    pathways = myGO_target_gene,
    stats = sort(stats, decreasing = TRUE),
    minSize = 2,
    maxSize = 200,
    nperm = 1000
  )
  res$gene_set <- gene_set_name  # Add the name of the gene set to the results
  return(res)
}

# Loop over df_to_dsea_list and run fgsea for each entry
fgsea_results <- lapply(names(df_to_dsea_list), function(gene_set_name) {
  stats <- df_to_dsea_list[[gene_set_name]]
  run_fgsea_for_all(gene_set_name, stats)
})

# Combine results into a single data frame
target_gene_fgsea_results <- do.call(rbind, fgsea_results) %>% separate(gene_set, into = c("ctype", "class","sample"), sep = "_")

# test enrichment for the pathway to which the target of the drug belong
pathway_targets <- unique(screened_compounds$TARGET_PATHWAY)

drug_name <- list()
for(pt in pathway_targets){
  for (r in 1:nrow(screened_compounds)){
    if(grepl(pt, screened_compounds[r,"TARGET_PATHWAY"])){
      drug_name[[paste0(r,":",pt)]] <- screened_compounds[r,"DRUG_NAME"]
    }
  }
}

drug_tab <- do.call(rbind, drug_name)
drug_tab <- as.data.frame(drug_tab) %>% rownames_to_column()
colnames(drug_tab) <- c("TARGET_PATHWAY","NAME")
idx <- seq(from = 1, to = 1382, by = 2)
drug_tab$TARGET_PATHWAY <- unlist(str_split(drug_tab$TARGET_PATHWAY,":"))[-idx]
gmt_target_pathway <- drug_tab
rWikiPathways::writeGMT(gmt_target_pathway, "~/OneDrive - Htechnopole/OD_Old/Documents/Project-Regorafenib/untitled folder/drug_target_pathway.gmt")
GO_file <- "~/OneDrive - Htechnopole/OD_Old/Documents/Project-Regorafenib/untitled folder/drug_target_pathway.gmt"
myGO_target_pathway <- fgsea::gmtPathways(GO_file)


run_fgsea_for_all <- function(gene_set_name, stats) {
  res <- fgsea(
    pathways = myGO_target_pathway,
    stats = sort(stats, decreasing = TRUE),
    minSize = 2,
    maxSize = 200,
    nperm = 1000
  )
  res$gene_set <- gene_set_name  # Add the name of the gene set to the results
  return(res)
}

# Loop over df_to_dsea_list and run fgsea for each entry
fgsea_results <- lapply(names(df_to_dsea_list), function(gene_set_name) {
  stats <- df_to_dsea_list[[gene_set_name]]
  run_fgsea_for_all(gene_set_name, stats)
})

# Combine results into a single data frame
target_pathway_fgsea_results <- do.call(rbind, fgsea_results) %>% separate(gene_set, into = c("ctype", "class","sample"), sep = "_")

RES_PATHWAY <- target_pathway_fgsea_results %>% 
  filter(padj < 0.05) %>%
  filter(ctype == "OPC-like")


RES_TARGET <- target_gene_fgsea_results %>% 
  filter(padj < 0.05) %>%
  filter(ctype == "OPC-like")


fgsea::plotEnrichment(myGO_target_pathway$`EGFR signaling`, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-00662`) + theme_classic() + 
  labs(title="C3N-00662 EGFR signaling", x ="rank", y = "ES") 

fgsea::plotEnrichment(myGO_target_pathway$`ERK MAPK signaling`, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-00662`) + theme_classic() + 
  labs(title="C3N-00662 ERK MAPK signaling", x ="rank", y = "ES") 


fgsea::plotEnrichment(myGO_target_pathway$`EGFR signaling`, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-02783`) + theme_classic() + 
  labs(title="C3N-02783 EGFR signaling", x ="rank", y = "ES") 

fgsea::plotEnrichment(myGO_target_pathway$`EGFR signaling`, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-03188`) + theme_classic() + 
  labs(title="C3N-03188 EGFR signaling", x ="rank", y = "ES") 


###

fgsea::plotEnrichment(myGO_target_gene$EGFR, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-00662`) + theme_classic() + 
  labs(title="C3N-00662 EGFR", x ="rank", y = "ES") 

fgsea::plotEnrichment(myGO_target_gene$ERK2, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-00662`) + theme_classic() + 
  labs(title="C3N-00662 ERK2", x ="rank", y = "ES") 

fgsea::plotEnrichment(myGO_target_gene$EGFR, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-01814`) + theme_classic() + 
  labs(title="C3N-01814 EGFR", x ="rank", y = "ES") 

fgsea::plotEnrichment(myGO_target_gene$ERBB2, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-01814`) + theme_classic() + 
  labs(title="C3N-01814 ERBB2", x ="rank", y = "ES") 

fgsea::plotEnrichment(myGO_target_gene$EGFR, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-02190`) + theme_classic() + 
  labs(title="C3N-02190 EGFR", x ="rank", y = "ES")

fgsea::plotEnrichment(myGO_target_gene$BRAF, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-02190`) + theme_classic() + 
  labs(title="C3N-02190 BRAF", x ="rank", y = "ES")

fgsea::plotEnrichment(myGO_target_gene$EGFR, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-02783`) + theme_classic() + 
  labs(title="C3N-02783 EGFR", x ="rank", y = "ES")

fgsea::plotEnrichment(myGO_target_gene$ERBB2, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-02783`) + theme_classic() + 
  labs(title="C3N-02783 ERBB2", x ="rank", y = "ES") 

fgsea::plotEnrichment(myGO_target_gene$PI3K, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-02783`) + theme_classic() + 
  labs(title="C3N-02783 PI3K", x ="rank", y = "ES") 

fgsea::plotEnrichment(myGO_target_gene$EGFR, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-03188`) + theme_classic() + 
  labs(title="C3N-03188 EGFR", x ="rank", y = "ES")

fgsea::plotEnrichment(myGO_target_gene$ERBB2, stats=df_to_dsea_list$`OPC-like_Downregulated_C3N-03188`) + theme_classic() + 
  labs(title="C3N-03188 ERBB2", x ="rank", y = "ES") 

