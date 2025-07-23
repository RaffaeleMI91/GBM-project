library(stringr)

Sanger_expr <- read.csv("~/Downloads/rnaseq_all_20220624/rnaseq_tpm_20220624.csv", stringsAsFactors = FALSE, header = TRUE)
tpm <- Sanger_expr %>% dplyr::filter(model_id != "SIDG06476" & model_id != "SIDG03156")
genes <- tpm$X[5:length(tpm$X)]
tpm <- tpm[,(tpm[2,] == "Sanger & Broad Cell Lines RNASeq")]
colnames(tpm) <- tpm[1,]
tpm <- tpm[-c(1:4),]
rownames(tpm) <- genes

#load("/Volumes/iorio/Raffaele/depFC/Sanger_Broad_higQ_depFC.RData")
#depFCv2 <-read.csv("/Volumes/iorio/Raffaele/depFC/filtered_scaled_DepMat_2024.csv")
write.csv(scaled_depFC, "/Volumes/iorio/Raffaele/depFC/Sanger_Broad_higQ_depFC.csv")
CEG <- read.csv("~/OneDrive - Htechnopole/OD_Old/Documents/Project-Regorafenib/data/AchillesCommonEssentialControls.csv")
IDX <- grep("\\(",unlist(strsplit(CEG$Gene," ")))
CEG <- unlist(strsplit(CEG$Gene," "))[-IDX]

dim(scaled_depFC)
scaled_depFC <- scaled_depFC[-which(rownames(scaled_depFC) %in% CEG),]
dim(scaled_depFC)

## deriving binary essentiality scores
bdep <- apply(scaled_depFC, 2, function(x){
  x[which(x >= -0.5)] <- 0
  x[which(x < -0.5)] <- 1
  x
})

NESS <- which(rowSums(bdep) == 0)
scaled_depFC <- scaled_depFC[-NESS,] # remove never essential genes
dim(scaled_depFC)

cmatch <- intersect(colnames(scaled_depFC),colnames(tpm))

tpm <- tpm[,cmatch]
SD_GEX <- apply(tpm, 1, sd)
threshold <- quantile(SD_GEX, 0.35, na.rm=T)
HVG <- names(SD_GEX[SD_GEX >= threshold])
tpm <- tpm[which(rownames(tpm) %in% HVG),]
genes <- rownames(tpm) 
tpm <- apply(tpm,2,as.numeric)
rownames(tpm) <- genes
dim(tpm)

gmatch <- intersect(rownames(scaled_depFC),rownames(tpm))

scaled_depFC <- t(scaled_depFC[gmatch, cmatch])
tpm <- t(tpm[gmatch, cmatch])

dim(tpm)
dim(scaled_depFC)

cor_dep_matrix <- cor(tpm,scaled_depFC)


#save(cor_dep_matrix, file="/Volumes/iorio/Raffaele/Project_ADD3/review/ESS_EXPR_CORR_MAT.RData")

load("OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/ESS_EXPR_CORR_MAT.RData")

sorted_list <- lapply(colnames(cor_dep_matrix),function(i){
  ref <- cor_dep_matrix[,i]
  idx <- order(ref,decreasing = T)
  vec <- ref[idx]
  names(vec) <- rownames(cor_dep_matrix)[idx]
  return(vec)
})

names(sorted_list) <- colnames(cor_dep_matrix)

sorted_list <- lapply(sorted_list, function(x) {x[is.finite(x)]})

bin_list <- list()
for (seed in names(sorted_list)){
  vec <- sorted_list[[seed]]
  vec <- as.numeric(names(vec)==seed)
  bin_list[[seed]] <- vec
}

bin_mat <- as.matrix(do.call(rbind, bin_list))
rownames(bin_mat) <- names(sorted_list)
colSums(bin_mat)

num_cores <- detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)

run_fgsea <- function(gene_name, rank, query_sets, minSize=2, maxSize=300, nperm=1000) {
  result <- fgsea(
    pathways=query_sets, 
    stats=rank, 
    minSize=minSize, 
    maxSize=maxSize,
    nperm=nperm)
  result$condition <- gene_name 
  return(result)
}

fgsea_results <- foreach(i=names(sorted_list), .combine = rbind, .packages = "fgsea") %dopar% {
  ranked_list <- sorted_list[[i]]
  res_list <- lapply(
    myGO_LIST,function(x){
      run_fgsea(
        gene_name=i, 
        rank=ranked_list, 
        query_sets=x)})
  return(res_list)
}

stopCluster(cl)

save(fgsea_results, file="/Volumes/iorio/Raffaele/Project_ADD3/review/DF_FGSEA_ESS_EXPR.RData")

load("/Volumes/iorio/Raffaele/Project_ADD3/review/data/DF_FGSEA_ESS_EXPR.RData")

df_fgsea <- do.call(rbind.data.frame, 
                    lapply(1:ncol(fgsea_results), function(col){
                      df <- do.call(rbind, fgsea_results[,col]) %>%
                        as.data.frame() %>% 
                        #filter(padj < 0.01) %>%
                        mutate(NES=ifelse(padj > 0.05, 0, NES)) %>% # same as in DREEP
                        mutate(TCGA_sample=colnames(fgsea_results)[col]) %>%
                        mutate(id=paste0(pathway,"_",TCGA_sample))
                      return(df)
                    }
                    ))


NES_dep_matrix <- df_fgsea %>% 
  dplyr::select(id,condition,NES) %>% 
  pivot_wider(names_from=id, values_from=NES)

RGENES <- NES_dep_matrix$condition
CTYPE <- colnames(NES_dep_matrix)[-1]
NES_dep_matrix$condition <- NULL
NES_dep_matrix <- t(NES_dep_matrix) 
colnames(NES_dep_matrix) <- RGENES

IDX <- grep("C3N",unlist(strsplit(rownames(NES_dep_matrix),"_")))

cell_type_ann <- data.frame(
  sample=unlist(strsplit(rownames(NES_dep_matrix),"_"))[IDX],
  ctype=unlist(strsplit(rownames(NES_dep_matrix),"_"))[-IDX])

rownames(cell_type_ann) <- rownames(NES_dep_matrix)

newCols <- colorRampPalette(grDevices::rainbow(length(unique(cell_type_ann$ctype))))
annoCol <- newCols(length(unique(cell_type_ann$ctype)))
names(annoCol) <- unique(cell_type_ann$ctype)
annoCol <- list(category = annoCol)

set.seed(666)
pheatmap(NES_dep_matrix,
         cluster_cols=T,
         cutree_rows=4,
         color=colorRampPalette(c("blue", "white", "red"))(50),
         #color=colorRampPalette(c("lightblue", "blue", "darkblue"))(50),
         annotation_colors=annoCol,
         annotation_row=cell_type_ann, 
         border_color = "black",
         show_rownames = F, 
         show_colnames = F)#,
         #filename="/Volumes/iorio/Raffaele/Project_ADD3/review/pheatmap_baseline_dep_upreg.pdf")

BiocManager::install("msigdbr")
library(msigdbr)
library(purrr)

pathways <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

test <- do.call(
  rbind.data.frame, 
  lapply(1:ncol(fgsea_results), function(col){
    df <- do.call(rbind, fgsea_results[,col]) %>%
      as.data.frame() %>% 
      filter(padj < 0.05) %>%
      mutate(TCGA_sample=colnames(fgsea_results)[col]) %>%
      mutate(id=paste0(pathway,"_",TCGA_sample)) %>%
      group_by(id) %>%
      arrange(desc(NES), .by.group=TRUE) %>%
      summarise(ranked_genes=list(setNames(NES,condition))) %>%
      ungroup()
    })
)

gsea_results <- test %>%
  mutate(
    fgsea_result = map(ranked_genes, ~ fgsea(
      pathways = pathways,
      stats = .x,
      minSize = 2,
      maxSize = 300,
      nperm = 10000
    ))
  ) %>% unnest(fgsea_result) 


res <- lapply(unique(gsea_results$id), function(ct){
  gsea_results %>% filter(id==ct) %>% 
    mutate(significance = ifelse(pval < 0.05, "pval<0.05", "ns")) %>%
    ggplot(., aes(x = reorder(pathway, NES), y = NES, fill = significance)) +
    geom_bar(stat="identity", width=0.7) +  # Map NES to size and padj to color
    scale_color_manual(values = c("pval<0.05" = "red", "ns" = "blue")) + 
    coord_flip() +
    labs(
      title = "",
      x = "Msigdbr pathway",
      y = "NES"
      ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
  )
})

names(res) <-   unique(gsea_results$id)
res$`Oligodendrocytes_C3N-01814`


