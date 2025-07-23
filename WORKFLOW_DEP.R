getwd()
library(parallel)
library(doParallel)
library(CoRe)
?CoRe.scale_to_essentials()
depMap <- read.csv("~/Downloads/CRISPRGeneEffect_24Q4.csv", stringsAsFactors = FALSE, header = TRUE, row.names=1)
genes <- lapply(strsplit(colnames(depMap), "\\.."), function(x) x[[1]])
colnames(depMap) <- unlist(genes)
depMap <- t(depMap)

dim(depMap)

data('curated_BAGEL_essential')
data('curated_BAGEL_nonEssential')

data("EssGenes.DNA_REPLICATION_cons")
data("EssGenes.HISTONES")
data("EssGenes.KEGG_rna_polymerase")
data("EssGenes.PROTEASOME_cons")
data("EssGenes.ribosomalProteins")
data("EssGenes.SPLICEOSOME_cons")

scaled_depMap <- CoRe.scale_to_essentials(depMap, curated_BAGEL_essential, curated_BAGEL_nonEssential)

CEG <- read.csv("~/OneDrive - Htechnopole/OD_Old/Documents/Project-Regorafenib/data/AchillesCommonEssentialControls.csv")
IDX <- grep("\\(",unlist(strsplit(CEG$Gene," ")))
CEG <- unlist(strsplit(CEG$Gene," "))[-IDX]

coreFitGenes <- Reduce(union, list(
  EssGenes.DNA_REPLICATION_cons,
  EssGenes.HISTONES,
  EssGenes.KEGG_rna_polymerase,
  EssGenes.PROTEASOME_cons,
  EssGenes.ribosomalProteins,
  EssGenes.SPLICEOSOME_cons,
  curated_BAGEL_essential,
  CEG))

length(coreFitGenes)

scaled_depMap <- scaled_depMap[-which(rownames(scaled_depMap) %in% coreFitGenes),]

binaryDep <- apply(scaled_depMap, 2, function(x){
  x[which(x >= -0.5)] <- 0
  x[which(x < -0.5)] <- 1
  x
}
)

NESS <- which(rowSums(binaryDep) == 0)
scaled_depMap <- scaled_depMap[-NESS,] # remove never essential genes
dim(scaled_depMap)

gex <- read.csv("~/Downloads/Batch_corrected_Expression_Public_24Q4_subsetted.csv", 
                stringsAsFactors = FALSE, 
                header = TRUE, 
                row.names=1)
sd <- apply(gex, 2, sd)
th <- quantile(sd, 0.35, na.rm=T)
hvg <- names(sd[sd >= th])
gex <- gex[,which(colnames(gex) %in% hvg)]
gex <- t(gex)


col_match <- intersect(colnames(scaled_depMap),colnames(gex))
row_match <- intersect(rownames(scaled_depMap),rownames(gex))

scaled_depMap <- scaled_depMap[row_match, col_match]
gex <- gex[row_match, col_match]

results <- mapply(function(x, y) 
  cor.test(x, y, method = "spearman", complete.obs=T),
  split(scaled_depMap, row(scaled_depMap)), 
  split(gex, row(gex)))

results <- as.data.frame(t(results))
rownames(results) <- row_match

filtered <- results %>% 
  select(p.value, estimate) %>%
  filter(p.value < 0.05, abs(as.numeric(estimate)) > 0.1)

filtered$p.value <- as.numeric(filtered$p.value)
filtered$estimate <- as.numeric(filtered$estimate)

filtered_scaled_depMap <- scaled_depMap[rownames(filtered),]
filtered_gex <- gex[rownames(filtered),]

dim(filtered_scaled_depMap)

dim(gex)
dim(scaled_depMap)
cor_dep_matrix <- cor(t(gex), t(scaled_depMap))

sorted_list <- lapply(colnames(cor_dep_matrix), function(i){
  ref <- cor_dep_matrix[,i]
  idx <- order(ref,decreasing = T)
  vec <- ref[idx]
  names(vec) <- rownames(cor_dep_matrix)[idx]
  return(vec)
})

names(sorted_list) <- colnames(cor_dep_matrix)

sorted_list <- lapply(sorted_list, function(x) {x[is.finite(x)]})

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

df_fgsea <- do.call(rbind.data.frame, 
                    lapply(1:ncol(fgsea_results), function(col){
                      df <- do.call(rbind, fgsea_results[,col]) %>%
                        as.data.frame() %>% 
                        #filter(padj < 0.01) %>%
                        mutate(NES=ifelse(padj > 0.05, 0, NES)) %>% # same as in DREEP
                        mutate(TCGA_sample=colnames(fgsea_results)[col]) %>%
                        mutate(id=paste0(pathway,"_",TCGA_sample)) %>%
                        separate(pathway, into=c("ctype","class"),sep = "_") %>%
                      return(df)
                    }
                    ))

df_fgsea_up <- df_fgsea %>% filter(class=="Upregulated")
df_fgsea_down <- df_fgsea %>% filter(class=="Downregulated")


NES_dep_matrix <- df_fgsea_down %>% 
  dplyr::select(id,condition,NES) %>% 
  pivot_wider(names_from=id, values_from=NES)

RGENES <- NES_dep_matrix$condition
CTYPE <- colnames(NES_dep_matrix)[-1]
NES_dep_matrix$condition <- NULL
NES_dep_matrix <- t(NES_dep_matrix) 
colnames(NES_dep_matrix) <- RGENES

cell_type_ann <- data.frame(idx=rownames(NES_dep_matrix)) %>% separate(idx, into = c("cell type", "class", "sample"), sep = "_")

rownames(cell_type_ann) <- rownames(NES_dep_matrix)

newCols <- colorRampPalette(grDevices::rainbow(length(unique(cell_type_ann$`cell type`))))
annoCol <- newCols(length(unique(cell_type_ann$ctype)))
names(annoCol) <- unique(cell_type_ann$ctype)
annoCol <- list(category = annoCol)
cell_type_ann <- cell_type_ann %>% dplyr::select(`cell type`)

set.seed(666)
pheatmap(NES_dep_matrix,
         cluster_cols=T,
         color=colorRampPalette(c("blue", "white", "red"))(50),
         #color=colorRampPalette(c("lightblue", "blue", "darkblue"))(50),
         annotation_colors=annoCol,
         annotation_row=cell_type_ann, 
         border_color = "black",
         show_rownames = F, 
         show_colnames = F)
         #filename="~/Downloads/pheatmap_neftel_sclevel_down_essential_NES.pdf",
         #height = 6,
         #width = 5)
dev.off()


# PLOT RANKED NES BY PATIENTS

plt <- lapply(names(myGO_LIST), function(tcga_code){
  df_fgsea_down %>%
    filter(TCGA_sample==tcga_code) %>%
    filter(padj < 0.05) %>%
    filter(ctype %in% "OPC-like") %>%
    mutate("type"=ifelse(NES < 0, "essential","non essential")) %>%
    arrange(NES) %>%
    ggplot(., aes(x=1:length(NES),y=NES)) + 
    geom_text_repel(aes(label=ifelse(NES < -1.5, condition, "")), force=10, max.overlaps =50, size=4) +
    geom_point(aes(color=factor(type))) +
    geom_hline(yintercept=0, color="grey") + 
    labs(title=paste0(tcga_code," OPC-like Downregulated DEGs"), x ="rank", y = "essentiality score", color="predicted") +
    ylim(-3, 3) +
    theme(axis.text.x = element_text(size = 14),  # Increase font size of x-axis labels
          axis.text.y = element_text(size = 14)) +
    theme_classic() +
    scale_color_manual(values = c("essential" = "blue", "non essential" = "grey"))
})

names(plt) <- names(myGO_LIST)

plt$`C3N-01814`
plt$`C3N-00662`
plt$`C3N-02190`
plt$`C3N-02783`
plt$`C3N-03186`
plt$`C3N-03188`

ess_up_list <- lapply(names(myGO_LIST), function(tcga_code){
  res <- df_fgsea_up %>%
    filter(TCGA_sample==tcga_code) %>%
    filter(padj < 0.05, NES <0) %>%
    filter(ctype %in% "OPC-like") %>%
    dplyr::select(condition) %>%
    pull() %>% 
    unique()
  return(res)
})

ess_down_list <- lapply(names(myGO_LIST), function(tcga_code){
  res <- df_fgsea_down %>%
    filter(TCGA_sample==tcga_code) %>%
    filter(padj < 0.05, NES <0) %>%
    filter(ctype %in% "OPC-like") %>%
    dplyr::select(condition) %>%
    pull() %>% 
    unique()
  return(res)
})


names(ess_down_list) <- names(ess_up_list) <- names(myGO_LIST)

intersect(Reduce(union, ess_down_list), labels)
intersect(Reduce(union, ess_up_list),labels)

write.table(Reduce(union, ess_up_list),  "~/Downloads/ess_up_union.txt", row.names=F)
write.table(Reduce(union, ess_down_list),  "~/Downloads/ess_down_union.txt", row.names=F)


write.table(ess_down_list$`C3N-01814`, "Downloads/c3n01814_opc_like_down_ess.txt", row.names=F)
write.table(ess_down_list$`C3N-00662`, "Downloads/c3n00662_opc_like_down_ess.txt", row.names=F)
write.table(ess_down_list$`C3N-02190`, "Downloads/c3n02190_opc_like_down_ess.txt", row.names=F)
write.table(ess_down_list$`C3N-02783`, "Downloads/c3n02783_opc_like_down_ess.txt", row.names=F)
write.table(ess_down_list$`C3N-03186`, "Downloads/c3n03186_opc_like_down_ess.txt", row.names=F)
write.table(ess_down_list$`C3N-03188`, "Downloads/c3n03188_opc_like_down_ess.txt", row.names=F)

write.table(ess_up_list$`C3N-01814`, "Downloads/c3n01814_opc_like_up_ess.txt", row.names=F)
write.table(ess_up_list$`C3N-00662`, "Downloads/c3n00662_opc_like_up_ess.txt", row.names=F)
write.table(ess_up_list$`C3N-02190`, "Downloads/c3n02190_opc_like_up_ess.txt", row.names=F)
write.table(ess_up_list$`C3N-02783`, "Downloads/c3n02783_opc_like_up_ess.txt", row.names=F)
write.table(ess_up_list$`C3N-03186`, "Downloads/c3n03186_opc_like_up_ess.txt", row.names=F)
write.table(ess_up_list$`C3N-03188`, "Downloads/c3n03188_opc_like_up_ess.txt", row.names=F)


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
      arrange((NES), .by.group=TRUE) %>%
      summarise(ranked_genes=list(setNames(NES,condition))) %>%
      ungroup()
  })
)


ids <- c("OPC-like_Downregulated_C3N-00662",
         "OPC-like_Upregulated_C3N-00662",
         "OPC-like_Downregulated_C3N-02190",
         "OPC-like_Upregulated_C3N-02190",
         "OPC-like_Downregulated_C3N-01814",
         "OPC-like_Upregulated_C3N-01814",
         "OPC-like_Downregulated_C3N-03188",
         "OPC-like_Upregulated_C3N-03188",
         "OPC-like_Downregulated_C3N-03186",
         "OPC-like_Upregulated_C3N-03186",
         "OPC-like_Downregulated_C3N-02783",
         "OPC-like_Upregulated_C3N-02783")

test <- test %>% filter(id %in% ids)

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
      title = ct,
      x = "Msigdbr pathway",
      y = "NES"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=20),
      axis.text.y = element_text(size = 15)
    )
})
names(res) <- ids
res$`OPC-like_Upregulated_C3N-00662`

ggsave("~/Downloads/patwhay_up_c3n00662.png",res$`OPC-like_Upregulated_C3N-00662`,width=10)
ggsave("~/Downloads/patwhay_up_c3n02190.png",res$`OPC-like_Upregulated_C3N-02190`,width=10)
ggsave("~/Downloads/patwhay_up_c3n01814.png",res$`OPC-like_Upregulated_C3N-01814`,width=10)
ggsave("~/Downloads/patwhay_up_c3n03188.png",res$`OPC-like_Upregulated_C3N-03188`,width=10)
ggsave("~/Downloads/patwhay_up_c3n03186.png",res$`OPC-like_Upregulated_C3N-03186`,width=10)
ggsave("~/Downloads/patwhay_up_c3n02783.png",res$`OPC-like_Upregulated_C3N-02783`,width=10)

ggsave("~/Downloads/patwhay_down_c3n00662.png",res$`OPC-like_Downregulated_C3N-00662`,width=10)
ggsave("~/Downloads/patwhay_down_c3n02190.png",res$`OPC-like_Downregulated_C3N-02190`,width=10)
ggsave("~/Downloads/patwhay_down_c3n01814.png",res$`OPC-like_Downregulated_C3N-01814`,width=10)
ggsave("~/Downloads/patwhay_down_c3n03188.png",res$`OPC-like_Downregulated_C3N-03188`,width=10)
ggsave("~/Downloads/patwhay_down_c3n03186.png",res$`OPC-like_Downregulated_C3N-03186`,width=10)
ggsave("~/Downloads/patwhay_down_c3n02783.png",res$`OPC-like_Downregulated_C3N-02783`,width=10)

plt_ess_down <- df_fgsea_down %>% filter(ctype == "OPC-like") %>% filter(padj < 0.05) %>%
  mutate("type"=ifelse(NES < 0, "essential","non essential")) %>%
  arrange(NES) %>%
  ggplot(., aes(x=1:length(NES),y=NES)) + 
  #geom_text_repel(aes(label=ifelse(NES > 1.5, condition, "")), force=10, max.overlaps =100, size=2) +
  geom_text_repel(aes(label=ifelse(NES < -1.5, condition, "")), force=10, max.overlaps =50, size=2) +
  geom_point(aes(color=factor(type))) +
  geom_hline(yintercept=0, color="grey") + 
  #facet_wrap(~ , scales = "free_y")+
  labs(title="OPC-like Downregulated DEGs", x ="rank", y = "essentiality score") +
  ylim(-3, 3) +
  theme_classic() +
  scale_color_manual(values = c("essential" = "blue", "non essential" = "grey"))
dev.off()

plot(plt_ess_down)

ggsave("~/Downloads/essentiality_OPClike_down.pdf",plt_ess_down, height = 5)


library(clusterProfiler)
library(enrichplot) 
library(org.Hs.eg.db)
library(topGO)

BPmapping <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol") 
genesUniverse <- unique(unlist(BPmapping))

GO <- enrichGO(gene = opc_upreg_ess, 
               keyType = "SYMBOL", 
               universe = genesUniverse, 
               ont="BP", 
               OrgDb = "org.Hs.eg.db")

#pdf(paste(resultsDir,'/FIG5_E.pdf', sep=''),8,6)

dotplot(GO, showCategory=10)
dev.off()

library(openxlsx)

genes <- read.xlsx("~/Downloads/oligoscore.xlsx", sheet = 1)
write.table(toupper(na.omit(genes$Gene)), "~/Downloads/oligo_genes.txt", row.names=F)

intersect(Reduce(union, ess_down_list), labels)
intersect(Reduce(union, ess_up_list),labels)

consensus <- Reduce(union, ess_down_list)
background <- genesUniverse
labels <- read.table("~/Downloads/oligo_genes.txt")
labels <- labels$V1
labels <- labels[-c(1,length(labels))]

TP <- length(intersect(consensus, labels)) 
FP <- length(intersect(consensus, setdiff(background, labels))) 
FN <- length(intersect(setdiff(background, consensus), labels)) 
TN <- length(intersect(setdiff(background, consensus),
                       setdiff(background,  labels))) 

dat <- data.frame(
  "col_1" = c(TP, FP),
  "col_2" = c(FN, TN),
  row.names = c("Consensus", "Other"),
  stringsAsFactors = FALSE
)

colnames(dat) <- c("is", "is not")

pval <- fisher.test(dat)$p.value


# data from survival analysis:
opc_downreg_ess_SignSurvival
opc_upreg_ess_SignSurvival
