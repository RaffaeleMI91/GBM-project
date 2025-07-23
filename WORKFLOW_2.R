psbulk_degs_df <- read.csv(
  file="/Volumes/iorio/Raffaele/Project_ADD3/DEGs_PSEUDOBULKS_CTYPE_TOTAL.csv",
  stringsAsFactors = FALSE, 
  header=TRUE, 
  row.names = 1)

psbulk_degs_df <- psbulk_degs_df %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 2)

LIST_PB_DEGs <- lapply(unique(psbulk_degs_df$target_ct), function(ct){
  deg_by_ct <- psbulk_degs_df %>% filter(target_ct==ct) %>% dplyr::select(level_1) %>% pull()
  my_gmt <- NULL
  my_gmt <- rbind(my_gmt, data.frame(
    CT=c(rep(ct, length(deg_by_ct))),
    NAME=deg_by_ct,
    stringsAsFactors = FALSE))
  return(my_gmt)
})

names(LIST_PB_DEGs) <- c("AC-like","G1-S","G2-M","MES1-like","MES2-like","NPC1-like","NPC2-like","OPC-like")

lapply(names(LIST_PB_DEGs), function(x){
  rWikiPathways::writeGMT(
    LIST_PB_DEGs[[x]], 
    paste0("/Users/raffaele.iannuzzi/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/",x,"_degs.gmt"))
})

file_list <- list.files("/Users/raffaele.iannuzzi/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/")
idx <- grep("_degs.gmt",file_list)
gmt_list <- file_list[idx]

myGO_LIST <- lapply(gmt_list, function(x) {
  GO_file <- paste0("/Users/raffaele.iannuzzi/OneDrive - Htechnopole/OD_Old/Documents/Project-ADD3/",x)
  myGO <- fgsea::gmtPathways(GO_file)
  return(myGO)
})

names(myGO_LIST) <- names(LIST_PB_DEGs)

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

run_fgsea <- function(drug_name, rank, query_sets, minSize=2, maxSize=3000, nperm=1000) {
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
    myGO_LIST,function(x){
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
      mutate(PB_CTYPE=colnames(fgsea_results)[col])
    return(df)
  }
  ))

NES_matrix <- df_fgsea %>% dplyr::select(PB_CTYPE,condition,NES) %>%
  pivot_wider(names_from = condition, values_from = NES, values_fn = median)

cnames <- NES_matrix$PB_CTYPE
NES_matrix$PB_CTYPE <- NULL
NES_matrix <- t(NES_matrix)
colnames(NES_matrix) <- cnames

cell_type_ann <- data.frame(ctype=colnames(NES_matrix))
rownames(cell_type_ann) <- cell_type_ann$ctype

set.seed(666)
pheatmap(NES_matrix,
         cluster_cols=T,
         cluster_rows=T,
         annotation_col=cell_type_ann,
         show_rownames=F, 
         show_colnames=F)


# UTILS for eventually do t-test when you have more samples

pVWilc_plots <- lapply(unique(df_fgsea$condition), function(x){ 
  res <- df_fgsea %>% filter(condition==x, class=="Upregulated") %>% 
    ggboxplot(., x="ctype", y="NES", color="ctype", short.panel.labs = T) +
    rotate_x_text(angle=90) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    stat_compare_means(
      label = "p.signif", 
      method = "wilcox.test",
      method.args=list(alternative="two.sided"),
      ref.group = NULL, 
      mu=0)
  output <- ggplot_build(res)
  return(output)
})

names(pVWilc_plots) <- unique(df_fgsea$condition)

pVWilc_plots$`PLX-4720`

pVWilc <- lapply(unique(df_fgsea$condition), function(x) {
  filtered_data <- df_fgsea %>%
    filter(condition == x, class == "Upregulated") %>%
    mutate(ctype = as.factor(ctype))
  if (length(levels(filtered_data$ctype)) > 2) {
    res <- filtered_data %>% 
      group_by(ctype) %>%
      rstatix::wilcox_test(NES~1, alternative="two.sided", mu=0, paired=FALSE, ref.group = NULL)
  } else {
    res <- NULL
  }
  return(res)
})

names(pVWilc) <- unique(df_fgsea$condition)
pVWilc_filtered <- pVWilc[!sapply(pVWilc, is.null)]

eSWilc <- lapply(unique(df_fgsea$condition), function(x) {
  filtered_data <- df_fgsea %>%
    filter(condition==x, class=="Upregulated") %>%
    mutate(ctype = as.factor(ctype))
  if (length(levels(filtered_data$ctype)) > 2) {
    res <- filtered_data %>% group_by(ctype) %>%
      rstatix::wilcox_effsize(NES~1, mu=0, alternative="two.sided", paired=FALSE, ref.group = NULL) %>%
      #rowwise() %>%  # Ensures per-row calculation
      mutate(
        # Compute direction by comparing medians dynamically for each group2
        direction = ifelse(
          median(filtered_data$NES[filtered_data$ctype == ctype], na.rm = TRUE) > 0,
          1, -1
        )
      ) %>%
      ungroup()  # Ungroup after rowwise operations
  } else {
    res <- NULL
  }
  return(res)
})

names(eSWilc) <- unique(df_fgsea$condition)
eSWilc_filtered <- eSWilc[!sapply(eSWilc, is.null)]

eSWilc_data <- do.call(
  rbind, lapply(
    names(eSWilc_filtered), function(x){
      eSWilc <- eSWilc_filtered[[x]] %>% mutate(condition = x)
      return(eSWilc)
    }
  )
)

pVWilc_data <- do.call(
  rbind, lapply(
    names(pVWilc_filtered), function(x){
      pVWilc <- pVWilc_filtered[[x]] %>% mutate(condition = x)
      return(pVWilc)
    }
  )
)

combined_res <- merge(pVWilc_data,eSWilc_data) %>% select(group1, group2, condition, p, effsize, magnitude,direction)
combined_res <- combined_res %>%
  mutate(
    signed_eS = direction * effsize,
    adjusted_pv = p.adjust(p, method = "BH"),
    is_significant = p < 0.05,
    class=(
      case_when(
        is_significant & signed_eS > 0 ~ "resistant",
        is_significant & signed_eS < 0 ~ "sensitive"
      )),
    neg_log10_p = -log10(p)
  )

# Volcano plot
ggplot(combined_res, aes(x=signed_eS, y=neg_log10_p, color=class)) +
  geom_point(alpha = 0.6) +
  #scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(
    data = combined_res %>% filter(is_significant),
    aes(label = condition),
    size = 2, max.overlaps = 50       # Avoid label clutter
  ) +
  labs(
    title = "",
    x = "signed ES",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal()

