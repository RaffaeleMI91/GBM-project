## Kaplan-Meier with TCGA BULK DATA
dirname <- '/Volumes/iorio/Raffaele/Project_ADD3/review/TCGA/Bulk/'
sample_sheet <- 'gdc_sample_sheet.2022-09-28.tsv'
dirs <- list.dirs(dirname)
idx <- grep("log",dirs)
dirs <- dirs[-idx]

files <- list()

for (i in 1:length(dirs)) {
  fn <- list.files(dirs[i])
  files[[i]] <- paste(dirs[i], fn, sep='/')
}

files <- unlist(files)
files <- grep('.tsv', files, value=TRUE)

META <- read.table(paste(dirname,sample_sheet, sep='/'), 
                   sep='\t', 
                   stringsAsFactors = FALSE, 
                   row.names = 1, 
                   header=TRUE)

strfile <- strsplit(files, split="/")
strfile <- unlist(strfile)
strfile <- grep('.tsv', strfile, value=TRUE)

META <- META[which(META$File.Name %in% strfile),]

files <- files[-1]

BULK <- lapply(
  files, function(x){
    data <- read.table(x, sep='\t', stringsAsFactors = FALSE, skip = 1, row.names = 1, header=TRUE)
    data <- data %>% filter(gene_type == 'protein_coding') %>% dplyr::select(gene_name, unstranded, tpm_unstranded)
    data <- data[complete.cases(data),]
  }
)

GENES <- Reduce('intersect', lapply(1:length(BULK), function(x) BULK[[x]][, 'gene_name']))

UNs <- Reduce(cbind, lapply(
  1:length(BULK), function(x) BULK[[x]][which(GENES %in% BULK[[x]]$gene_name), 'unstranded'])
)

TPMs <- Reduce(cbind, lapply(
  1:length(BULK), function(x) BULK[[x]][which(GENES %in% BULK[[x]]$gene_name), 'tpm_unstranded'])
)

DATA <- cbind.data.frame(GENES,TPMs)

rownames(DATA) <- GENES
DATA$GENES <- NULL
colnames(DATA) <- META$Case.ID
dim(DATA)

DATA <- t(DATA)

# Implement method from Lin et al (2024)

UP_DEF <- deg_df %>% filter(class=="Upregulated") %>% dplyr::select(cell_type, gene)
UP_DEF$Sample <- NULL
UP_DEF$class <- NULL

UP_DEF <- UP_DEF %>% filter(gene %in% colnames(DATA))

###############################

clinical <- data.table::fread("/Volumes/iorio/Raffaele/Project_ADD3/review/clinical.project-tcga-gbm.2024-07-18/clinical.tsv",sep="\t",header = T)
clinical <- as.data.frame(clinical) %>% dplyr::select(case_submitter_id, days_to_death, age_at_index, vital_status)
which(duplicated(names(clinical)))
clinical <- clinical %>% filter(case_submitter_id %in% META$Case.ID)
clinical <- (clinical[!duplicated(clinical),])
clinical$time <- as.numeric(clinical$days_to_death)
clinical$status <- ifelse(clinical$vital_status == "Dead", 1, 0)
clinical <- clinical %>% dplyr::select(case_submitter_id,vital_status,time,status)
DATA_subset <- DATA[intersect(clinical$case_submitter_id, rownames(DATA)), ]

MERGED_DATA <- cbind.data.frame(clinical, DATA_subset)

library(survival)
library(survminer)

results <- list()
significant_genes <- data.frame(Gene = character(), PValue = numeric())

opc_upreg_ess[-c(30,33,135)]
opc_downreg_ess[-c(668)]
intersect(labels, opc_downreg_ess)

SET <- Reduce(union, ess_up_list)

for (gene in SET[-125]) {
  print(grep(gene,SET))
  # Check if the gene exists in the MERGED_DATA columns
  if (!(gene %in% colnames(MERGED_DATA))) {
    next  # Skip the gene if it is not in the dataset
  }
  MERGED_DATA$expression_group <- ifelse(MERGED_DATA[, gene] > median(MERGED_DATA[, gene]), "High", "Low")
  # Check if the gene results in only one group (either High or Low)
  if (length(unique(MERGED_DATA$expression_group)) == 1) {
    next  # Skip this gene if there is only one group
  }
  surv_object <- Surv(time = MERGED_DATA$time, event = MERGED_DATA$status)
  fit <- survfit(surv_object ~ expression_group, data = MERGED_DATA)
  surv_test <- survdiff(surv_object ~ expression_group, data = MERGED_DATA)
  p_value <- 1 - pchisq(surv_test$chisq, length(surv_test$n) - 1)
  # Check if the result is significant
  if (p_value < 0.05) {  # Adjust threshold as needed
    # Save results for significant genes
    significant_genes <- rbind(significant_genes, data.frame(Gene = gene, PValue = p_value))
    
    # Save the plot to the results list
    plot <- ggsurvplot(fit,
                       data = MERGED_DATA,
                       pval = TRUE,
                       risk.table = TRUE,
                       title = paste("Kaplan-Meier Curve for", gene),
                       xlab = "Time",
                       ylab = "Survival Probability")
    results[[gene]] <- plot
    
    # Optionally save plots to files
    ggsave(filename = paste0("~/Downloads/KM_Curve_UpregEss_final", gene, ".png"), plot = plot$plot, dpi = 300, width = 5, height = 6)
  }
}


opc_upreg_ess_SignSurvival <- c("CCS", "GPR3", "MAP2K4", "RHOQ", "SKI", "WLS","ZNF404")
opc_downreg_ess_SignSurvival <- setdiff(names(results), opc_upreg_ess_SignSurvival)

# Check the result
