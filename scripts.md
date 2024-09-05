Project Path: /Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts

Source Tree:

```
Scripts
├── ExomeClassifier
│   ├── ExomeHRDClassification
│   │   ├── SMC_HRDclassification.R
│   │   ├── TCGA_geneCNAenrichment.R
│   │   ├── decoupleR_full.R
│   │   ├── TCGA_chromosomeArmAlterations_HRDonly.R
│   │   ├── HRD_Hypoxia.R
│   │   ├── TCGA_chromosomeArmAlterations.R
│   │   ├── dNdS_Analysis.R
│   │   ├── TCGA_HRDclassification.R
│   │   ├── TCGA_HRDhallmarks.R
│   │   ├── TCGA_geneCNAenrichment_HRDonly.R
│   │   ├── TCGA_compareHRDclassifiers.R
│   │   └── TCGA_HRDclassification_thresholds.R
│   └── ICGC_PhenotypeClusteringAndSimulationAnalysis
│       ├── Likelihood_Cluster_Generation.R
│       ├── ICGC_deconstructSigs_genome2exome.R
│       ├── ICGC_BRCA_sigProfilerContributions.R
│       ├── ICGC_simulations_indelProportions.R
│       ├── ICGC_BRCA_UKandEU_MatrixGeneration.R
│       ├── ICGC_simulations_indelLikelihoodAlterations.R
│       ├── ICGC_simulations_clusterReassignment.R
│       ├── ICGC_PhenotypeClusterDevelopment_normalised.R
│       └── ICGC_indelCounts_ExomeVsGenome.R
├── TCGA_HRDclassificationHeatmap.jpg
├── README.md
└── TranscriptionalSignature
    ├── TranscriptionalSignatureDevelopment
    │   ├── RNAseqTransformationScript.R
    │   ├── MultinomialWEN_alpha0.25_1to250.R
    │   ├── ReduxSig_adjacencyMatrixFormation.R
    │   ├── BayesPrism_TumorPurityComparison.R
    │   ├── CollateAlternativeSignatures.R
    │   ├── Qian2020_signatureDropOuts.R
    │   ├── CentroidModelFormation.R
    │   ├── MultinomialElasticNet_alpha0.25_1to100.R
    │   └── TCGA_BRCA.RNASeq_prep.R
    ├── SingleCellAnalysis
    │   ├── Bassez2021_preprocessing.R
    │   ├── Qian2020_jointComparison.R
    │   ├── Bassez2021_analysis.R
    │   ├── Bassez2021_HRDprofiling.R
    │   ├── Qian2020_analysis.R
    │   ├── Chung2017_analysis.R
    │   ├── Qian2020_preprocessing.R
    │   ├── Qian2020_HRDprofiling.R
    │   ├── cpdb_interactionCounts.R
    │   └── cpdb_individualInteractions.R
    └── SignatureValidation
        ├── ISPY2_HRDscoring.R
        ├── CCLE_jointComparisons.R
        ├── GSEA_enrichR.R
        ├── TCGA_HRDscoreByERstatus.R
        ├── GSEA_pathfindR.R
        ├── TCGA_testScoring.R
        ├── SMC_validation.R
        ├── SMC_signatureComparisons.R
        ├── TCGA_trainHeatmap.R
        └── TCGA_testScoring_signatureComparisons.R

```

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/SMC_HRDclassification.R`:

```````R
#####
## Apply HRD exome classifier to SMC BRCA cohort
#####

setwd('~/Data/SMC_BRCA/')

# Load libraries
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(readxl)

# Load SMC data and tally mutation contributions
smc_mutations <- read.delim('data_mutations.txt')

mut.maf <- read.maf(smc_mutations, 
                    vc_nonSyn = names(table(smc_mutations$Variant_Classification)))

mt_tally.smc <- sig_tally(
  mut.maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ALL',
  use_syn = TRUE
)

smc_muts <- as.data.frame(cbind(mt_tally.smc$SBS_96, mt_tally.smc$ID_83))

## Load relevant data for classifier

setwd('~/Projects/HRD_MutationalSignature/')

# Prior cluster mean distributions (cluster_distributions = mut.dists_mean)
load('Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.Rdata')

# Signature Phenotype Assignment (cluster_assign = pheno_assigned)
load('Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.Rdata')
pheno_assigned <- ann$Phenotype

# Likelihood function that aligns a dataset with the designated mean distributions
likelihood_calc <- function(input_data, cluster_distributions, cluster_assign) {
  
  # This function:
  #   Takes a dataset of mutations as its input (rows = samples, cols = mutation types)
  #   For now, this function does not limit mutation types: SBS and indels will be included
  #   Applies a likelihood approach to calculate posterior probabilities of cluster assignment
  #   NB, For now, we ensure that our cluster distributions have the correct mutation types
  #   Output: posterior distributions
  
  # Calculate likelihood: P(mutation spectrum | cluster)
  log_likelihoods <- matrix(NA, nrow = nrow(input_data), ncol = nrow(cluster_distributions))
  rownames(log_likelihoods) <- rownames(input_data); colnames(log_likelihoods) <- rownames(cluster_distributions)
  print('Calculating log-likelihoods...')
  for (i in 1:nrow(input_data)) {
    print(paste0('Calculating log likelihoods for sample ', i, ' of ', nrow(input_data), ': ', rownames(input_data)[i]))
    log_likelihoods[i, ] <- apply(cluster_distributions, 1, 
                                  function(x) sum(log10(x) * input_data[i, ]))
  }
  
  # Set prior: P(cluster)
  marginal.probs <- table(cluster_assign)/length(cluster_assign)
  marginal.probs <- marginal.probs[colnames(log_likelihoods)]
  
  # Calculate log_posteriors and return final posteriors
  log_posteriors <- data.frame(t(apply(log_likelihoods, 1, function(x) 
    x + log10(marginal.probs)
  )))
  
  # Generate final posteriors
  final_probs <- log_posteriors
  for (i in 1:ncol(final_probs)) {
    final_probs[,i] <- 10^(apply(log_posteriors, 1, function(x) x[i] - max(x[1:ncol(final_probs)])))
  }
  
  final_probs <- data.frame(t(apply(final_probs, 1, function(x) x/sum(x))))
  
  return(final_probs)
  
}

# Apply log-likelihood approach
results.smc_loglik <- likelihood_calc(input_data = smc_muts, 
                                       cluster_distributions = mut.dists_mean,
                                       cluster_assign = pheno_assigned)

results.smc_df <- data.frame(
  Patient = rownames(results.smc_loglik),
  Phenotype_Assigned = apply(results.smc_loglik, 1,
                             function(x) names(x)[x==max(x)]),
  Phenotype_Assigned.prob = apply(results.smc_loglik, 1, max),
  HRD_prob = apply(results.smc_loglik[,grepl(pattern = 'HRD', names(results.smc_loglik))],
                   1, sum)
)
results.smc_df$HRD <- sapply(results.smc_df$HRD_prob,
                              function(x) ifelse(x >= 0.79, 'HRD', 'HR-proficient'))
save(results.smc_df, file = 'Results/SMC_HRD_resultsSummary.Rdata')

# Match with BRCA status
data.brca <- read_excel('~/Data/SMC_BRCA/SMC_BRCA.BrcaStatus.xlsx', skip=2)

results.smc_df$sample_id <- sapply(results.smc_df$Patient, 
                                   function(x) substr(x, 15, nchar(x)))
results.smc_df <- merge(x = results.smc_df, y = data.brca[,c('sample_id', 'gene_symbol')],
                        all.x = TRUE)

# Match with additional clinical data, specifically BRCA subtype
smc_clinic <- read.delim('~/Data/SMC_BRCA/data_clinical_sample.txt', skip=4)

results.smc_df <- merge(x = results.smc_df, y = smc_clinic[,c('PATIENT_ID', 'SUBTYPE_CONSENSUS')],
                        by.x = 'Patient', by.y = 'PATIENT_ID')
names(results.smc_df)[ncol(results.smc_df)] <- 'Subtype'
names(results.smc_df)[ncol(results.smc_df)-1] <- 'BRCA_defect'

# Plot results
ann_smc <- results.smc_df[,c('BRCA_defect','HRD_prob', 'HRD', 'Phenotype_Assigned', 'Subtype')]
rownames(ann_smc) <- results.smc_df$Patient
ann_smc <- ann_smc[order(ann_smc$Phenotype_Assigned), ]

results.smc_plot <- as.data.frame(t(results.smc_loglik[rownames(ann_smc), order(colnames(results.smc_loglik))]))

# Sort colours
cols <- colorRampPalette(brewer.pal(9,'Set1'))
cols_pheno <- cols(length(unique(ann_smc$Phenotype_Assigned)))
names(cols_pheno) <- unique(ann_smc$Phenotype_Assigned)

ann_smc_colours <- list(
  Phenotype_Assigned = cols_pheno,
  HRD = c('HRD' = 'black', 'HR-proficient' = 'white'),
  BRCA_defect = c('BRCA1' = 'blue', 'BRCA2' = 'red'),
  Subtype = c('ER+' = 'navy', 'HER2+' = 'darkgreen',
              'ER+HER2+' = 'gray', 'TN' = 'yellow')
)
ann_smc <- ann_smc[,ncol(ann_smc):1]

# Use white -> navy scale
cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)

pheatmap(results.smc_plot, 
         show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_col = ann_smc, 
         annotation_colors = ann_smc_colours,
         color = cols_scale, fontsize = 8, fontsize_row = 10,
         filename = 'Figures/Supp_SMCheatmap.pdf')

# Compare with BRCA subtype
table(ann_smc$HRD_prob >= 0.79, ann_smc$BRCA_defect, useNA = 'always')

# Plotting barplots
ann_smc$BRCA_status <- sapply(ann_smc$BRCA_defect, function(x)
  ifelse(!is.na(x), 'BRCA-defective', 'BRCA+'))
ann_smc$BRCA_status <- factor(ann_smc$BRCA_status,
                              levels = c('BRCA-defective','BRCA+'))

ann_smc$HRDgroup <- 'HR-proficient'
ann_smc$HRDgroup[ann_smc$HRD_prob >= 0.5] <- 'HRD > 0.5'
ann_smc$HRDgroup[ann_smc$HRD_prob >= 0.79] <- 'HRD > 0.79'
ann_smc$HRDgroup <- factor(ann_smc$HRDgroup,
                           levels = c('HR-proficient','HRD > 0.5', 'HRD > 0.79'))

ann_smc.plot1 <- ann_smc %>%
  group_by(BRCA_status, HRDgroup) %>%
  summarise(n = n())

g_brca <- ggplot(ann_smc.plot1, aes(x = BRCA_status, y = n, fill = HRDgroup)) +
  geom_bar(stat = 'identity', position = 'fill') + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylab('% Samples') +
  scale_fill_brewer(palette = 'Blues')
ggsave(filename = 'Figures/Supp_SMCbrcaClassification.pdf', plot = g_brca,
       height = 4, width = 4)

ann_smc.plot2 <- ann_smc %>%
  group_by(Subtype, HRDgroup) %>%
  summarise(n = n())

g_subtype <- ggplot(ann_smc.plot2, aes(x = Subtype, y = n, fill = HRDgroup)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylab('% Samples') +
  scale_fill_brewer(palette = 'Blues')
ggsave(filename = 'Figures/Supp_SMCsubtypeClassification.pdf', plot = g_subtype,
       height = 4, width = 5)


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/TCGA_geneCNAenrichment.R`:

```````R
#####
## Fisher's tests to determine differential chromosome Gene alterations in HRD vs HR-proficient samples
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(ggplot2)
library(wesanderson)
library(ggrepel)

# Load Gene-level CNA data, and organise to remove metadata
cna <- read.delim('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_cna.txt')
cna <- cna[!duplicated(cna$Hugo_Symbol), ]

# Extract cancer genes with alterations in >5% cases
cgc_data <- read.csv('~/Data/Census_allMon Jul  3 15_56_40 2023.csv')
cgc_data.genes <- as.character(cgc_data$Gene.Symbol)

cna <- cna[cna$Hugo_Symbol %in% cgc_data.genes, ]

# Sort dataframe
rownames(cna) <- cna$Hugo_Symbol
cna <- cna[,-c(1,2)]
cna <- as.data.frame(t(cna))
cna$Patient <- sapply(rownames(cna), function(x)
  paste(strsplit(x,split='[.]')[[1]][1:3], collapse='-'))

# Load HRD/HR-proficiency labels and merge with CNA data
load('Results/TCGA_HRDclassification_BRCAannotation.Rdata')
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga$group <- factor(ann_tcga$HRD, levels = c('HR-proficient','HRD'))

df.ann <- data.frame(
  Patient = ann_tcga$Patient,
  group = ann_tcga$group
)

df <- merge(x = df.ann, y = cna)

# Initialise data frames to track gain/loss enrichments
fishers.gain = fishers.loss <- data.frame(
  Gene = names(df)[-c(1,2)], Estimate = NA, pVal = NA
)

# For each chromosome Gene:
#   Conduct a Fisher's exact test for HRD and HR-proficient samples against:
#     1) Gain vs not-Gain
#     2) Loss vs not-Loss
#   Save odds ratios and p-values in initialised datasets
for (i in 1:nrow(fishers.gain)) {
  Gene.i <- fishers.gain$Gene[i]
  
  df.i <- data.frame(
    group = df$group,
    Gain = df[,Gene.i] > 0,
    Loss = df[,Gene.i] < 0
  )
  df.i <- df.i[!is.na(df.i[,'Gain']), ]
  
  # df.i$group <- factor(df.i$group, levels = c('BRCA-defective', 'BRCA+'))
  df.i$Gain <- factor(df.i$Gain, levels = c(FALSE,TRUE))
  df.i$Loss <- factor(df.i$Loss, levels = c(FALSE,TRUE))
  
  # Gains
  table.gain.i <- table(df.i$group, df.i$Gain)
  fishers.gain.i <- fisher.test(table.gain.i)
  
  fishers.gain$Estimate[i] <- fishers.gain.i$estimate
  fishers.gain$pVal[i] <- fishers.gain.i$p.value
  
  # Losses
  table.loss.i <- table(df.i$group, df.i$Loss)
  fishers.loss.i <- fisher.test(table.loss.i)
  
  fishers.loss$Estimate[i] <- fishers.loss.i$estimate
  fishers.loss$pVal[i] <- fishers.loss.i$p.value
  
}

## Add information to gain/loss results:
#   - Adjust p-values
#   - If Estimate == 'Inf', set to maximum estimate
#   - If Estimate == 0, set to minimum estimate
#   - log2-normalise estimates and -log10(p-adjust)
#   - Add labels for Genes with enrichment significance < .01
#   - 'Volcano plot' displaying results

# Gains
fishers.gain$padj <- p.adjust(fishers.gain$pVal)
# fishers.gain$padj <- fishers.gain$pVal
fishers.gain$Estimate[fishers.gain$Estimate == 'Inf'] <- max(fishers.gain$Estimate)
fishers.gain$Estimate[fishers.gain$Estimate == 0] <- min(fishers.gain$Estimate)
fishers.gain$l2fc <- log2(fishers.gain$Estimate)
fishers.gain$logp <- -log10(fishers.gain$padj)
fishers.gain$label <- fishers.gain$Gene
fishers.gain$label[fishers.gain$padj > .05] <- NA

g_gains <- ggplot(fishers.gain, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Gain vs !Gain')
# ggsave(filename = 'Figures/Supplementary/ChrGeneEnrich_GainVsNotGain.pdf',
#        plot = g_gains)

# Losses
fishers.loss$padj <- p.adjust(fishers.loss$pVal)
# fishers.loss$padj <- fishers.loss$pVal
fishers.loss$Estimate[fishers.loss$Estimate == 'Inf'] <- max(fishers.loss$Estimate)
fishers.loss$Estimate[fishers.loss$Estimate == 0] <- min(fishers.loss$Estimate)
fishers.loss$l2fc <- log2(fishers.loss$Estimate)
fishers.loss$logp <- -log10(fishers.loss$padj)
fishers.loss$label <- fishers.loss$Gene
fishers.loss$label[fishers.loss$padj > .05] <- NA

g_losses <- ggplot(fishers.loss, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Loss vs !Loss')
# ggsave(filename = 'Figures/Supplementary/ChrGeneEnrich_LossVsNotLoss.pdf',
#        plot = g_losses)

# Merge log-fold changes (+ labels, showing significant Genes for each test)
df.fishers <- merge(x = fishers.gain[,c('Gene','l2fc','label')],
                    y = fishers.loss[,c('Gene','l2fc','label')], by = 'Gene')
names(df.fishers) <- c('Gene','l2fc_GAIN', 'label_GAIN','l2fc_LOSS','label_LOSS')

# Add GAIN/LOSS/GAIN+LOSS labels for which tests are significant
df.fishers$Label <- 'none'
df.fishers$Label[!is.na(df.fishers$label_GAIN)] <- 'GAIN'
df.fishers$Label[!is.na(df.fishers$label_LOSS)] <- 'LOSS'
df.fishers$Label[!is.na(df.fishers$label_GAIN) & !is.na(df.fishers$label_LOSS)] <- 'GAIN+LOSS'
df.fishers$Gene_Label <- NA
df.fishers$Gene_Label[df.fishers$Label != 'none'] <- df.fishers$Gene[df.fishers$Label != 'none']

# df.fishers$Label <- factor(df.fishers$Label, levels = c('GAIN','GAIN+LOSS','LOSS','none'))

g_GainvsLoss <- ggplot(df.fishers, aes(x = l2fc_GAIN, y = l2fc_LOSS, color = Label, label = Gene_Label)) +
  geom_point() + geom_text_repel() + theme_minimal() +
  scale_color_manual(values = c(wes_palette('Zissou1')[1], wes_palette('Zissou1')[5],
                                wes_palette('Zissou1')[4], 'gray90')) +
  xlab('log2(Fold Change) - Gain') + ylab('log2(Fold Change) - Loss') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'top',
        legend.title = element_blank())
ggsave(filename = 'Figures/Supp_GeneCnaEnrich_GainVsLoss.pdf',
       plot = g_GainvsLoss, height = 4, width = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/decoupleR_full.R`:

```````R
setwd('~/Data/TCGA')

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tibble)
library(decoupleR)
library(ggplot2)

query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)
# GDCdownload(query)
data <- GDCprepare(query = query)
data <- data[, data$sample_type == 'Primary Tumor']
data <- data[, !duplicated(data$patient)]

# Match with HRD classification
load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
# ann_tcga$Patient <- rownames(ann_tcga)
ann_tcga$HRD <- sapply(ann_tcga$HRD_prob, function(x) ifelse(x>=0.79,'HRD','HR-proficient'))
ann_tcga <- ann_tcga[!is.na(ann_tcga$ER_status), ]
# ann_tcga <- ann_tcga[ann_tcga$ER_status == 'Negative', ]
ann_tcga <- ann_tcga[,c('Patient','HRD')]

patients.intersect <- intersect(ann_tcga$Patient, data$patient)
ann_tcga <- ann_tcga[match(patients.intersect, ann_tcga$Patient), ]
data <- data[, match(patients.intersect, data$patient)]

data$HRD_status <- factor(ann_tcga$HRD, levels = c('HR-proficient','HRD'))

# Construct a DESeqDataSet data object
dds <- DESeqDataSet(data, design = ~ HRD_status)

# Gene count filtering
dds <- dds[rowSums(counts(dds)) > 10, ]

# Normalisation
dds <- estimateSizeFactors(dds)

# DIFFERENTIAL GENE EXPRESSION ANALYSIS
dds_DGE <- DESeq(dds)
dds_DGE_results <- results(dds_DGE)

dds_statVals <- data.frame(
  EnsemblID = rownames(dds_DGE_results),
  stat = dds_DGE_results$stat
)
dds_statVals$ID <- rowData(dds)$gene_name[match(rownames(rowData(dds)), dds_statVals$EnsemblID)]
dds_statVals <- dds_statVals[!duplicated(dds_statVals$ID), ]
rownames(dds_statVals) <- NULL

deg <- dds_statVals %>%
  select(ID, stat) %>%
  column_to_rownames(var = 'ID') %>%
  as.matrix()

counts <- assay(dds)
counts <- counts[dds_statVals$EnsemblID, ]
rownames(counts) <- dds_statVals$ID

counts_logNorm <- log2(counts + 1)
colnames(counts_logNorm) <- sapply(colnames(counts_logNorm), function(x) substr(x,1,12))

design <- ann_tcga
names(design) <- c('sample','condition')

# Run fucking progeny
net <- get_progeny(organism = 'human', top = 100)

# Run mlm
contrast_acts <- run_mlm(mat = deg, net = net, .source = 'source',
                         .target = 'target', .mor = 'weight', minsize = 5)
contrast_acts

# Plot
g_progeny <- ggplot(contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab('Pathways') + ylab('Enrichment in HRD Samples')
ggsave(filename = '~/Projects/HRD_MutationalSignature/Figures/Figure2/HRD_decoupleR.pdf',
       plot = g_progeny, height = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/TCGA_chromosomeArmAlterations_HRDonly.R`:

```````R
#####
## Fisher's tests to determine differential chromosome arm alterations in HRD vs HR-proficient samples
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(ggplot2)
library(wesanderson)
library(ggrepel)

# Load arm-level CNA data, and organise to remove metadata
cna <- read.delim('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_armlevel_cna.txt')
rownames(cna) <- cna$NAME
cna <- cna[,-c(1,2,3)]
cna <- as.data.frame(t(cna))
cna$Patient <- sapply(rownames(cna), function(x)
  paste(strsplit(x,split='[.]')[[1]][1:3], collapse='-'))

# Load HRD/HR-proficiency labels and merge with CNA data
load('Results/TCGA_HRDclassification_BRCAannotation.Rdata')
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[ann_tcga$HRD == 'HRD', ]
ann_tcga$group <- sapply(ann_tcga$BRCA_status,
                         function(x) ifelse(x == 'none', 'BRCA+', 'BRCA-defective'))

df.ann <- data.frame(
  Patient = ann_tcga$Patient,
  group = ann_tcga$group
)

df <- merge(x = df.ann, y = cna)

# Initialise data frames to track gain/loss enrichments
fishers.gain = fishers.loss <- data.frame(
  Arm = names(df)[-c(1,2)], Estimate = NA, pVal = NA
)

# For each chromosome arm:
#   Conduct a Fisher's exact test for HRD and HR-proficient samples against:
#     1) Gain vs not-Gain
#     2) Loss vs not-Loss
#   Save odds ratios and p-values in initialised datasets
for (i in 1:nrow(fishers.gain)) {
  arm.i <- fishers.gain$Arm[i]
  
  df.i <- data.frame(
    group = df$group,
    Gain = df[,arm.i] == 'Gain',
    Loss = df[,arm.i] == 'Loss'
  )
  df.i <- df.i[!is.na(df.i[,'Gain']), ]
  
  df.i$group <- factor(df.i$group, levels = c('BRCA-defective', 'BRCA+'))
  df.i$Gain <- factor(df.i$Gain, levels = c(FALSE,TRUE))
  df.i$Loss <- factor(df.i$Loss, levels = c(FALSE,TRUE))
  
  # Gains
  table.gain.i <- table(df.i$group, df.i$Gain)
  fishers.gain.i <- fisher.test(table.gain.i)
  
  fishers.gain$Estimate[i] <- fishers.gain.i$estimate
  fishers.gain$pVal[i] <- fishers.gain.i$p.value
  
  # Losses
  table.loss.i <- table(df.i$group, df.i$Loss)
  fishers.loss.i <- fisher.test(table.loss.i)
  
  fishers.loss$Estimate[i] <- fishers.loss.i$estimate
  fishers.loss$pVal[i] <- fishers.loss.i$p.value
  
}

## Add information to gain/loss results:
#   - Adjust p-values
#   - If Estimate == 'Inf', set to maximum estimate
#   - If Estimate == 0, set to minimum estimate
#   - log2-normalise estimates and -log10(p-adjust)
#   - Add labels for arms with enrichment significance < .01
#   - 'Volcano plot' displaying results

# Gains
fishers.gain$padj <- p.adjust(fishers.gain$pVal)
fishers.gain$Estimate[fishers.gain$Estimate == 'Inf'] <- max(fishers.gain$Estimate)
fishers.gain$Estimate[fishers.gain$Estimate == 0] <- min(fishers.gain$Estimate)
fishers.gain$l2fc <- log2(fishers.gain$Estimate)
fishers.gain$logp <- -log10(fishers.gain$padj)
fishers.gain$label <- fishers.gain$Arm
fishers.gain$label[fishers.gain$padj > .05] <- NA

g_gains <- ggplot(fishers.gain, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Gain vs !Gain')
# ggsave(filename = 'Figures/Supplementary/ChrArmEnrich_GainVsNotGain.pdf',
#        plot = g_gains)

# Losses
fishers.loss$padj <- p.adjust(fishers.loss$pVal)
fishers.loss$Estimate[fishers.loss$Estimate == 'Inf'] <- max(fishers.loss$Estimate)
fishers.loss$Estimate[fishers.loss$Estimate == 0] <- min(fishers.loss$Estimate)
fishers.loss$l2fc <- log2(fishers.loss$Estimate)
fishers.loss$logp <- -log10(fishers.loss$padj)
fishers.loss$label <- fishers.loss$Arm
fishers.loss$label[fishers.loss$padj > .05] <- NA

g_losses <- ggplot(fishers.loss, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Loss vs !Loss')
# ggsave(filename = 'Figures/Supplementary/ChrArmEnrich_LossVsNotLoss.pdf',
#        plot = g_losses)

# Merge log-fold changes (+ labels, showing significant arms for each test)
df.fishers <- merge(x = fishers.gain[,c('Arm','l2fc','label')],
                    y = fishers.loss[,c('Arm','l2fc','label')], by = 'Arm')
names(df.fishers) <- c('Arm','l2fc_GAIN', 'label_GAIN','l2fc_LOSS','label_LOSS')

# Add GAIN/LOSS/GAIN+LOSS labels for which tests are significant
df.fishers$Label <- 'none'
df.fishers$Label[!is.na(df.fishers$label_GAIN)] <- 'GAIN'
df.fishers$Label[!is.na(df.fishers$label_LOSS)] <- 'LOSS'
df.fishers$Label[!is.na(df.fishers$label_GAIN) & !is.na(df.fishers$label_LOSS)] <- 'GAIN+LOSS'
df.fishers$Arm_Label <- NA
df.fishers$Arm_Label[df.fishers$Label != 'none'] <- df.fishers$Arm[df.fishers$Label != 'none']

# df.fishers$Label <- factor(df.fishers$Label, levels = c('GAIN','GAIN+LOSS','LOSS','none'))

g_GainvsLoss <- ggplot(df.fishers, aes(x = l2fc_GAIN, y = l2fc_LOSS, color = Label, label = Arm_Label)) +
  geom_point() + geom_text_repel() + theme_minimal() +
  scale_color_manual(values = c(wes_palette('Zissou1')[1], 'gray90')) +
  xlab('log2(Fold Change) - Gain') + ylab('log2(Fold Change) - Loss') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'top',
        legend.title = element_blank())
ggsave(filename = 'Figures/Supp_ChrArmEnrich_GainVsLoss_HRDonly.pdf',
       plot = g_GainvsLoss, height = 4, width = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/HRD_Hypoxia.R`:

```````R
setwd('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/')

library(introdataviz)
library(wesanderson)

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
hypoxia <- read.table('data_clinical_supp_hypoxia.txt', h=T)

ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status) &
                       !(ann_tcga$BRCA_status %in% c('PALB2','RAD51C')), ]
ann_tcga <- ann_tcga[ann_tcga$ER_status %in% c('Negative','Positive'), ]

df <- merge(x = ann_tcga[,c('Patient','HRD','BRCA_status','ER_status')],
            y = hypoxia[,c('PATIENT_ID','BUFFA_HYPOXIA_SCORE')],
            by.x = 'Patient', by.y = 'PATIENT_ID')
names(df)[5] <- 'Hypoxia_Buffa'
df$ER_status <- sapply(df$ER_status, function(x) 
  ifelse(x == 'Negative', 'ER-', 'ER+'))
df$ER_status <- factor(df$ER_status, levels = c('ER+', 'ER-'))

library(ggpubr)

df$Group <- df$BRCA_status
df$Group[df$BRCA_status == 'none' & df$HRD == 'HRD'] <- 'HRD_BRCA+'
df$Group[df$Group == 'none'] <- 'HR-proficient'
df$Group <- factor(df$Group, levels = c('HR-proficient','HRD_BRCA+',
                                        'BRCA1','BRCA2'))

comps <- list(c('HRD_BRCA+','BRCA2'),c('HR-proficient','HRD_BRCA+'),c('HRD_BRCA+','BRCA1'))

g_hypoxia <- ggplot(df, aes(x = Group, y = Hypoxia_Buffa, fill = ER_status)) +
  geom_split_violin(alpha = 0.4) + geom_boxplot() +
  stat_compare_means(comparisons = comps) +
  scale_fill_manual(values = wes_palette('Royal1')[c(3,2)]) +
  ylab('Hypoxia Score (Buffa)') + theme_minimal() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())
# ggsave(filename = '~/Projects/HRD_MutationalSignature/Figures/Figure2/HRD_Hypoxia.pdf',
#        plot = g_hypoxia, height = 4, width = 7)
ggsave(filename = '~/Projects/Thesis/Chapter 4/HRD_Hypoxia.pdf',
       plot = g_hypoxia)

ggplot(df[df$ER_status == 'ER+', ], aes(x = Group, y = Hypoxia_Buffa)) +
  geom_boxplot() +
  stat_compare_means(comparisons = comps)

lm_hypoxia <- lm(Hypoxia_Buffa ~ ER_status + BRCA_status*HRD, data = df)
anova(lm_hypoxia)


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/TCGA_chromosomeArmAlterations.R`:

```````R
#####
## Fisher's tests to determine differential chromosome arm alterations in HRD vs HR-proficient samples
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(ggplot2)
library(wesanderson)
library(ggrepel)

# Load arm-level CNA data, and organise to remove metadata
cna <- read.delim('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_armlevel_cna.txt')
rownames(cna) <- cna$NAME
cna <- cna[,-c(1,2,3)]
cna <- as.data.frame(t(cna))
cna$Patient <- sapply(rownames(cna), function(x)
  paste(strsplit(x,split='[.]')[[1]][1:3], collapse='-'))

# Load HRD/HR-proficiency labels and merge with CNA data
load('Results/TCGA_HRD_resultsSummary.Rdata')
df.ann <- data.frame(
  Patient = results.tcga_df$Patient,
  HRD = results.tcga_df$HRD
)

df <- merge(x = df.ann, y = cna)

# Initialise data frames to track gain/loss enrichments
fishers.gain = fishers.loss <- data.frame(
  Arm = names(df)[-c(1,2)], Estimate = NA, pVal = NA
)

# For each chromosome arm:
#   Conduct a Fisher's exact test for HRD and HR-proficient samples against:
#     1) Gain vs not-Gain
#     2) Loss vs not-Loss
#   Save odds ratios and p-values in initialised datasets
for (i in 1:nrow(fishers.gain)) {
  arm.i <- fishers.gain$Arm[i]
  
  df.i <- data.frame(
    HRD = df$HRD,
    Gain = df[,arm.i] == 'Gain',
    Loss = df[,arm.i] == 'Loss'
  )
  df.i <- df.i[!is.na(df.i[,'Gain']), ]
  
  df.i$HRD <- factor(df.i$HRD, levels = c('HR-proficient', 'HRD'))
  
  # Gains
  table.gain.i <- table(df.i$HRD, df.i$Gain)
  fishers.gain.i <- fisher.test(table.gain.i)
  
  fishers.gain$Estimate[i] <- fishers.gain.i$estimate
  fishers.gain$pVal[i] <- fishers.gain.i$p.value
  
  # Losses
  table.loss.i <- table(df.i$HRD, df.i$Loss)
  fishers.loss.i <- fisher.test(table.loss.i)
  
  fishers.loss$Estimate[i] <- fishers.loss.i$estimate
  fishers.loss$pVal[i] <- fishers.loss.i$p.value
  
}

## Add information to gain/loss results:
#   - Adjust p-values
#   - If Estimate == 'Inf', set to maximum estimate
#   - If Estimate == 0, set to minimum estimate
#   - log2-normalise estimates and -log10(p-adjust)
#   - Add labels for arms with enrichment significance < .01
#   - 'Volcano plot' displaying results

# Gains
fishers.gain$padj <- p.adjust(fishers.gain$pVal)
fishers.gain$Estimate[fishers.gain$Estimate == 'Inf'] <- max(fishers.gain$Estimate)
fishers.gain$Estimate[fishers.gain$Estimate == 0] <- min(fishers.gain$Estimate)
fishers.gain$l2fc <- log2(fishers.gain$Estimate)
fishers.gain$logp <- -log10(fishers.gain$padj)
fishers.gain$label <- fishers.gain$Arm
fishers.gain$label[fishers.gain$padj > .01] <- NA

g_gains <- ggplot(fishers.gain, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Gain vs !Gain')
# ggsave(filename = 'Figures/Supplementary/ChrArmEnrich_GainVsNotGain.pdf',
#        plot = g_gains)

# Losses
fishers.loss$padj <- p.adjust(fishers.loss$pVal)
fishers.loss$Estimate[fishers.loss$Estimate == 'Inf'] <- max(fishers.loss$Estimate)
fishers.loss$Estimate[fishers.loss$Estimate == 0] <- min(fishers.loss$Estimate)
fishers.loss$l2fc <- log2(fishers.loss$Estimate)
fishers.loss$logp <- -log10(fishers.loss$padj)
fishers.loss$label <- fishers.loss$Arm
fishers.loss$label[fishers.loss$padj > .01] <- NA

g_losses <- ggplot(fishers.loss, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Loss vs !Loss')
# ggsave(filename = 'Figures/Supplementary/ChrArmEnrich_LossVsNotLoss.pdf',
#        plot = g_losses)

# Merge log-fold changes (+ labels, showing significant arms for each test)
df.fishers <- merge(x = fishers.gain[,c('Arm','l2fc','label')],
                    y = fishers.loss[,c('Arm','l2fc','label')], by = 'Arm')
names(df.fishers) <- c('Arm','l2fc_GAIN', 'label_GAIN','l2fc_LOSS','label_LOSS')

# Add GAIN/LOSS/GAIN+LOSS labels for which tests are significant
df.fishers$Label <- 'none'
df.fishers$Label[!is.na(df.fishers$label_GAIN)] <- 'GAIN'
df.fishers$Label[!is.na(df.fishers$label_LOSS)] <- 'LOSS'
df.fishers$Label[!is.na(df.fishers$label_GAIN) & !is.na(df.fishers$label_LOSS)] <- 'GAIN+LOSS'
df.fishers$Arm_Label <- NA
df.fishers$Arm_Label[df.fishers$Label != 'none'] <- df.fishers$Arm[df.fishers$Label != 'none']

g_GainvsLoss <- ggplot(df.fishers, aes(x = l2fc_GAIN, y = l2fc_LOSS, color = Label, label = Arm_Label)) +
  geom_point() + geom_text_repel() + theme_minimal() +
  scale_color_manual(values = c(wes_palette('Zissou1')[1], wes_palette('Zissou1')[5],
                                wes_palette('Zissou1')[4], 'gray90')) +
  xlab('log2(Fold Change) - Gain') + ylab('log2(Fold Change) - Loss') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'top',
        legend.title = element_blank())
ggsave(filename = 'Figures/Figure2/ChrArmEnrich_GainVsLoss.pdf',
       plot = g_GainvsLoss, height = 4, width = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/dNdS_Analysis.R`:

```````R
#####
## dN/dS analysis to identify mutations under positive selection in HRD/BRCA groups
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(TCGAbiolinks)
library(maftools)
library(dndscv)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(wesanderson)

# Load data and dN/dS references
load('~/Data/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda')

setwd('~/Data/TCGA')
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Simple Nucleotide Variation',
  data.type = 'Masked Somatic Mutation'
)
# GDCdownload(query)
mutations <- GDCprepare(query = query)

data <- read.maf(maf = mutations, isTCGA = TRUE)
data <- rbind(data@data, data@maf.silent)
dat_small <- data[,c('Tumor_Sample_Barcode', 'Chromosome',
                     'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')]
names(dat_small) <- c('sampleID', 'chr', 'pos', 'ref', 'mut')

# Load HRD/BRCA groups and arrange data accordingly (note, samples can overlap)
setwd('~/Projects/HRD_MutationalSignature/')

load('Results/TCGA_HRDclassification_BRCAannotation.Rdata')
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ] # remove samples with missing BRCA_status
# ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('RAD51C','PALB2')), ] # remove samples with deficiencies in other HR genes
dat_small <- merge(x = dat_small, y = ann_tcga[,c('Patient', 'HRD', 'BRCA_status')],
                   by.x = 'sampleID', by.y = 'Patient')

# HRD vs HR-proficient
data_HRD <- dat_small[dat_small$HRD == 'HRD', 1:5]
data_HRproficient <- dat_small[dat_small$HRD == 'HR-proficient', 1:5]

# BRCA-defect categories
data_BRCA1 <- dat_small[dat_small$BRCA_status == 'BRCA1', 1:5]
data_BRCA2 <- dat_small[dat_small$BRCA_status == 'BRCA2', 1:5]
data_RAD51C <- dat_small[dat_small$BRCA_status == 'RAD51C', 1:5]
data_HRDBRCApos <- dat_small[dat_small$HRD == 'HRD' &
                               dat_small$BRCA_status == 'none', 1:5]

## Run dN/dS analysis on each group separately

# BRCA1
dndsoutBRCA1 <- dndscv(data_BRCA1, cv = NULL,
                       refdb = RefCDS)
sel_cvBRCA1 <- dndsoutBRCA1$sel_cv

# BRCA2
dndsoutBRCA2 <- dndscv(data_BRCA2, cv = NULL,
                       refdb = RefCDS)
sel_cvBRCA2 <- dndsoutBRCA2$sel_cv

# RAD51C
dndsoutRAD51C <- dndscv(data_RAD51C, cv = NULL,
                        refdb = RefCDS)
sel_cvRAD51C <- dndsoutRAD51C$sel_cv

# HRD_BRCA+
dndsoutHRDBRCApos <- dndscv(data_HRDBRCApos, cv = NULL,
                            refdb = RefCDS)
sel_cvHRDBRCApos <- dndsoutHRDBRCApos$sel_cv

# HRD full
dndsoutHRD <- dndscv(data_HRD, cv = NULL,
                     refdb = RefCDS)
sel_cvHRD <- dndsoutHRD$sel_cv

# HR-proficient
dndsoutHRprof <- dndscv(data_HRproficient, cv = NULL,
                        refdb = RefCDS)
sel_cvHRprof <- dndsoutHRprof$sel_cv


## Plot HRD vs HR-proficient gene selections

# Find genes under positive selection in at least one group
sig_genes <- unique(c(sel_cvHRD$gene_name[sel_cvHRD$qind_cv < .1],
                      sel_cvHRprof$gene_name[sel_cvHRprof$qind_cv < .1]))

sel_cvSigGenes <- merge(
  x = sel_cvHRD[sel_cvHRD$gene_name %in% sig_genes, c('gene_name', 'wind_cv', 'qind_cv')],
  y = sel_cvHRprof[sel_cvHRprof$gene_name %in% sig_genes, c('gene_name', 'wind_cv', 'qind_cv')],
  by = 'gene_name'
)
names(sel_cvSigGenes)[-1] <- c('dNdS_HRD', 'FDR_HRD',
                               'dNdS_HRprof', 'FDR_HRprof')

# Prepare labels for plotting
sel_cvSigGenes$Significant <- 'Both'
sel_cvSigGenes$Significant[sel_cvSigGenes$FDR_HRD > .1] <- 'HR-proficient ONLY'
sel_cvSigGenes$Significant[sel_cvSigGenes$FDR_HRprof > .1] <- 'HRD ONLY'

sel_cvSigGenes$dNdS_HRD <- log2(sel_cvSigGenes$dNdS_HRD + 1)
sel_cvSigGenes$dNdS_HRprof <- log2(sel_cvSigGenes$dNdS_HRprof + 1)

# save(sel_cvSigGenes, file = 'Results/ExomeClassifier/TCGA_BRCA/dNdS_HRDvsHRprof.Rdata')

ggplot(sel_cvSigGenes, aes(x = dNdS_HRD, y = dNdS_HRprof, 
                           color = Significant, label = gene_name)) +
  geom_point() + theme_minimal() + geom_text_repel()

  
# Plot dN/dS of indsense variants across all groups
sel_cvFull <- rbind(sel_cvHRprof, sel_cvHRD, sel_cvBRCA1,
                    sel_cvBRCA2, sel_cvRAD51C, sel_cvHRDBRCApos)
sel_cvFull <- sel_cvFull[,c('gene_name','wmis_cv','qmis_cv')]
sel_cvFull$group <- rep(c('HR-proficient','HRD full','BRCA1',
                          'BRCA2','RAD51C','HRD BRCA+'),
                        each = nrow(sel_cvFull)/6)

sel_cvFull <- sel_cvFull[sel_cvFull$qmis_cv < 0.05, ]

sel_cvFull$log_dNdS <- log2(sel_cvFull$wmis_cv)
sel_cvFull$log_FDR <- -log10(sel_cvFull$qmis_cv)
sel_cvFull$log_FDR[sel_cvFull$log_FDR == 'Inf' | sel_cvFull$log_FDR >= 9] <- 9

sel_cvFull$group <- factor(sel_cvFull$group,
                           levels = c('HR-proficient','HRD full','BRCA1',
                                      'BRCA2','RAD51C','HRD BRCA+'))

g_dnds <- ggballoonplot(sel_cvFull, x = 'gene_name', y = 'group',
              size = 'log_dNdS', fill = 'log_FDR') +
  scale_fill_continuous(low = 'lightblue', high = 'navy') +
  geom_hline(yintercept = 2.5, col = 'red', linetype = 'dashed') +
  theme(legend.position = 'top')
# ggsave(filename = 'Figures/Figure2/HRD_dNdSballoonPlot.pdf',
#        plot = g_dnds, width = 4, height = 7)
ggsave(filename = '~/Projects/Thesis/Chapter 4/HRD_dNdSballoonPlot.pdf',
       plot = g_dnds, width = 7, height = 3)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/TCGA_HRDclassification.R`:

```````R
#####
## Apply HRD exome classifier to TCGA-BRCA cohort
#####

setwd('~/Data/TCGA/')

# Load libraries
library(plyr)
library(TCGAbiolinks)
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# Load TCGA data and tally mutation contributions
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Simple Nucleotide Variation',
  data.type = 'Masked Somatic Mutation'
)
# GDCdownload(query)
tcga_mutations <- GDCprepare(query = query)

# Exclude mutations from non-primary tumours
tcga_mutations$sample_type_code <- sapply(tcga_mutations$Tumor_Sample_Barcode,
                                          function(x) substr(x,14,15))
tcga_mutations <- tcga_mutations[tcga_mutations$sample_type_code == '01', ]

mut.maf <- read.maf(tcga_mutations, isTCGA = TRUE,
                    vc_nonSyn = names(table(tcga_mutations$Variant_Classification)))

mt_tally.tcga <- sig_tally(
  mut.maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg38',
  mode = 'ALL',
  use_syn = TRUE
)

# Check proportion of samples with SBS/ID loads >= 50 and plot distributions
sampleLoads <- data.frame(
  Sample = rownames(mt_tally.tcga$SBS_96),
  SBS = log2(apply(mt_tally.tcga$SBS_96, 1, sum)+1),
  ID = log2(apply(mt_tally.tcga$ID_83, 1, sum) + 1)
)

sampleLoads <- sampleLoads %>%
  pivot_longer(cols = -Sample, names_to = 'MutationType', values_to = 'log_count')
sampleLoads$MutationType <- factor(sampleLoads$MutationType,
                                   levels = c('SBS','ID'))

g_mutLoads <- gghistogram(data = sampleLoads, x = 'log_count',
                           fill = 'lightgrey') +
  geom_vline(xintercept = log2(50+1), col = 'red', linetype = 'dashed') +
  geom_vline(data = ddply(sampleLoads, "MutationType", summarize, wavg = median(log_count)), aes(xintercept=wavg),
             col = 'blue', linetype = 'dashed') +
  facet_wrap(~MutationType, scales = 'free')

ggsave(filename = '~/Projects/HRD_MutationalSignature/Figures/SupplementaryFigures/Supp_TCGAmutationLoads.pdf',
       plot = g_mutLoads, width = 6, height = 3)

table(sampleLoads$log_count[sampleLoads$MutationType == 'SBS'] >= log2(50+1))
table(sampleLoads$log_count[sampleLoads$MutationType == 'ID'] >= log2(50+1))

# Collate SBS_96 and ID_83 contributions and save results
tcga_muts <- as.data.frame(cbind(mt_tally.tcga$SBS_96, mt_tally.tcga$ID_83))

save(tcga_muts, file = 'TCGA_BRCA_mutContributions.Rdata')

load('TCGA_BRCA_mutContributions.Rdata') # if running without data loading

## Load relevant data for classifier

setwd('~/Projects/HRD_MutationalSignature/')

# Prior cluster mean distributions (cluster_distributions = mut.dists_mean)
load('Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.Rdata')

# Signature Phenotype Assignment (cluster_assign = pheno_assigned)
load('Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.Rdata')
pheno_assigned <- ann$Phenotype

# Likelihood function that aligns a dataset with the designated mean distributions
likelihood_calc <- function(input_data, cluster_distributions, cluster_assign) {
  
  # This function:
  #   Takes a dataset of mutations as its input (rows = samples, cols = mutation types)
  #   For now, this function does not limit mutation types: SBS and indels will be included
  #   Applies a likelihood approach to calculate posterior probabilities of cluster assignment
  #   NB, For now, we ensure that our cluster distributions have the correct mutation types
  #   Output: posterior distributions
  
  # Calculate likelihood: P(mutation spectrum | cluster)
  log_likelihoods <- matrix(NA, nrow = nrow(input_data), ncol = nrow(cluster_distributions))
  rownames(log_likelihoods) <- rownames(input_data); colnames(log_likelihoods) <- rownames(cluster_distributions)
  print('Calculating log-likelihoods...')
  for (i in 1:nrow(input_data)) {
    print(paste0('Calculating log likelihoods for sample ', i, ' of ', nrow(input_data), ': ', rownames(input_data)[i]))
    log_likelihoods[i, ] <- apply(cluster_distributions, 1, 
                                  function(x) sum(log10(x) * input_data[i, ]))
  }
  
  # Set prior: P(cluster)
  marginal.probs <- table(cluster_assign)/length(cluster_assign)
  marginal.probs <- marginal.probs[colnames(log_likelihoods)]
  
  # Calculate log_posteriors and return final posteriors
  log_posteriors <- data.frame(t(apply(log_likelihoods, 1, function(x) 
    x + log10(marginal.probs)
  )))
  
  # Generate final posteriors
  final_probs <- log_posteriors
  for (i in 1:ncol(final_probs)) {
    final_probs[,i] <- 10^(apply(log_posteriors, 1, function(x) x[i] - max(x[1:ncol(final_probs)])))
  }
  
  final_probs <- data.frame(t(apply(final_probs, 1, function(x) x/sum(x))))
  
  return(final_probs)
  
}

# Apply log-likelihood approach
results.tcga_loglik <- likelihood_calc(input_data = tcga_muts, 
                                       cluster_distributions = mut.dists_mean,
                                       cluster_assign = pheno_assigned)

results.tcga_df <- data.frame(
  Patient = rownames(results.tcga_loglik),
  Phenotype_Assigned = apply(results.tcga_loglik, 1,
                             function(x) names(x)[x==max(x)]),
  Phenotype_Assigned.prob = apply(results.tcga_loglik, 1, max),
  HRD_prob = apply(results.tcga_loglik[,grepl(pattern = 'HRD', names(results.tcga_loglik))],
                   1, sum)
)
results.tcga_df$HRD <- sapply(results.tcga_df$HRD_prob,
                              function(x) ifelse(x > .79, 'HRD', 'HR-proficient'))
# results.tcga_df$HRD_lenient <- sapply(results.tcga_df$HRD_prob,
#                                       function(x) ifelse(x >= .27, 'HRD', 'HR-proficient'))
# results.tcga_df$HRD_strict <- sapply(results.tcga_df$HRD_prob,
#                                      function(x) ifelse(x >= .79, 'HRD','HR-proficient'))
save(results.tcga_df, file = 'Results/TCGA_HRD_resultsSummary.Rdata')

# Compare with BRCA defects and clinical features

valieris <- read_excel('~/Data/TCGA/cancers-12-03687-s001/BRCA.xlsx', sheet = 'class-original')
valieris <- valieris[valieris$BRCA1_somatic_null != 'NA', ]
valieris <- valieris[,c('sample','event.BRCA1','event.BRCA2','event.RAD51C','event.PALB2')]
valieris$BRCA1 <- valieris$event.BRCA1 != 0
valieris$BRCA2 <- !(valieris$event.BRCA2 %in% c(0,'Mono-allelic-inactivation'))
valieris$RAD51C <- valieris$event.RAD51C != 0
valieris$PALB2 <- valieris$event.PALB2 != 0

valieris$BRCA_status <- 'none'
valieris$BRCA_status[valieris$PALB2] <- 'PALB2'
valieris$BRCA_status[valieris$RAD51C] <- 'RAD51C'
valieris$BRCA_status[valieris$BRCA2] <- 'BRCA2'
valieris$BRCA_status[valieris$BRCA1] <- 'BRCA1'

# Clinical subtypes
clin <- read.delim('~/Data/TCGA/TCGA_clinicalStatus.txt', header = TRUE)
tcga.clinical <- merge(x = valieris[,c('sample','BRCA_status')], y = clin[,c('bcr_patient_barcode','er_status_by_ihc')],
                       by.x = 'sample', by.y = 'bcr_patient_barcode')
names(tcga.clinical)[3] <- 'ER_status'

# Create annotation
ann_tcga <- merge(x = results.tcga_df[,c('Patient','Phenotype_Assigned', 'HRD', 'HRD_prob')],
                  y = tcga.clinical,
                  by.x = 'Patient', by.y = 'sample', all.x = TRUE)
ann_tcga <- ann_tcga[!duplicated(ann_tcga$Patient), ]
save(ann_tcga, file = 'Results/TCGA_HRDclassification_BRCAannotation.Rdata')

rownames(ann_tcga) <- ann_tcga$Patient; ann_tcga <- ann_tcga[,-1]
ann_tcga <- ann_tcga[order(ann_tcga$Phenotype_Assigned), ]

# Plots showing BRCA-defect classification as HRD
ann_tcga.plot <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga.plot$HRD_prob.plot <- ann_tcga.plot$HRD_prob+0.005
ann_tcga.plot <- ann_tcga.plot[order(ann_tcga.plot$HRD_prob), ]
ann_tcga.plot$index <- 1:nrow(ann_tcga.plot)

g_waterfall <- ggplot(ann_tcga.plot, aes(x = index, y = HRD_prob.plot, fill = BRCA_status)) +
  geom_bar(stat = 'identity') + theme_minimal() + ylab('p(HRD)') +
  scale_fill_manual(values = c('blue','red','lightgray','gold','darkgreen')) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'top') +
  geom_hline(yintercept = 0.79, col = 'red', linetype = 'dashed')
ggsave(filename = 'Figures/Figure1/TCGA_BRCAWaterfallPlot.pdf',
       plot = g_waterfall, height = 3, width = 7.5)

## Plot heatmap of results

# Reorder results
results.tcga_plot <- as.data.frame(t(results.tcga_loglik[rownames(ann_tcga), order(colnames(results.tcga_loglik))]))

# Sort colours
cols <- colorRampPalette(brewer.pal(9,'Set1'))
cols_pheno <- cols(length(unique(ann_tcga$Phenotype_Assigned)))
names(cols_pheno) <- unique(ann_tcga$Phenotype_Assigned)

ann_tcga_colours <- list(
  Phenotype_Assigned = cols_pheno,
  HRD = c('HRD' = 'black', 'HR-proficient' = 'white'),
  BRCA_status = c('BRCA1' = 'blue', 'BRCA2' = 'red',
                  'PALB2' = 'gold', 'RAD51C' = 'darkgreen','none' = 'white'),
  ER_status = c('[Not Evaluated]' = 'white', 'Indeterminate' = 'navy',
                'Negative' = 'gold', 'Positive' = 'darkgreen')
)

# Use white -> navy scale
cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)

pheatmap(results.tcga_plot, 
         show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_col = ann_tcga[,c(1,5,2:4)], 
         annotation_colors = ann_tcga_colours,
         color = cols_scale, fontsize = 8, fontsize_row = 10,
         filename = 'Figures/Figure1/TCGA_HRDclassificationHeatmapExtended.pdf')

table(ann_tcga$HRD, ann_tcga$BRCA_status != 'none')
table(ann_tcga$HRD, ann_tcga$BRCA_status %in% c('BRCA1','BRCA2'))
table(ann_tcga$HRD, ann_tcga$BRCA_status)

# Plot BRCA-defect/HRD status
ann_tcga$BRCA_status_broad <- sapply(ann_tcga$BRCA_status,
                                     function(x) ifelse(x == 'none', 'BRCA+', 'BRCA-defective'))

ann_tcga.plot <- ann_tcga[!is.na(ann_tcga$BRCA_status_broad), ] %>%
  group_by(HRD, BRCA_status_broad) %>% 
  summarise(n = n())
g_hrdBrca <- ggplot(ann_tcga.plot, aes(x = BRCA_status_broad, y = n, fill = HRD)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme(legend.position = 'top', axis.title.x = element_blank(),
        legend.title = element_blank()) +
  scale_fill_brewer(palette = 'Paired')
ggsave(filename = 'Figures/Figure1/TCGA_BRCASensitivity.pdf',
       plot = g_hrdBrca, width = 3.5, height = 3.5)

# BRCA type-specific HRD classification
hrd_brca1type <- c('HRD_APOBEC', 'HRD_ID6mid', 'HRD_ID8', 'HRD_SBS8')
hrd_brca2type <- c('HRD_ID6high')

ann_tcga$HRD_BRCAgroup <- ann_tcga$HRD
ann_tcga$HRD_BRCAgroup[ann_tcga$Phenotype_Assigned %in% hrd_brca1type] <- 'BRCA1-type HRD'
ann_tcga$HRD_BRCAgroup[ann_tcga$Phenotype_Assigned %in% hrd_brca2type] <- 'BRCA2-type HRD'
ann_tcga$HRD_BRCAgroup[ann_tcga$HRD_BRCAgroup == 'HRD'] <- 'HRD unassigned'

ann_tcga.plot2 <- ann_tcga[!is.na(ann_tcga$BRCA_status), ] %>%
  group_by(HRD_BRCAgroup, BRCA_status) %>%
  summarise(n=n())
ann_tcga.plot2$BRCA_status <- factor(ann_tcga.plot2$BRCA_status,
                                     levels = c('BRCA1','BRCA2',
                                                'RAD51C','PALB2','none'))
ann_tcga.plot2$HRD_BRCAgroup <- factor(ann_tcga.plot2$HRD_BRCAgroup,
                                       levels = c('HR-proficient','HRD unassigned',
                                                  'BRCA2-type HRD', 'BRCA1-type HRD'))

g_brcaHRDgroup <- ggplot(ann_tcga.plot2, aes(x = BRCA_status, y = n, fill = HRD_BRCAgroup)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme(legend.position = 'top', axis.title.x = element_blank(),
        legend.title = element_blank(), axis.title.y = element_blank()) +
  scale_fill_manual(values = c('gray90','gray50','red', 'blue'))
ggsave(filename = 'Figures/Supp_BRCAspecificHRD_BRCASensitivity.pdf',
       plot = g_brcaHRDgroup, width = 6, height = 4.5)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/TCGA_HRDhallmarks.R`:

```````R
#####
## Plotting hallmarks of HRD
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(wesanderson)

# Load TCGA HRD classifications and hallmarks, and combine
load('Results/TCGA_HRDclassification_BRCAannotation.Rdata')
tcga_HRDHallmarks <- read.table('~/Data/TCGA/HRD_hallmarks.txt', h=T, sep='\t')

resultsHRD <- merge(x = ann_tcga[,c('Patient', 'Phenotype_Assigned', 'HRD')], 
                    y = tcga_HRDHallmarks, all.x = TRUE)
resultsHRD$HRD <- factor(resultsHRD$HRD,
                         levels = c('HR-proficient', 'HRD'))

# # Plot 4 main HRD hallmarks
# resultsHRD.plot <- resultsHRD[,c('Patient','HRD','HRD_index',
#                                  'CX3','log_POLQ_FPKM','ProliferativeCapacity')] %>%
#   pivot_longer(cols = -c(Patient,HRD), names_to = 'Hallmark')
# resultsHRD.plot$Hallmark <- factor(resultsHRD.plot$Hallmark,
#                                    levels = c('HRD_index', 'CX3',
#                                               'log_POLQ_FPKM','ProliferativeCapacity'))
# g_hallmarks <- ggboxplot(data = resultsHRD.plot, x = 'HRD', y = 'value', fill = 'HRD') + 
#   stat_compare_means() + 
#   scale_fill_manual(values = wes_palette('GrandBudapest1')) +
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
#         legend.position = 'top', legend.title = element_blank(),
#         axis.text.x = element_blank()) +
#   facet_wrap(~Hallmark, scales = 'free', nrow = 1)
# ggsave(filename = 'Figures/Figure2/HRDHallmarks.pdf', height = 4)

# MYC amplification
resultsHRD$MYC_status <- factor(resultsHRD$MYC_status, 
                                levels = c('Deletion', 'Normal','Amplification'))
res_MYC <- resultsHRD %>%
  group_by(HRD, MYC_status) %>% summarise(n = n()) %>% drop_na()
res_MYC

g_myc <- ggplot(data = res_MYC, aes(x = HRD, y = n, fill = MYC_status)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme_minimal() + ylab('% Samples') +
  theme(legend.position = 'top', axis.title.y = element_blank()) + coord_flip() +
  scale_fill_manual(values = c(wes_palette('Zissou1')[3], 'gray90', wes_palette('Zissou1')[1]))
ggsave('Figures/Figure2/HRD_vs_MYCamplification.pdf', g_myc,
       width = 6, height = 3)

# POLQ expression
g_polq <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'log_POLQ_FPKM', fill = 'HRD') + 
  stat_compare_means() + 
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank()) 
ggsave('Figures/Figure2/HRD_vs_POLQexpression.pdf', g_polq,
       width = 4, height = 4)

# HRD index scores
g_hrdIndex <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'HRD_index', fill = 'HRD') + 
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank())
ggsave('Figures/Figure2/HRD_vs_HRDscore.pdf', g_hrdIndex,
       width = 4, height = 4)

# Individual HRD index scores
df_hrdIndex <- resultsHRD[,c('HRD', 'NtAI', 'LST', 'HRD.LOH')] %>%
  pivot_longer(cols = -HRD, names_to = 'HRD_index', values_to = 'score')
g_hrdIndex_individual <- ggboxplot(data = df_hrdIndex, x = 'HRD', y = 'score', fill = 'HRD') +
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~HRD_index)
ggsave('Figures/Supp_HRDvsHRDscore_ind.pdf', g_hrdIndex_individual,
       width = 8, height = 4)

# CX3 Copy Number Signature Exposure
g_cx3 <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'CX3', fill = 'HRD') + 
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank())
ggsave('Figures/Figure2/HRD_vs_CX3Exposure.pdf', g_cx3,
       width = 4, height = 4)

# # CX2/5 Copy Number Signature Exposure
# resultsHRD_CXs <- resultsHRD[,c('Patient','HRD','CX2','CX5')] %>%
#   pivot_longer(cols = c(CX2,CX5), names_to = 'CX_Signature', values_to = 'Exposure')
# g_cxs <- ggboxplot(data = resultsHRD_CXs, x = 'HRD', y = 'Exposure', fill = 'HRD') +
#   stat_compare_means() +
#   scale_fill_manual(values = wes_palette('GrandBudapest1')) +
#   theme(axis.title.x = element_blank()) +
#   facet_wrap(~CX_Signature)
# ggsave('Figures/Supplementary/HRD_vs_CXExposure.pdf', g_cxs,
#        width = 5, height = 3.5)

# Quiescence Score Comparison
g_prol <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'ProliferativeCapacity', fill = 'HRD') +
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank())
ggsave('Figures/Figure2/HRD_vs_Quiescence.pdf', g_prol,
       width = 4, height = 4)


# Breast Cancer Subtypes
res.BRCA <- resultsHRD[,c('Patient', 'HRD', 'er_status_by_ihc',
                          'her2_status_by_ihc', 'pr_status_by_ihc')]

res.BRCA$BRCA_subtype <- NA
res.BRCA$BRCA_subtype[res.BRCA$er_status_by_ihc == 'Positive'] <- 'ER+'
res.BRCA$BRCA_subtype[res.BRCA$her2_status_by_ihc == 'Positive'] <- 'HER2+'
res.BRCA$BRCA_subtype[res.BRCA$er_status_by_ihc == 'Positive' &
                        res.BRCA$her2_status_by_ihc == 'Positive'] <- 'ER+HER2+'
res.BRCA$BRCA_subtype[res.BRCA$er_status_by_ihc == 'Negative' &
                        res.BRCA$her2_status_by_ihc == 'Negative' &
                        res.BRCA$pr_status_by_ihc == 'Negative'] <- 'TNBC'

res.BRCA.plot <- res.BRCA %>%
  drop_na(BRCA_subtype) %>%
  group_by(HRD, BRCA_subtype) %>% summarise(n=n())

g_resBRCA <- ggplot(res.BRCA.plot, aes(x = BRCA_subtype, y = n, fill = HRD)) +
  geom_bar(stat = 'identity', position = 'fill') + theme_minimal() +
  xlab('') + ylab('% Samples') + theme(legend.position = 'top') +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = 'Figures/Figure2/HRDvsBRCAsubtype.pdf',
       plot = g_resBRCA, width = 4, height = 4)

g_resBRCA2 <- ggplot(res.BRCA.plot, aes(x = HRD, y = n, fill = BRCA_subtype)) +
  geom_bar(stat = 'identity', position = 'fill') + theme_minimal() +
  xlab('') + ylab('% Samples') + theme(legend.position = 'top') + 
  coord_flip() +
  scale_fill_manual(values = wes_palette('Royal2')[4:1])
ggsave(filename = 'Figures/Figure2/BRCAsubtypevsHRD.pdf',
       plot = g_resBRCA2, width = 6, height = 3)

## Plot supplementary figures

# 1. HRD hallmarks: HR-proficient vs HRD_BRCA+ vs HRD_BRCA-
ann_tcga2 <- ann_tcga
ann_tcga2 <- ann_tcga2[!is.na(ann_tcga$BRCA_status), ]
ann_tcga2$group <- ann_tcga2$HRD
ann_tcga2$group[ann_tcga2$HRD == 'HRD' & 
                  ann_tcga2$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga2$group[ann_tcga2$HRD == 'HRD' & 
                  ann_tcga2$BRCA_status != 'none'] <- 'HRD_BRCA-'
ann_tcga2$group <- factor(ann_tcga2$group,
                          levels = c('HR-proficient','HRD_BRCA+','HRD_BRCA-'))

df_supp1 <- merge(x = ann_tcga2[,c('Patient','group')],
                  y = tcga_HRDHallmarks[,c('Patient','HRD_index','CX3',
                                           'log_POLQ_FPKM','ProliferativeCapacity')])
df_supp1 <- df_supp1 %>%
  pivot_longer(cols = -c(Patient,group), names_to = 'Hallmark', values_to = 'score')

comps_supp1 <- list(c('HR-proficient','HRD_BRCA+'),
                    c('HRD_BRCA+','HRD_BRCA-'),
                    c('HR-proficient','HRD_BRCA-'))

g_brcaPos <- ggboxplot(df_supp1, x = 'group', y = 'score', fill = 'group') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top') +
  stat_compare_means(comparisons = comps_supp1) +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  facet_wrap(~Hallmark, scales = 'free', nrow=1)

ggsave(filename = 'Figures/Supp_TCGAhallmarksBRCAPositive.pdf',
       plot = g_brcaPos, width = 10, height = 5)


# 1. HRD hallmarks: HR-proficient vs p(HRD) > 0.5 vs p(HRD) > 0.79
ann_tcga3 <- ann_tcga[,c(1,3,4)]
ann_tcga3$group <- sapply(ann_tcga3$HRD_prob, function(x)
  ifelse(x > 0.79, 'p(HRD) > 0.79',
         ifelse(x > 0.5, 'p(HRD) > 0.5', 'HR-proficient')))
ann_tcga3$group <- factor(ann_tcga3$group,
                          levels = c('HR-proficient',
                                     'p(HRD) > 0.5', 'p(HRD) > 0.79'))

df_supp2 <- merge(x = ann_tcga3[,c('Patient','group')],
                  y = tcga_HRDHallmarks[,c('Patient','HRD_index','CX3',
                                           'log_POLQ_FPKM','ProliferativeCapacity')])
df_supp2 <- df_supp2 %>%
  pivot_longer(cols = -c(Patient,group), names_to = 'Hallmark', values_to = 'score')

comps_supp2 <- list(c('HR-proficient','p(HRD) > 0.5'),
                    c('p(HRD) > 0.5','p(HRD) > 0.79'),
                    c('HR-proficient','p(HRD) > 0.79'))

g_hrdThresholds <- ggboxplot(df_supp2, x = 'group', y = 'score', fill = 'group') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top') +
  stat_compare_means(comparisons = comps_supp2) +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  facet_wrap(~Hallmark, scales = 'free', nrow=1)

ggsave(filename = 'Figures/Supp_TCGAhallmarksHRDthresholds.pdf',
       plot = g_hrdThresholds, width = 10, height = 5)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/TCGA_geneCNAenrichment_HRDonly.R`:

```````R
#####
## Fisher's tests to determine differential chromosome Gene alterations in HRD vs HR-proficient samples
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(ggplot2)
library(wesanderson)
library(ggrepel)

# Load Gene-level CNA data, and organise to remove metadata
cna <- read.delim('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_cna.txt')
cna <- cna[!duplicated(cna$Hugo_Symbol), ]

# Extract cancer genes with alterations in >5% cases
cgc_data <- read.csv('~/Data/Census_allMon Jul  3 15_56_40 2023.csv')
cgc_data.genes <- as.character(cgc_data$Gene.Symbol)

cna <- cna[cna$Hugo_Symbol %in% cgc_data.genes, ]

# Sort dataframe
rownames(cna) <- cna$Hugo_Symbol
cna <- cna[,-c(1,2)]
cna <- as.data.frame(t(cna))
cna$Patient <- sapply(rownames(cna), function(x)
  paste(strsplit(x,split='[.]')[[1]][1:3], collapse='-'))

# Load HRD/HR-proficiency labels and merge with CNA data
load('Results/TCGA_HRDclassification_BRCAannotation.Rdata')
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[ann_tcga$HRD == 'HRD', ]
ann_tcga$group <- sapply(ann_tcga$BRCA_status,
                         function(x) ifelse(x == 'none', 'BRCA+', 'BRCA-defective'))

df.ann <- data.frame(
  Patient = ann_tcga$Patient,
  group = ann_tcga$group
)

df <- merge(x = df.ann, y = cna)

# Initialise data frames to track gain/loss enrichments
fishers.gain = fishers.loss <- data.frame(
  Gene = names(df)[-c(1,2)], Estimate = NA, pVal = NA
)

# For each chromosome Gene:
#   Conduct a Fisher's exact test for HRD and HR-proficient samples against:
#     1) Gain vs not-Gain
#     2) Loss vs not-Loss
#   Save odds ratios and p-values in initialised datasets
for (i in 1:nrow(fishers.gain)) {
  Gene.i <- fishers.gain$Gene[i]
  
  df.i <- data.frame(
    group = df$group,
    Gain = df[,Gene.i] > 0,
    Loss = df[,Gene.i] < 0
  )
  df.i <- df.i[!is.na(df.i[,'Gain']), ]
  
  df.i$group <- factor(df.i$group, levels = c('BRCA-defective', 'BRCA+'))
  df.i$Gain <- factor(df.i$Gain, levels = c(FALSE,TRUE))
  df.i$Loss <- factor(df.i$Loss, levels = c(FALSE,TRUE))
  
  # Gains
  table.gain.i <- table(df.i$group, df.i$Gain)
  fishers.gain.i <- fisher.test(table.gain.i)
  
  fishers.gain$Estimate[i] <- fishers.gain.i$estimate
  fishers.gain$pVal[i] <- fishers.gain.i$p.value
  
  # Losses
  table.loss.i <- table(df.i$group, df.i$Loss)
  fishers.loss.i <- fisher.test(table.loss.i)
  
  fishers.loss$Estimate[i] <- fishers.loss.i$estimate
  fishers.loss$pVal[i] <- fishers.loss.i$p.value
  
}

## Add information to gain/loss results:
#   - Adjust p-values
#   - If Estimate == 'Inf', set to maximum estimate
#   - If Estimate == 0, set to minimum estimate
#   - log2-normalise estimates and -log10(p-adjust)
#   - Add labels for Genes with enrichment significance < .01
#   - 'Volcano plot' displaying results

# Gains
fishers.gain$padj <- p.adjust(fishers.gain$pVal)
# fishers.gain$padj <- fishers.gain$pVal
fishers.gain$Estimate[fishers.gain$Estimate == 'Inf'] <- max(fishers.gain$Estimate)
fishers.gain$Estimate[fishers.gain$Estimate == 0] <- min(fishers.gain$Estimate)
fishers.gain$l2fc <- log2(fishers.gain$Estimate)
fishers.gain$logp <- -log10(fishers.gain$padj)
fishers.gain$label <- fishers.gain$Gene
fishers.gain$label[fishers.gain$padj > .05] <- NA

g_gains <- ggplot(fishers.gain, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Gain vs !Gain')
# ggsave(filename = 'Figures/Supplementary/ChrGeneEnrich_GainVsNotGain.pdf',
#        plot = g_gains)

# Losses
fishers.loss$padj <- p.adjust(fishers.loss$pVal)
# fishers.loss$padj <- fishers.loss$pVal
fishers.loss$Estimate[fishers.loss$Estimate == 'Inf'] <- max(fishers.loss$Estimate)
fishers.loss$Estimate[fishers.loss$Estimate == 0] <- min(fishers.loss$Estimate)
fishers.loss$l2fc <- log2(fishers.loss$Estimate)
fishers.loss$logp <- -log10(fishers.loss$padj)
fishers.loss$label <- fishers.loss$Gene
fishers.loss$label[fishers.loss$padj > .05] <- NA

g_losses <- ggplot(fishers.loss, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Loss vs !Loss')
# ggsave(filename = 'Figures/Supplementary/ChrGeneEnrich_LossVsNotLoss.pdf',
#        plot = g_losses)

# Merge log-fold changes (+ labels, showing significant Genes for each test)
df.fishers <- merge(x = fishers.gain[,c('Gene','l2fc','label')],
                    y = fishers.loss[,c('Gene','l2fc','label')], by = 'Gene')
names(df.fishers) <- c('Gene','l2fc_GAIN', 'label_GAIN','l2fc_LOSS','label_LOSS')

# Add GAIN/LOSS/GAIN+LOSS labels for which tests are significant
df.fishers$Label <- 'none'
df.fishers$Label[!is.na(df.fishers$label_GAIN)] <- 'GAIN'
df.fishers$Label[!is.na(df.fishers$label_LOSS)] <- 'LOSS'
df.fishers$Label[!is.na(df.fishers$label_GAIN) & !is.na(df.fishers$label_LOSS)] <- 'GAIN+LOSS'
df.fishers$Gene_Label <- NA
df.fishers$Gene_Label[df.fishers$Label != 'none'] <- df.fishers$Gene[df.fishers$Label != 'none']

# df.fishers$Label <- factor(df.fishers$Label, levels = c('GAIN','GAIN+LOSS','LOSS','none'))

g_GainvsLoss <- ggplot(df.fishers, aes(x = l2fc_GAIN, y = l2fc_LOSS, color = Label, label = Gene_Label)) +
  geom_point() + geom_text_repel() + theme_minimal() +
  scale_color_manual(values = c('gray90')) +
  xlab('log2(Fold Change) - Gain') + ylab('log2(Fold Change) - Loss') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'top',
        legend.title = element_blank())
ggsave(filename = 'Figures/Supp_GeneCnaEnrich_GainVsLoss_HRDonly.pdf',
       plot = g_GainvsLoss, height = 4, width = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/TCGA_compareHRDclassifiers.R`:

```````R
#####
## Compare HRD exome classifiers
#####

setwd('~/Projects/HRD_MutationalSignature/Results')

# Load libraries
library(tidyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(pROC)

# 1. My classifier
load('TCGA_HRD_resultsSummary.Rdata')
# results.tcga_df <- data.frame(
#   Tumor = results.tcga_df$Patient,
#   HRD79 = ifelse(results.tcga_df$HRD_prob >= 0.79,'HRD','HR-proficient'),
#   HRD50 = ifelse(results.tcga_df$HRD_prob >= 0.50,'HRD','HR-proficient'),
#   HRD32 = ifelse(results.tcga_df$HRD_prob >= 0.32,'HRD','HR-proficient')
# )
results.tcga_df <- data.frame(
  Tumor = results.tcga_df$Patient,
  p_HRDi = results.tcga_df$HRD_prob,
  HRDi = results.tcga_df$HRD
)

# 2. SigMA (run using webtool)
sigMA <- read.csv('SigMA_TCGA_output.csv')
sigMA <- data.frame(
  Tumor = sapply(sigMA$tumor, function(x) substr(x, 1, 12)),
  p_HRD_SigMA = sigMA$Signature_3_mva,
  SigMA = sigMA$pass_mva
)

# 3. deconstructSigs (run manually)
SBS3 <- read.table('TCGA_BRCA_deconstructSigs.txt', header = TRUE)
SBS3 <- data.frame(
  Tumor = sapply(rownames(SBS3), function(x) substr(x, 1, 12)),
  p_HRD_SBS3 = SBS3$SBS3,
  SBS3_dominant = apply(SBS3, 1, function(x) names(x)[x==max(x)]) == 'SBS3'
)

# 4. CX3 Copy Number Signatures (taken from Drews et al. Nature 2022)
mac.full <- read.delim('~/Data/TCGA/Macintyre_TCGAExposures.txt', sep = '\t')
mac.full <- mac.full[mac.full$Cancer == 'BRCA',c('Sample','Signature','Exposure')] %>% 
  pivot_wider(names_from = Signature, values_from = Exposure)
mac.full <- data.frame(
  Tumor = mac.full$Sample,
  p_HRD_CX3 = mac.full$CX3,
  CX3_dominant = apply(mac.full[,-1], 1, function(x) names(x)[x==max(x)]) == 'CX3'
)
names(mac.full) <- c('Tumor','p_HRD_CX3','CX3')

# 5. HRD index scores (taken from Marquard et al. 2015)
marq <- read.table('~/Data/TCGA/Marquard_HRDScores.txt', h=T)
marq$HRD_index <- apply(marq[,c('NtAI','LST','HRD.LOH')], 1, sum)
marq <- data.frame(
  Tumor = marq$Tumor,
  HRD_index = marq$HRD_index,
  HRDindex_42 = marq$HRD_index >= 42,
  HRDindex_63 = marq$HRD_index > 63
)

# Join all
df.full <- merge(x = results.tcga_df, y = sigMA, all.x = TRUE)
df.full <- merge(x = df.full, y = SBS3, all.x = TRUE)
df.full <- merge(x = df.full, y = mac.full, all.x = TRUE)
df.full <- merge(x = df.full, y = marq, all.x = TRUE)

# Add valieris BRCA defects
valieris <- read_excel('~/Data/TCGA/cancers-12-03687-s001/BRCA.xlsx', sheet = 'class-original')
valieris <- valieris[valieris$BRCA1_somatic_null != 'NA', ]
valieris <- valieris[,c('sample','event.BRCA1','event.BRCA2','event.RAD51C','event.PALB2')]
valieris$BRCA1 <- valieris$event.BRCA1 != 0
valieris$BRCA2 <- !(valieris$event.BRCA2 %in% c(0,'Mono-allelic-inactivation'))
valieris$RAD51C <- valieris$event.RAD51C != 0
valieris$PALB2 <- valieris$event.PALB2 != 0

valieris$BRCA_status <- 'none'
valieris$BRCA_status[valieris$PALB2] <- 'PALB2'
valieris$BRCA_status[valieris$RAD51C] <- 'RAD51C'
valieris$BRCA_status[valieris$BRCA2] <- 'BRCA2'
valieris$BRCA_status[valieris$BRCA1] <- 'BRCA1'

valieris <- data.frame(
  Tumor = valieris$sample,
  BRCA_defective = valieris$BRCA_status != 'none'
)

df.full <- merge(x = df.full, y = valieris, all.x = TRUE)

# Calculate sensitivities, specificities, and F-scores
break_table.func <- function(column) {
  t <- table(df.full[,column], df.full$BRCA_defective)
  
  true_negative = t[1,1]
  false_negative = t[1,2]
  false_positive = t[2,1]
  true_positive = t[2,2]
  
  return(c(true_positive, true_negative, false_positive, false_negative))
}

results <- data.frame(
  Method = c('HRDi','SigMA','SBS3_dominant','CX3','HRDindex_42','HRDindex_63'),
  TP = NA, TN = NA, FP = NA, FN = NA
)

# For each method: calculate total TP, TN, FP, FN 
for (i in 1:length(results$Method)) {
  method = results$Method[i]
  results[i,2:5] <- break_table.func(method)
}

# Calculate sensitivity, specificity, and balanced F-score of each method
weight = 1

results$recall = results$TP/(results$TP + results$FN)
results$precision = results$TP/(results$TP + results$FP)
results$F_score = (1+weight^2)*results$precision*results$recall/((weight^2)*results$precision + results$recall)
results$Method <- factor(results$Method, levels = results$Method[order(results$F_score, decreasing = TRUE)])

# Plot results
g_Fscore <- ggplot(results, aes(x = Method, y = F_score, fill = Method)) + 
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_fill_brewer(palette = 'Paired')
# ggsave(filename = '../Figures/Figure1/TCGA_HRDcomparisonsFscores.pdf', plot = g_Fscore,
#        width = 4.5, height = 4.5)

# Plot sensitivity and specificity
res.plot <- results[,c('Method','recall','precision')] %>%
  pivot_longer(-Method, names_to = 'Measure', values_to = 'value')
g_SensSpec <- ggplot(res.plot, aes(x = Method, y = value, fill = Method)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_fill_brewer(palette = 'Paired') +
  facet_wrap(~Measure, ncol = 1)

# ggsave(filename = '../Figures/Supp_TCGAcomparisonsSensSpec.pdf', plot = g_SensSpec,
#        width = 6, height = 3.6)
ggsave(filename = '~/Projects/Thesis/Chapter 3/TCGAcomparisons_SensSpec.pdf', plot = g_SensSpec,
       width = 4.5, height = 4.5)

# Plot comparative AUC curves for respective measures
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$p_HRDi), print.auc = TRUE)
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$p_HRD_SigMA), print.auc = TRUE,
                col = 'blue', print.auc.y = .45, add = TRUE)
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$p_HRD_SBS3), print.auc = TRUE,
                col = 'darkgreen', print.auc.y = .4, add = TRUE)
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$p_HRD_CX3), print.auc = TRUE,
                col = 'red', print.auc.y = .35, add = TRUE)
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$HRD_index), print.auc = TRUE,
                col = 'purple', print.auc.y = .3, add = TRUE)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ExomeHRDClassification/TCGA_HRDclassification_thresholds.R`:

```````R
#####
## Try different HRD thresholds to maximise F-score
#####

# Load libaries
library(tidyr)
library(ggplot2)
library(pROC)

# Function to calculate sensitivity/specificity/F-score

f_score <- function(t1, weight = 1) {
  recall <- t1['HRD','defective']/sum(t1[,'defective'])
  precision <- t1['HRD','defective']/sum(t1['HRD',])
  f_score <- (1 + weight^2)*(precision*recall)/((weight^2)*precision+recall)
  return(c(recall,precision,f_score))
}

hrd_thresholds <- data.frame()

# Load TCGA assignments with BRCA-defect and ER status labels
load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')

# In incremenets of 0.01, calculate F-score and add to the hrd_thresholds data frame
for (i in seq(0,1,.01)) {
  tcga_results <- ann_tcga
  tcga_results$HR_defect <- factor(ifelse(tcga_results$BRCA_status != 'none',
                                   'defective','proficient'), levels = c('defective','proficient'))
  tcga_results$HRD <- factor(ifelse(tcga_results$HRD_prob > i,
                                    'HRD','HR-proficient'), levels = c('HRD','HR-proficient'))
  
  t.i <- table(tcga_results$HRD, tcga_results$HR_defect)
  f_scores.i <- c(i, f_score(t.i,weight = 1))
  
  hrd_thresholds <- rbind(hrd_thresholds, f_scores.i)
}

names(hrd_thresholds) <- c('p_HRD', 'recall', 'precision', 'F-score')

# Plot results
hrd_thresholds.plot <- hrd_thresholds %>%
  pivot_longer(cols = -p_HRD, values_to = 'value', names_to = 'measure')
g_thres <- ggplot(hrd_thresholds.plot, aes(x = p_HRD, y = value, colour = measure)) +
  geom_line() + theme_minimal() +
  geom_vline(xintercept = hrd_thresholds$p_HRD[hrd_thresholds$`F-score` == max(hrd_thresholds$`F-score`, na.rm = TRUE)][1], 
             linetype = 'dashed') +
  xlab('p(HRD)') + theme(axis.title.y = element_blank())
ggsave(filename = '~/Projects/HRD_MutationalSignature/Figures/Supp_TCGAoptimalHRDthreshold.pdf',
       plot = g_thres, width = 5, height = 4)

# Plot AUC curve here
ann_tcga.auc <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga.auc$HR_geneDefective <- ann_tcga.auc$BRCA_status != 'none'
pdf('~/Projects/HRD_MutationalSignature/Figures/Supp_TCGAHRDthresholdAUC.pdf',
    width = 5, height = 4)
pROC_HRD <- roc(HR_geneDefective ~ HRD_prob, data = ann_tcga.auc,
                plot = TRUE, print.auc = TRUE)
dev.off()

# Repeat above analysis with subsetting for Positive and Negative ER status
er.sub <- 'Positive'
ann_tcga.er <- ann_tcga[ann_tcga$ER_status == er.sub, ]

hrd_thresholds.er <- data.frame()
for (i in seq(0,1,.01)) {
  tcga_results <- ann_tcga.er
  tcga_results$HR_defect <- factor(ifelse(tcga_results$BRCA_status != 'none',
                                          'defective','proficient'), levels = c('defective','proficient'))
  tcga_results$HRD <- factor(ifelse(tcga_results$HRD_prob > i,
                                    'HRD','HR-proficient'), levels = c('HRD','HR-proficient'))
  
  t.i <- table(tcga_results$HRD, tcga_results$HR_defect)
  f_scores.i <- c(i, f_score(t.i,weight = 1))
  
  hrd_thresholds.er <- rbind(hrd_thresholds.er, f_scores.i)
}

names(hrd_thresholds.er) <- c('p_HRD', 'recall', 'precision', 'F-score')

hrd_thresholds.er.plot <- hrd_thresholds.er %>%
  pivot_longer(cols = -p_HRD, values_to = 'value', names_to = 'measure')
g_thres.er <- ggplot(hrd_thresholds.er.plot, aes(x = p_HRD, y = value, colour = measure)) +
  geom_line() + theme_minimal() +
  geom_vline(xintercept = hrd_thresholds$p_HRD[hrd_thresholds$`F-score` == max(hrd_thresholds$`F-score`, na.rm = TRUE)][1], 
             linetype = 'dashed') +
  xlab('p(HRD)') + theme(axis.title.y = element_blank()) +
  ggtitle(paste0('ER-',er.sub, ' samples only'))
ggsave(filename = paste0('~/Projects/HRD_MutationalSignature/Figures/Supp_TCGAoptimalHRDthreshold_',er.sub,'.pdf'),
       plot = g_thres.er, width = 5, height = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ICGC_PhenotypeClusteringAndSimulationAnalysis/Likelihood_Cluster_Generation.R`:

```````R
#####
## Generate prior distributions for each signature phenotype cluster
#####

setwd('~/Projects/HRD_MutationalSignature/Results/')

# Load libraries
library(tidyr)
library(deconstructSigs) # for context ordering
library(ggplot2)

# 1. Load mutation tallies
load('~/Projects/HRD_MutationalSignature/Data/BRCA_UKEU_mt_tally.Rdata')
mut_complete <- as.data.frame(cbind(mt_tally.brca_wgs$SBS_96, mt_tally.brca_wgs$ID_83))

# 2. Load annotation (ann)
load('ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.Rdata')
mut_complete <- mut_complete[rownames(ann), ]

# Separate mutation counts into clusters
mut_byClust <- list()
for (pheno in unique(ann$Phenotype)) {
  mut_byClust[[pheno]] <- mut_complete[ann$Phenotype == pheno, ]
}

# For each input, generate a complete probability distribution, and collate them
#   NB, whilst summing seemed nice, SigMA used a mean. Let's try both!

collate_function <- function(input, variant_type = 'ALL', collation = 'total') {
  
  sbs.index <- sapply(names(input), function(x) grepl('>',x))
  id.index <- sapply(names(input), function(x) grepl(':',x))
  
  if (variant_type == 'SBS') {
    input_final <- input[,sbs.index]
  } else if (variant_type == 'ID') {
    input_final <- input[,id.index]
  } else {
    input_final <- input
  }
  
  # Add pseudocount
  input_final <- input_final + 1 # original
  if (collation == 'mean') {
    dist_temp = t(apply(input_final, 1, function(x) x/sum(x)))
    dist = apply(dist_temp, 2, mean)
  } else {
    if (collation != 'total') print('Set collation to mean or total. Setting to total...')
    dist_temp = apply(input_final, 2, sum)
    dist = dist_temp/sum(dist_temp)
  }
  
  return(dist)
}

mut.dists_total = mut.dists_mean <- data.frame()

for (pheno in names(mut_byClust)) {
  mut.dists_total <- rbind(mut.dists_total, collate_function(input = mut_byClust[[pheno]]))
  mut.dists_mean <- rbind(mut.dists_mean, collate_function(mut_byClust[[pheno]], collation = 'mean'))
}

colnames(mut.dists_total) = colnames(mut.dists_mean) <- colnames(mut_byClust[[1]])
rownames(mut.dists_total) = rownames(mut.dists_mean) <- names(mut_byClust)

# Save prior clusters
save(mut.dists_mean, file = '../Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.Rdata')
save(mut.dists_total, file = '../Data/ClusterLikelihoods/ICGC_clust20_mclust_totalCont.Rdata')

# Plot prior clusters
mut.dists_mean.plot <- cbind(mut.dists_mean[,colnames(signatures.cosmic)],
                             mut.dists_mean[,97:ncol(mut.dists_mean)])
mut.dists_mean.plot$Phenotype <- rownames(mut.dists_mean)
mut.dists_mean.plot <- mut.dists_mean.plot %>%
  pivot_longer(cols = -Phenotype, names_to = 'Context', values_to = 'Contribution')
mut.dists_mean.plot$Phenotype <- factor(mut.dists_mean.plot$Phenotype,
                                        levels = sort(rownames(mut.dists_mean)))

# Sort indel context order
signatures.id83 <- read.table('~/Data/COSMIC_v3.3_ID_GRCh37.txt', h=T)
mut.dists_mean.plot$Context <- factor(mut.dists_mean.plot$Context,
                                      levels = c(colnames(signatures.cosmic), signatures.id83$Type))
mut.dists_mean.plot$Type <- ifelse(
  grepl(pattern = '>', mut.dists_mean.plot$Context),
  'SBS','indel')

g_meanPlot <- ggplot(mut.dists_mean.plot, aes(x = Context, y = Contribution, fill = Type)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = 'top') +
  facet_wrap(~Phenotype, scales = 'free')
ggsave(filename = '../Figures/SupplementaryFigures/Supp_LikelihoodDistributionsMeans.pdf',
       plot = g_meanPlot, height = 5, width = 10)



```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ICGC_PhenotypeClusteringAndSimulationAnalysis/ICGC_deconstructSigs_genome2exome.R`:

```````R
#####
## Contribution of BRCA-associated signatures in ICGC using deconstructSigs
#####

# Load libraries
library(deconstructSigs)

# Load ICGC-BRCA signature contributions
#   Only signatures appearing in >1% samples will be included
brca_final.sigs <- read.table('~/Data/ICGC/ICGC_BRCA_sigProfilerCont.txt', header = TRUE)
brca_final.sigs <- brca_final.sigs[brca_final.sigs$Prop_present > 0.01, ]
sigs.sbs <- brca_final.sigs$Sigs[grepl(pattern = 'SBS', brca_final.sigs$Sigs)]
sigs.id <- brca_final.sigs$Sigs[grepl(pattern = 'ID', brca_final.sigs$Sigs)]

# Load hg19 signatures (for ICGC samples)
signatures.sbs96_hg19 <- read.table('~/Data/COSMIC_v3.3.1_SBS_GRCh37.txt', h=T)
sig_sbs96_types <- signatures.sbs96_hg19$Type
signatures.sbs96_hg19 <- signatures.sbs96_hg19[, 2:ncol(signatures.sbs96_hg19)]
signatures.sbs96_hg19 <- data.frame(t(signatures.sbs96_hg19))
colnames(signatures.sbs96_hg19) <- sig_sbs96_types
signatures.sbs96_hg19 <- signatures.sbs96_hg19[, colnames(signatures.cosmic)]

# Load COSMIC indel signatures
signatures.id83 <- read.table('~/Data/COSMIC_v3.3_ID_GRCh37.txt', h=T)
sig_id83_types <- signatures.id83$Type
signatures.id83 <- signatures.id83[, 2:ncol(signatures.id83)]
signatures.id83 <- data.frame(t(signatures.id83))
colnames(signatures.id83) <- sig_id83_types

# Load in datasets and tweak into deconstructSigs inputs:
#   rows = samples, cols = signature contexts
#   order the same as signature data frames
load('~/Projects/HRD_MutationalSignature/Data/BRCA_UKEU_mt_tally.Rdata')
sigs.sbs96_input <- as.data.frame(mt_tally.brca_wgs$SBS_96)
sigs.sbs96_input <- sigs.sbs96_input[,colnames(signatures.sbs96_hg19)]

sigs.id83_input <- as.data.frame(mt_tally.brca_wgs$ID_83)
sigs.id83_input <- sigs.id83_input[,colnames(signatures.id83)]

# # Remove samples with fewer than 50 indels
# index.lowIDcount <- which(apply(sigs.id83_input, 1, sum) < 30)
# sigs.sbs96_input <- sigs.sbs96_input[-index.lowIDcount, ]
# sigs.id83_input <- sigs.id83_input[-index.lowIDcount, ]

# Normalise ID83 counts
icgc_idCounts <- read.table('~/Data/ICGC/ICGC_BRCA_indelCounts.txt')
icgc_idCounts$genome_to_exome <- icgc_idCounts$exome/icgc_idCounts$wgs
icgc_idCounts <- icgc_idCounts[colnames(sigs.id83_input), ]

sigs.id83_input <- as.data.frame(t(apply(sigs.id83_input, 1, function(x) x * icgc_idCounts$genome_to_exome)))

# Run deconstructSigs on each input tp calculate signature proportions
run_deconstructSigs <- function(sigs.input, sig_type = 'SBS') {
  
  if (sig_type == 'SBS') {
    sig_ref = signatures.sbs96_hg19[sigs.sbs, ] # only SBS signatures in ICGC-BRCA
    print('Calculating SBS signature contributions...')
  } else if (sig_type == 'ID') {
    sig_ref = signatures.id83[sigs.id, ] # only ID signatures in ICGC-BRCA
    print('Calculating ID signature contributions...')
  } else {
    print('Set sig_type to SBS or ID. Using SBS signatures...')
    sig_ref = signatures.sbs96_hg19[sigs.sbs, ]
  }
  
  sigs_out <- NULL
  for (i in 1:nrow(sigs.input)) {
    sample.i <- rownames(sigs.input)[i]
    print(paste0(sig_type, ' contributions, Sample ', i, ' of ', nrow(sigs.input), ': ', sample.i))
    
    sigs_i <- whichSignatures(
      tumor.ref = sigs.input,
      signatures.ref = sig_ref,
      sample.id = sample.i,
      contexts.needed = TRUE,
      signature.cutoff = 0,
      tri.counts.method = ifelse(sig_type == 'SBS', 'genome2exome', 'default')
    )
    sigs_out <- rbind(sigs_out, sigs_i$weights)
  }
  
  return(sigs_out)
  
}
sigs.sbs96 <- run_deconstructSigs(sigs.sbs96_input, 'SBS')
sigs.id83 <- run_deconstructSigs(sigs.id83_input, 'ID')
sigs_complete <- cbind(sigs.sbs96, sigs.id83)

save(sigs_complete,
     file = '~/Projects/HRD_MutationalSignature/Results/ICGC_BRCA_deconstructSigs_Cutoff0.01_SBSandIDnormalised.Rdata')

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ICGC_PhenotypeClusteringAndSimulationAnalysis/ICGC_BRCA_sigProfilerContributions.R`:

```````R
# Extract SBS and ID signatures in the ICGC-BRCA cohort

setwd('~/Data/ICGC/')

library(tidyr)

# Load SBS data
sbs_pcawg <- read.csv('PCAWG_sigProfiler_SBS_signatures_in_samples.csv', row.names = 2)
sbs_nopcawg <- read.csv('nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv', row.names = 2)

sbs <- rbind(sbs_pcawg, sbs_nopcawg)
sbs_brca <- sbs[grepl(pattern = 'Breast', sbs$Cancer.Types), -c(1,2)]

sbs_brca_props <- data.frame(
  Sigs = colnames(sbs_brca),
  Prop_present = apply(sbs_brca, 2, function(x) sum(x > 0)/length(x))
)

# Load ID data
id <- read.csv('PCAWG_SigProfiler_ID_signatures_in_samples.csv', row.names = 2)
id_brca <- id[grepl(pattern = 'Breast', id$Cancer.Types), -c(1,2)]

id_brca_props <- data.frame(
  Sigs = colnames(id_brca),
  Prop_present = apply(id_brca, 2, function(x) sum(x > 0)/length(x))
)

# Collate and save
brca_sigs <- rbind(sbs_brca_props, id_brca_props)
write.table(brca_sigs, file = 'ICGC_BRCA_sigProfilerCont.txt',
            row.names = FALSE, col.names = TRUE)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ICGC_PhenotypeClusteringAndSimulationAnalysis/ICGC_simulations_indelProportions.R`:

```````R
#####
## Simulation analysis to determine utility of indels in HRD classification
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(pheatmap)
library(pROC)
library(ggpubr)

# Aim: demonstrate that, in subsampled cohorts, the inclusion of indels improves classification
#   In this case, they improve the probability of reclassifying SBS3-enriched samples as SBS3
#   The ICGC cohort easily clusters into 3 based on SBS signature contributions:
#     SBS2/13 (APOBEC), SBS3 (HRD), SBS5 (Ageing)

## 1. CREATE SBS MUTATIONAL SIGNATURE CLUSTERS

# Signature contributions have already been calculated. Extract SBS signatures
load('Results/ICGC_BRCA_deconstructSigs_Cutoff0.01_SBSandIDnormalised.Rdata')
sigs_complete <- sigs_complete[, grepl(pattern = 'SBS', names(sigs_complete))]

cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)
pheatmap(sigs_complete, cutree_rows = 3,
         show_rownames = FALSE, color = cols_scale,
         file = 'Figures/Supp_ICGCSignaturesSBSHeatmap.pdf')

# Create groups using hierarchical clustering
sigs.dist_mat <- dist(sigs_complete, method = 'euclidean')
sigs.hclust3 <- hclust(sigs.dist_mat, method = 'complete')
sigs.clust3 <- cutree(sigs.hclust3, k = 3)
sigs.clust3 <- data.frame(Cluster.num = sigs.clust3)
sigs.clust3$Phenotype <- sapply(sigs.clust3$Cluster.num, function(x)
  ifelse(x == 1, 'SBS5',
         ifelse(x == 2, 'SBS3', 'APOBEC')))

# Load in ICGC mutation tallies, which have already been processed
load('Data/BRCA_UKEU_mt_tally.Rdata')
mut_complete <- cbind(mt_tally.brca_wgs$SBS_96,
                      mt_tally.brca_wgs$ID_83)

## 2. GENERATE CLUSTER-SPECIFIC MUTATIONAL SPECTRA

# Separate out mut_complete based on cluster assignment
#   Then use collate_function() to generate the SBS+ID spectrum for each

mut_apobec <- mut_complete[sigs.clust3$Phenotype == 'APOBEC', ]
mut_sbs3 <- mut_complete[sigs.clust3$Phenotype == 'SBS3', ]
mut_sbs5 <- mut_complete[sigs.clust3$Phenotype == 'SBS5', ]

collate_function <- function(input, variant_type = 'ALL') {
  
  # Spectra can be generate for SBS/ID-only, or ALL
  sbs.index <- sapply(names(input), function(x) grepl('>',x))
  id.index <- sapply(names(input), function(x) grepl(':',x))
  
  if (variant_type == 'SBS') {
    input_final <- input[,sbs.index]
  } else if (variant_type == 'ID') {
    input_final <- input[,id.index]
  } else {
    input_final <- input
  }
  
  # Return average spectrum
  total = apply(input_final, 2, sum)
  dist = total/sum(total)
  return(dist)
  
}

mut.prob_apobec <- collate_function(mut_apobec)
mut.prob_sbs3 <- collate_function(mut_sbs3)
mut.prob_sbs5 <- collate_function(mut_sbs5)

mut.prob <- rbind(mut.prob_apobec, mut.prob_sbs3, mut.prob_sbs5)
rownames(mut.prob) <- c('APOBEC','SBS3','SBS5')

## 3. LIKELIHOOD FUNCTION

# This function aligns a dataset with the designated mean distribution
likelihood_calc <- function(input_data, cluster_distributions, cluster_assign) {
  
  # This function:
  #   Takes a dataset of mutations as its input (rows = samples, cols = mutation types)
  #   Limits the data to to designated mutation types
  #   Applies a likelihood approach to calculate posterior probabilities of cluster assignment
  #   NB, For now, we ensure that our cluster distributions have the correct mutation types
  #   Output: posterior distributions, note that for WES data a comparison column is not possible
  
  # Calculate likelihood: P(mutation spectrum | cluster)
  likelihoods <- matrix(NA, nrow = nrow(input_data), ncol = nrow(cluster_distributions))
  rownames(likelihoods) <- rownames(input_data); colnames(likelihoods) <- rownames(cluster_distributions)
  for (i in 1:nrow(input_data)) {
    likelihoods[i, ] <- apply(cluster_distributions, 1, 
                              function(x) prod(x ^ input_data[i, ]))
  }
  
  # Set prior: P(cluster)
  marginal.probs <- table(cluster_assign)/length(cluster_assign)
  
  # Calculate posteriors and return
  posteriors <- data.frame(t(apply(likelihoods, 1, function(x) 
    (x * marginal.probs)/sum(x * marginal.probs))))
  return(posteriors)
  
}

## 4. INDEL-BASED SUBSAMPLING

# This function applies the likelihood function to a newly generated subsample
#   This subsample consists of sample_size mutation events extracted with replacement from the original dataset
#   Of these events, indel_prop * sample_size will be indel events, and the rest are SBS (initially indel_prop = 0)

# Apply likelihood fucntion to subsampled data
simulate_likelihood_calc <- function(input_data, cluster_distributions,
                                     cluster_assign, sample_size, indel_prop = 0) {
  
  # Separate data into SBS and ID, for separate sampling
  data_sbs <- input_data[, grepl(pattern = '>', colnames(input_data))]
  data_id  <- input_data[, grepl(pattern = ':', colnames(input_data))]
  
  # Develop simulation matrix (same dimensions as input_data)
  sims <- matrix(0, nrow = nrow(input_data), ncol = ncol(input_data))
  colnames(sims) <- colnames(input_data); rownames(sims) <- rownames(input_data)
  
  for (i in 1:nrow(input_data)) {
    s_sbs <- table(sample(colnames(data_sbs), size = sample_size * (1 - indel_prop),
                          replace = TRUE, prob = data_sbs[i, ]/sum(data_sbs[i, ])))
    sims[i, names(s_sbs)] <- s_sbs
    
    s_id <- table(sample(colnames(data_id), size = sample_size * indel_prop,
                         replace = TRUE, prob = data_id[i, ]/sum(data_id[i, ])))
    sims[i, names(s_id)] <- s_id
  }
  
  # Apply likelihood function to the simulated data
  posteriors <- likelihood_calc(sims, cluster_distributions,
                                cluster_assign)
  
  # As this is simulated data, we can also add the comparisons
  posteriors$PhenoTrue <- cluster_assign
  
  # Calculate AUCs and return
  roc_apobec <- roc(posteriors$PhenoTrue == 'APOBEC', posteriors$APOBEC, quiet = TRUE)
  roc_sbs3 <- roc(posteriors$PhenoTrue == 'SBS3', posteriors$SBS3, quiet = TRUE)
  roc_sbs5 <- roc(posteriors$PhenoTrue == 'SBS5', posteriors$SBS5, quiet = TRUE)
  
  auc.values <- c(roc_apobec$auc, roc_sbs3$auc, roc_sbs5$auc)
  
  return(list(auc.values = auc.values, posteriors = posteriors))
  
}


## 5. COMPLETE SIMULATION RUNNING

# This function runs n simulations, returning a data frame of AUCs and plottable posteriors
run_likelihood_sims <- function(input_data, sample_size, indel_prop,
                                cluster_distributions, cluster_assign,
                                n_simulations) {
  
  # Initialise data frames
  results.mat <- matrix(NA, nrow=0, ncol=3); colnames(results.mat) <- rownames(cluster_distributions)
  posterior.df <- data.frame()
  
  # Run n_simulations of simulate_likelihood_calc()
  for (i in 1:n_simulations) {
    
    set.seed(i)
    
    print(paste0('Sample size = ', sample_size,
                 ' indel_prop = ', indel_prop,
                 ' Running simulation ', i, ' of ', n_simulations, '...'))
    run.i <- simulate_likelihood_calc(
      input_data, cluster_distributions, cluster_assign,
      sample_size, indel_prop)
    results.mat <- rbind(results.mat, run.i$auc.values)
    
    run.i$posteriors$Run <- i
    posterior.df <- rbind(posterior.df, run.i$posteriors)
    
  }
  
  # Further posterior info
  posterior.df$sample_size = sample_size
  posterior.df$indel_prop = indel_prop
  
  results <- data.frame(Pheno = rep(colnames(results.mat), each = nrow(results.mat)),
                        AUC = c(results.mat[,1], results.mat[,2], results.mat[,3]),
                        sample_size = sample_size, indel_prop = indel_prop)
  
  return(list(results = results, posteriors = posterior.df))
  
}

nSim = 100
res.full_25 <- res.full_50 <- res.full_100 <- data.frame()
full_posteriors <- data.frame()

for (indel_props in seq(0, 0.5, by = .05)) {
  
  res.i_25 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = 25, indel_prop = indel_props,
    cluster_distributions = mut.prob, cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_25 <- rbind(res.full_25, res.i_25$results)
  full_posteriors <- rbind(full_posteriors, res.i_25$posteriors)
  
  print(paste0('Indel proportion = ', indel_props,
               ' sample size = 50'))
  
  res.i_50 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = 50, indel_prop = indel_props,
    cluster_distributions = mut.prob, cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_50 <- rbind(res.full_50, res.i_50$results)
  full_posteriors <- rbind(full_posteriors, res.i_50$posteriors)
  
  print(paste0('Indel proportion = ', indel_props,
               ' sample size = 100'))
  
  res.i_100 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = 100, indel_prop = indel_props,
    cluster_distributions = mut.prob, cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_100 <- rbind(res.full_100, res.i_100$results)
  full_posteriors <- rbind(full_posteriors, res.i_100$posteriors)
  
}

myComp <- list(c('0', '0.05'), c('0', '0.1'),
               c('0', '0.15'), c('0', '0.2'), c('0', '0.25'))

res.full_25$best_AUC <- res.full_25$indel_prop == .2
g_sim25 <- ggboxplot(res.full_25[res.full_25$Pheno == 'SBS3', ],
          x = 'indel_prop', y = 'AUC', add = 'jitter', color = 'best_AUC') +
  xlab('Indel Proportions') +
  geom_vline(xintercept = 1+(1/0.05)*0.0675, linetype = 'dashed', color = 'blue') +
  guides(color = 'none') + scale_color_manual(values = c('black','red')) +
  stat_compare_means(comparisons = myComp) + ggtitle('sample size = 25')
# ggsave(filename = 'Figures/Supp_ICGCsimulations_indelProps_sim25.pdf', plot = g_sim25,
#        width = 8, height = 4)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_indelProps_sim25.pdf', plot = g_sim25,
       width = 12, height = 4)

g_sim50 <- ggboxplot(res.full_50[res.full_50$Pheno == 'SBS3', ],
          x = 'indel_prop', y = 'AUC', add = 'jitter') +
  geom_vline(xintercept = 1+(1/0.05)*0.0675, linetype = 'dashed', color = 'red') +
  stat_compare_means(comparisons = myComp)
# ggsave(filename = 'Figures/Figure1/ICGCsimulations_indelProps_sim50.pdf', plot = g_sim50,
#        width = 8, height = 4)

res.full_50$best_AUC <- res.full_50$indel_prop == .1
g_sim50 <- ggboxplot(res.full_50[res.full_50$Pheno == 'SBS3', ],
                     x = 'indel_prop', y = 'AUC', add = 'jitter', color = 'best_AUC') +
  xlab('Indel Proportions') +
  geom_vline(xintercept = 1+(1/0.05)*0.0675, linetype = 'dashed', color = 'blue') +
  guides(color = 'none') + scale_color_manual(values = c('black','red')) +
  stat_compare_means(comparisons = myComp) + ggtitle('sample size = 50')
# ggsave(filename = 'Figures/Supp_ICGCsimulations_indelProps_sim50.pdf', plot = g_sim50,
#        width = 8, height = 4)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_indelProps_sim50.pdf', plot = g_sim50,
       width = 12, height = 4)

res.full_100$best_AUC <- res.full_100$indel_prop == .05
g_sim100 <- ggboxplot(res.full_100[res.full_100$Pheno == 'SBS3', ],
          x = 'indel_prop', y = 'AUC', add = 'jitter', color = 'best_AUC') +
  xlab('Indel Proportions') +
  geom_vline(xintercept = 1+(1/0.05)*0.0675, linetype = 'dashed', color = 'blue') +
  guides(color = 'none') + scale_color_manual(values = c('black','red')) +
  stat_compare_means(comparisons = myComp) + ggtitle('sample size = 100')
# ggsave(filename = 'Figures/Supp_ICGCsimulations_indelProps_sim100.pdf', plot = g_sim100,
#        width = 8, height = 4)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_indelProps_sim100.pdf', plot = g_sim100,
       width = 12, height = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ICGC_PhenotypeClusteringAndSimulationAnalysis/ICGC_BRCA_UKandEU_MatrixGeneration.R`:

```````R
#####
## Processing of ICGC WGS breast cancer samples to infer SBS and indel mutation type contributions
####

# Load libraries
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)

## Load UK and EU BRCA data (US is exome, FR contains no indels)
# Extract WGS strategies, remove duplicated mutations

brca_uk <- read.table('~/Data/ICGC/simple_somatic_mutation.open.BRCA-UK.tsv.gz',
                      sep = '\t', header = TRUE)
brca_uk_wgs <- brca_uk[brca_uk$sequencing_strategy == 'WGS', ]
brca_uk_wgs <- brca_uk_wgs[!duplicated(brca_uk_wgs$icgc_mutation_id), ]

brca_eu <- read.table('~/Data/ICGC/simple_somatic_mutation.open.BRCA-EU.tsv.gz',
                      sep = '\t', header = TRUE)
brca_eu_wgs <- brca_eu[brca_eu$sequencing_strategy == 'WGS', ]
brca_eu_wgs <- brca_eu_wgs[!duplicated(brca_eu_wgs$icgc_mutation_id), ]

# Collate ICGC data and organise into MAF-readable format

brca_wgs <- rbind(brca_uk_wgs, brca_eu_wgs)

brca_wgs_input <- data.frame(
  Tumor_Sample_Barcode = brca_wgs$icgc_donor_id,
  Hugo_Symbol = NA,
  Chromosome = brca_wgs$chromosome,
  Start_position = brca_wgs$chromosome_start,
  End_position = brca_wgs$chromosome_end,
  Variant_Classification = brca_wgs$consequence_type,
  Variant_Type = sapply(brca_wgs$mutation_type,
                        function(x) ifelse(x == 'single base substitution', 'SNP',
                                           ifelse(x == 'insertion of <=200bp', 'INS',
                                                  ifelse(x == 'deletion of <=200bp', 'DEL',NA)))),
  Reference_Allele = brca_wgs$reference_genome_allele,
  Tumor_Seq_Allele2 = brca_wgs$mutated_to_allele
)

brca_maf <- read.maf(maf = brca_wgs_input,
                     vc_nonSyn = names(table(brca_wgs_input$Variant_Classification)))

# Run mt_tally() from sigminer package to collate mutation type contributions

mt_tally.brca_wgs <- sig_tally(
  object = brca_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ALL',
  useSyn = TRUE
)

save(mt_tally.brca_wgs, 
     file = '~/Projects/HRD_MutationalSignature/Data/BRCA_UKEU_mt_tally.Rdata')

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ICGC_PhenotypeClusteringAndSimulationAnalysis/ICGC_simulations_indelLikelihoodAlterations.R`:

```````R
#####
## Simulation analysis to determine utility of indels in HRD classification
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(pheatmap)
library(pROC)
library(ggpubr)

# Aim: demonstrate that, in subsampled cohorts, the inclusion of indels improves classification
#   In this case, they improve the probability of reclassifying SBS3-enriched samples as SBS3
#   The ICGC cohort easily clusters into 3 based on SBS signature contributions:
#     SBS2/13 (APOBEC), SBS3 (HRD), SBS5 (Ageing)

## 1. CREATE SBS MUTATIONAL SIGNATURE CLUSTERS

# Signature contributions have already been calculated. Extract SBS signatures
load('Results/ICGC_BRCA_deconstructSigs_Cutoff0.01_SBSandIDnormalised.Rdata')
sigs_complete <- sigs_complete[, grepl(pattern = 'SBS', names(sigs_complete))]

pheatmap(sigs_complete, cutree_rows = 3, method = 'average',
         show_rownames = FALSE)

# Create groups using hierarchical clustering
sigs.dist_mat <- dist(sigs_complete, method = 'euclidean')
sigs.hclust3 <- hclust(sigs.dist_mat, method = 'complete')
sigs.clust3 <- cutree(sigs.hclust3, k = 3)
sigs.clust3 <- data.frame(Cluster.num = sigs.clust3)
sigs.clust3$Phenotype <- sapply(sigs.clust3$Cluster.num, function(x)
  ifelse(x == 1, 'SBS5',
         ifelse(x == 2, 'SBS3', 'APOBEC')))

# Load in ICGC mutation tallies, which have already been processed
load('Data/BRCA_UKEU_mt_tally.Rdata')
mut_complete <- cbind(mt_tally.brca_wgs$SBS_96,
                      mt_tally.brca_wgs$ID_83)

## 2. GENERATE CLUSTER-SPECIFIC MUTATIONAL SPECTRA

# Separate out mut_complete based on cluster assignment
#   Then use collate_function() to generate the SBS+ID spectrum for each

mut_apobec <- mut_complete[sigs.clust3$Phenotype == 'APOBEC', ]
mut_sbs3 <- mut_complete[sigs.clust3$Phenotype == 'SBS3', ]
mut_sbs5 <- mut_complete[sigs.clust3$Phenotype == 'SBS5', ]

collate_function <- function(input, variant_type = 'ALL') {
  
  # Spectra can be generate for SBS/ID-only, or ALL
  sbs.index <- sapply(names(input), function(x) grepl('>',x))
  id.index <- sapply(names(input), function(x) grepl(':',x))
  
  if (variant_type == 'SBS') {
    input_final <- input[,sbs.index]
  } else if (variant_type == 'ID') {
    input_final <- input[,id.index]
  } else {
    input_final <- input
  }
  
  # Return average spectrum
  total = apply(input_final, 2, sum)
  dist = total/sum(total)
  return(dist)
  
}

mut.prob_apobec <- collate_function(mut_apobec)
mut.prob_sbs3 <- collate_function(mut_sbs3)
mut.prob_sbs5 <- collate_function(mut_sbs5)

mut.prob <- rbind(mut.prob_apobec, mut.prob_sbs3, mut.prob_sbs5)
rownames(mut.prob) <- c('APOBEC','SBS3','SBS5')

## 3. LIKELIHOOD FUNCTION

# This function aligns a dataset with the designated mean distribution
likelihood_calc <- function(input_data, cluster_distributions, cluster_assign) {
  
  # This function:
  #   Takes a dataset of mutations as its input (rows = samples, cols = mutation types)
  #   Limits the data to to designated mutation types
  #   Applies a likelihood approach to calculate posterior probabilities of cluster assignment
  #   NB, For now, we ensure that our cluster distributions have the correct mutation types
  #   Output: posterior distributions, note that for WES data a comparison column is not possible
  
  # Calculate likelihood: P(mutation spectrum | cluster)
  likelihoods <- matrix(NA, nrow = nrow(input_data), ncol = nrow(cluster_distributions))
  rownames(likelihoods) <- rownames(input_data); colnames(likelihoods) <- rownames(cluster_distributions)
  for (i in 1:nrow(input_data)) {
    likelihoods[i, ] <- apply(cluster_distributions, 1, 
                              function(x) prod(x ^ input_data[i, ]))
  }
  
  # Set prior: P(cluster)
  marginal.probs <- table(cluster_assign)/length(cluster_assign)
  
  # Calculate posteriors and return
  posteriors <- data.frame(t(apply(likelihoods, 1, function(x) 
    (x * marginal.probs)/sum(x * marginal.probs))))
  return(posteriors)
  
}

## 4. WRITE SIMULATION FUNCTION

# Apply likelihood function to simulated subsampled data
simulate_likelihood_calc <- function(input_data, cluster_distributions,
                                     cluster_assign, sample_size) {
  
  # Develop simulation matrix (same dimensions as input_data)
  sims <- matrix(0, nrow = nrow(input_data), ncol = ncol(input_data))
  colnames(sims) <- colnames(input_data); rownames(sims) <- rownames(input_data)
  
  for (i in 1:nrow(input_data)) {
    s.i <- table(sample(colnames(input_data), size = sample_size,
                        replace = TRUE, prob = input_data[i,]/sum(input_data)))
    sims[i, names(s.i)] <- s.i
  }
  
  # Apply likelihood function to the simulated data
  posteriors <- likelihood_calc(sims, cluster_distributions,
                                cluster_assign)
  
  # As this is simulated data, we can also add the comparisons
  posteriors$PhenoTrue <- cluster_assign
  
  # Calculate AUCs and return
  roc_apobec <- roc(posteriors$PhenoTrue == 'APOBEC', posteriors$APOBEC, quiet = TRUE)
  roc_sbs3 <- roc(posteriors$PhenoTrue == 'SBS3', posteriors$SBS3, quiet = TRUE)
  roc_sbs5 <- roc(posteriors$PhenoTrue == 'SBS5', posteriors$SBS5, quiet = TRUE)
  
  auc.values <- c(roc_apobec$auc, roc_sbs3$auc, roc_sbs5$auc)
  
  return(list(auc.values = auc.values, posteriors = posteriors))
  
}

## 5. Save likelihood distributions with altered indel contributions

# Generate a list of likelihood distributions with indel contribution * indel_alter
#   of which indel_alter = seq(0.2, 2, by=.2)

cluster_distribution.list <- list()
indel_alterations <- c(1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5)
for (i in 1:length(indel_alterations)) {
  mut.prob_sbs <- mut.prob[,grepl(pattern = '>', colnames(mut.prob))]
  mut.prob_id <- mut.prob[, grepl(pattern = ':', colnames(mut.prob))]
  
  indel_alter <- indel_alterations[i]
  
  mut.prob_id <- mut.prob_id * indel_alter
  for (j in 1:3) {
    mut.prob_sbs[j,] <- mut.prob_sbs[j,]*(1-sum(mut.prob_id[j,]))/sum(mut.prob_sbs[j,])
  }
  
  mut.prob_altered <- cbind(mut.prob_sbs, mut.prob_id)
  cluster_distribution.list[[i]] <- mut.prob_altered
  
}

## 6. COMPLETE SIMULATION RUNNING

# This function runs n simulations, returning a data frame of AUCs and plottable posteriors
run_likelihood_sims <- function(input_data, sample_size, 
                                cluster_distributions_index, cluster_assign,
                                n_simulations) {
  
  cluster_distributions <- cluster_distribution.list[[cluster_distributions_index]]
  
  # Initialise data frames
  results.mat <- matrix(NA, nrow=0, ncol=3); colnames(results.mat) <- rownames(cluster_distributions)
  posterior.df <- data.frame()
  
  
  # Run n_simulations of simulate_likelihood_calc()
  for (i in 1:n_simulations) {
    
    set.seed(i)
    
    print(paste0('Sample size = ', sample_size,
                 ' mut.prob_index = ', indel_alterations[cluster_distributions_index],
                 ' Running simulation ', i, ' of ', n_simulations, '...'))
    
    run.i <- simulate_likelihood_calc(
      input_data, cluster_distributions, cluster_assign,
      sample_size)
    results.mat <- rbind(results.mat, run.i$auc.values)
    
    run.i$posteriors$Run <- i
    posterior.df <- rbind(posterior.df, run.i$posteriors)
    
  }
  
  # Further posterior info
  posterior.df$sample_size = sample_size
  posterior.df$indel_alter = indel_alterations[cluster_distributions_index]
  
  results <- data.frame(Pheno = rep(colnames(results.mat), each = nrow(results.mat)),
                        AUC = c(results.mat[,1], results.mat[,2], results.mat[,3]),
                        sample_size = sample_size, 
                        indel_alteration = indel_alterations[cluster_distributions_index])
  
  return(list(results = results, posteriors = posterior.df))
  
}

nSim = 100
sampleSize = 25
res.full_25 <- data.frame()

set.seed(123)
for (i in 1:length(cluster_distribution.list)) {
  
  print(paste0('Indel alteration = ', indel_alterations[i]))
  
  res.i_25 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = sampleSize, 
    cluster_distributions_index = i, 
    cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_25 <- rbind(res.full_25, res.i_25$results)

}

# Plot results for HRD Pheno
res_25.sbs3 <- res.full_25[res.full_25$Pheno == 'SBS3', ]

my_comparisons <- list(c('0.2', '1'), c('1', '5'))

g25 <- ggboxplot(res_25.sbs3, x = 'indel_alteration', y = 'AUC', add = 'jitter') +
  stat_compare_means(comparisons = my_comparisons) + ylim(c(0.65,1.05)) +
  ggtitle('sample size = 25') +
  scale_x_discrete(labels = c('1/5', '1/4', '1/3', '1/2', '1',
                              '2', '3', '4', '5'))
# ggsave(filename = 'Figures/Supp_SimulationsIndelLikelihoods_size25.pdf', plot = g25,
#        height = 5, width = 5)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_IndelWeighting_sim25.pdf', plot = g25,
       height = 7, width = 4)


# Change sample size
sampleSize = 50
res.full_50 <- data.frame()

set.seed(123)
for (i in 1:length(cluster_distribution.list)) {
  
  print(paste0('Indel alteration = ', indel_alterations[i]))
  
  res.i_50 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = sampleSize, 
    cluster_distributions_index = i, 
    cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_50 <- rbind(res.full_50, res.i_50$results)
  
}

# Plot results for HRD Pheno
res_50.sbs3 <- res.full_50[res.full_50$Pheno == 'SBS3', ]

my_comparisons <- list(c('0.2', '1'), c('1', '5'))

g50 <- ggboxplot(res_50.sbs3, x = 'indel_alteration', y = 'AUC', add = 'jitter') +
  stat_compare_means(comparisons = my_comparisons) + ylim(c(0.65, 1.05)) +
  ggtitle('sample size = 50') +
  scale_x_discrete(labels = c('1/5', '1/4', '1/3', '1/2', '1',
                              '2', '3', '4', '5'))
# ggsave(filename = 'Figures/Supp_SimulationsIndelLikelihoods_size50.pdf', plot = g50,
#        height = 5, width = 5)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_IndelWeighting_sim50.pdf', plot = g50,
       height = 7, width = 4)

# Change sample size
sampleSize = 100
res.full_100 <- data.frame()

set.seed(123)
for (i in 1:length(cluster_distribution.list)) {
  
  print(paste0('Indel alteration = ', indel_alterations[i]))
  
  res.i_100 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = sampleSize, 
    cluster_distributions_index = i, 
    cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_100 <- rbind(res.full_100, res.i_100$results)
  
}

# Plot results for HRD Pheno
res_100.sbs3 <- res.full_100[res.full_100$Pheno == 'SBS3', ]

my_comparisons <- list(c('0.2', '1'), c('1', '5'))

g100 <- ggboxplot(res_100.sbs3, x = 'indel_alteration', y = 'AUC', add = 'jitter') +
  stat_compare_means(comparisons = my_comparisons) + ylim(c(0.65, 1.05)) +
  ggtitle('sample size = 100') +
  scale_x_discrete(labels = c('1/5', '1/4', '1/3', '1/2', '1',
                              '2', '3', '4', '5'))
# ggsave(filename = 'Figures/Supp_SimulationsIndelLikelihoods_size100.pdf', plot = g100,
#        height = 5, width = 5)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_IndelWeighting_sim100.pdf', plot = g100,
       height = 7, width = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ICGC_PhenotypeClusteringAndSimulationAnalysis/ICGC_simulations_clusterReassignment.R`:

```````R
#####
## Simulation analysis to determine overall subsampling reclassification in ICGC-BRCA cohort
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(pROC)
library(ggpubr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# Complete SBS/indel counts, prior probabilities, and cluster likelihood spectra
load('Data/BRCA_UKEU_mt_tally.Rdata')
mut_complete <- as.data.frame(cbind(mt_tally.brca_wgs$SBS_96, mt_tally.brca_wgs$ID_83))

load('Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.Rdata')

load('Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.Rdata')
mut_complete <- mut_complete[rownames(ann), ]
pheno_assigned <- ann$Phenotype


# Write a likelihood function that aligns a dataset with the designated mean distribution
likelihood_calc <- function(input_data, cluster_distributions, cluster_assign) {
  
  # This function:
  #   Takes a dataset of mutations as its input (rows = samples, cols = mutation types)
  #   Limits the data to to designated mutation types
  #   Applies a likelihood approach to calculate posterior probabilities of cluster assignment
  #   NB, For now, we ensure that our cluster distributions have the correct mutation types
  #   Output: posterior distributions, note that for WES data a comparison column is not possible
  
  # Calculate likelihood: P(mutation spectrum | cluster)
  log_likelihoods <- matrix(NA, nrow = nrow(input_data), ncol = nrow(cluster_distributions))
  rownames(log_likelihoods) <- rownames(input_data); colnames(log_likelihoods) <- rownames(cluster_distributions)
  for (i in 1:nrow(input_data)) {
    # print(i)
    log_likelihoods[i, ] <- apply(cluster_distributions, 1, 
                                  function(x) sum(log10(x) * input_data[i, ]))
  }
  
  # Set prior: P(cluster)
  marginal.probs <- table(cluster_assign)/length(cluster_assign)
  marginal.probs <- marginal.probs[colnames(log_likelihoods)]
  
  # Calculate log_posteriors and return final posteriors
  log_posteriors <- data.frame(t(apply(log_likelihoods, 1, function(x) 
    x + log10(marginal.probs)
  )))
  
  # Generate final posteriors
  final_probs <- log_posteriors
  Sys.time(); for (i in 1:ncol(final_probs)) {
    final_probs[,i] <- 10^(apply(log_posteriors, 1, function(x) x[i] - max(x[1:ncol(final_probs)])))
  }
  final_probs <- data.frame(t(apply(final_probs, 1, function(x) x/sum(x)))); Sys.time()
  
  return(final_probs)
  
}



# Apply likelihood function to simulated data
simulate_likelihood_calc <- function(input_data, cluster_distributions, 
                                     cluster_assign, mutation_types = 'ALL', sample_size) {
  
  # input_data: data frame of samples (rows) and 96 trinucloetide contexts (cols)
  # sample_size: number of mutations to be sampled from patient with replacement
  
  # Limit input data to designated mutation types
  index.sbs <- sapply(names(input_data), function(x) grepl(pattern = '>', x))
  index.id <- sapply(names(input_data), function(x) grepl(pattern = ':', x)) 
  
  if (mutation_types == 'SBS') input_data <- input_data[, index.sbs]
  if (mutation_types == 'ID') input_data <- input_data[, index.id]
  
  # Develop simulation matrix (same dimensions as input_data, but with sampled data)
  sims <- matrix(0, nrow = nrow(input_data), ncol = ncol(input_data))
  colnames(sims) <- colnames(input_data); rownames(sims) <- rownames(input_data)
  
  for (i in 1:nrow(input_data)) {
    s <- table(sample(colnames(input_data), size = sample_size, replace = TRUE,
                      prob = input_data[i, ]/sum(input_data[i, ])))
    sims[i,names(s)] <- s
  }
  
  # Apply likelihood function to the simulated data
  posteriors <- likelihood_calc(sims, cluster_distributions, 
                                cluster_assign)
  
  # As this is simulated data, we can add the comparisons
  posteriors$PhenoTrue <- cluster_assign
  
  # Save and return AUC values
  auc.df <- data.frame(Phenotype = unique(posteriors$PhenoTrue),
                       AUC = NA)
  
  for (pheno in posteriors$PhenoTrue) {
    roc.full <- roc(posteriors$PhenoTrue == pheno, posteriors[,pheno], quiet = TRUE)
    auc.val <- roc.full$auc
    auc.df$AUC[auc.df$Phenotype == pheno] <- auc.val
  }
  
  return(list(auc.df = auc.df,
              posteriors = posteriors))
  
}


# Write function to run n simulations and return a data frame of AUC values
run_likelihood_sims <- function(input_data, sample_size, 
                                cluster_distributions, cluster_assign, 
                                mutation_types = 'ALL', n_simulations) {
  
  # All of the same inputs as simulate_likelihood_calc() + n_simulations
  #   Again, this function is about comparing simulated data with true results
  
  # Initialise matrix
  results.mat <- matrix(NA, nrow = 0, ncol = nrow(cluster_distributions)); colnames(results.mat) <- rownames(cluster_distributions)
  
  # Save final posterior distributions in list
  posteriors_list <- list()
  
  # Run n_simulations of simulate_likelihood_calc():
  for (i in 1:n_simulations) {
    print(paste0('Running simulation ', i, ' of ', n_simulations, '...'), quote = FALSE)
    post.i <- simulate_likelihood_calc(input_data = input_data, sample_size = sample_size,
                                       cluster_distributions = cluster_distributions, cluster_assign = cluster_assign,
                                       mutation_types = mutation_types)
    posteriors_list[[i]] <- post.i$posteriors[,1:(ncol(post.i$posteriors)-1)]
    
    auc.i <- post.i$auc.df$AUC
    results.mat <- rbind(results.mat, auc.i)
  }
  
  results <- data.frame(Pheno = rep(colnames(results.mat), each = nrow(results.mat)),
                        AUC = as.vector(results.mat))
  
  # Return collated results and final posterior distribution
  # return(list(results = results, posterior_n = post.i$posteriors))
  return(list(results = results, posteriors_list = posteriors_list))
  
}


# Run likelihood simulations and save/plot output
set.seed(123)
sig.type = 'ALL' # one of c('ALL','SBS','ID)
n_simulations = 100

for (sampleSize in c(25, 50, 100)) {
  
  print(paste0('Running simulations: Sample Size = ', sampleSize, '...'), quote = FALSE)
  
  # Run simulation
  res.simSampleSize <- run_likelihood_sims(input_data = mut_complete, sample_size = sampleSize,
                                           cluster_distributions = mut.dists_mean, cluster_assign = pheno_assigned,
                                           mutation_types = sig.type, n_simulations = n_simulations)
  res.simSampleSize$results$sample_size = sampleSize
  
  # Separate results and plot AUCs
  res_results <- res.simSampleSize$results
  write.table(res_results, file = paste0('Results/ICGC_simulations_AUCs_sims', sampleSize, '.txt'),
              quote = FALSE, sep = '\t', row.names = FALSE)
  
  res_results$Group <- sapply(as.character(res_results$Pheno),
                              function(x) strsplit(x,split='_')[[1]][1])
  res_results$Group[grepl(pattern = 'ID', res_results$Group)] <- 'ID_enriched'
  res_results$Pheno <- factor(res_results$Pheno,
                              levels = names(table(res_results$Pheno))[length(unique(res_results$Pheno)):1])
  
  g_auc <- ggboxplot(data = res_results, x = 'Pheno', y = 'AUC',
                     color = 'Group', orientation = 'horizontal')
  ggsave(filename = paste0('Figures/Supp_ICGCsimulations_PhenoReassign_sim', sampleSize, '_AUCs.pdf'), plot = g_auc)
         
  # Plot overall posterior distributions
  res_totalPosterior <- do.call(rbind, res.simSampleSize$posteriors_list)
  res_totalPosterior$Pheno_Assigned <- apply(res_totalPosterior, 1, function(x) names(x)[x==max(x)])
  res_totalPosterior$Pheno_True <- rep(ann$Phenotype, n_simulations)
  
  res_totalPosterior_summary <- as.matrix(table(res_totalPosterior$Pheno_True, res_totalPosterior$Pheno_Assigned))
  res_totalPosterior_summary <- apply(res_totalPosterior_summary, 1, function(x) x/sum(x))
  
  cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)
  pheatmap(res_totalPosterior_summary, color = cols_scale,
           cluster_rows = FALSE, cluster_cols = FALSE,
           filename = paste0('Figures/Supp_ICGCsimulations_PhenoReassign_sim', sampleSize, '_posteriorHeatmap.pdf'))
  
}


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ICGC_PhenotypeClusteringAndSimulationAnalysis/ICGC_PhenotypeClusterDevelopment_normalised.R`:

```````R
##### 
## Creation of signature phenotypes in 614 ICGC-BRCA samples
#####

setwd('~/Projects/HRD_MutationalSignature/Results/')

# Load libraries
library(mclust)
library(pheatmap)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# Load ICGC deconstructSigs data
load('ICGC_BRCA_deconstructSigs_Cutoff0.01_SBSandIDnormalised.Rdata')

# Run mixture modelling using mclust
sigs.BIC <- mclustBIC(sigs_complete, G = 10:25,
                      modelNames = c('EEI','VII','EII'))
summary(sigs.BIC) # Top clusters: VEI, 19, 20, 16

pdf(file = '../Figures/Supp_ICGCSignaturesMixtureModelling.pdf',
    width = 5, height = 5)
plot(sigs.BIC) # save this plot
dev.off()

# Apply optimal clustering in sigs.BIC
mod_sigs.BIC <- Mclust(sigs_complete, x = sigs.BIC)
table(mod_sigs.BIC$classification) # classification

## Form annotation document for visualisation

# CHORD
chord <- read_excel('~/Data/ICGC/CHORD_output.xlsx', sheet = 'CHORD')
chord <- chord[!grepl(pattern = 'HMF', chord$group), ]
chord <- chord[,c('sample','response','hr_status','hrd_type')]
names(chord)[3:4] <- c('CHORD', 'CHORD_type')

# HRDetect
hrdetect.samples <- read_excel('~/Data/ICGC/HRDetect_sampleInfo.xlsx', sheet = 'Data', skip = 2)
hrdetect.pred <- read_excel('~/Data/ICGC/HRDetect_predictions.xlsx', sheet = 'b.Predictor', skip = 2)
hrdetect.pred <- hrdetect.pred[match(hrdetect.samples$Sample, hrdetect.pred$sample), ]

hrdetect <- merge(x = hrdetect.samples[,c(1,3:8)],
                  y = hrdetect.pred[,c('sample','predictorProb')],
                  by.x = 'Sample', by.y = 'sample')
hrdetect$Gene[hrdetect$isBrcaMonoallelic] <- NA
hrdetect <- hrdetect[,c('Sample','ER status','Gene','predictorProb')]
names(hrdetect) <- c('sample','ER_status','BRCA_defect', 'HRDetect')
hrdetect$sample <- sapply(hrdetect$sample,
                          function(x) substr(x, 1, nchar(x)-1)) # remove last letter

# HRDetect validation: sensitiity = 76/77 = 98.7%
table(hrdetect$BRCA_defect, hrdetect$HRDetect > .7,
      dnn = c('BRCA_defect', 'score > .7'), useNA = 'always')

# Load ICGC sample data (for sampleID matching)
samples.eu <- read.table('~/Data/ICGC/sample.BRCA-EU.tsv.gz', h=T, sep='\t')
samples.uk <- read.table('~/Data/ICGC/sample.BRCA-UK.tsv.gz', h=T, sep='\t')
samples.icgc <- rbind(samples.eu, samples.uk)
samples.icgc <- samples.icgc[,c('project_code','submitted_sample_id','icgc_donor_id')]
samples.icgc <- samples.icgc[!duplicated(samples.icgc$icgc_donor_id), ]

samples.icgc$final_letter <- sapply(samples.icgc$submitted_sample_id,
                                    function(x) substr(x,nchar(x),nchar(x)))
samples.icgc <- samples.icgc[samples.icgc$final_letter %in% c('a','b'),]
samples.icgc$sample_id <- sapply(samples.icgc$submitted_sample_id,
                                 function(x) substr(x,1,nchar(x)-1))
samples.icgc <- samples.icgc[!duplicated(samples.icgc$sample_id), ]

# Combine sample data with HRDetect and CHORD
ann.icgc <- merge(x = samples.icgc, y = hrdetect,
                  by.x = 'sample_id', by.y = 'sample', all.x = TRUE)
ann.icgc <- merge(x = ann.icgc, y = chord,
                  by.x = 'sample_id', by.y = 'sample', all.x = TRUE)
rownames(ann.icgc) <- ann.icgc$icgc_donor_id
ann.icgc <- ann.icgc[,c('ER_status', 'BRCA_defect',
                        'HRDetect', 'CHORD', 'CHORD_type')]

# Create annotation with finite mixture model clusters
ann <- data.frame(BIC_clust = factor(mod_sigs.BIC$classification,
                                     levels = 1:length(unique(mod_sigs.BIC$classification))))
ann <- merge(x = ann, y = ann.icgc,
             by = 0, all.x = TRUE)
rownames(ann) <- ann$Row.names; ann <- ann[,-1]
ann <- ann[order(ann$BIC_clust), ]

# Order samples by classification
sigs_order <- as.data.frame(t(sigs_complete[rownames(ann), ]))

# Set colours
cols <- colorRampPalette(brewer.pal(9,'Set1'))
cols_BIC <- cols(length(unique(ann$BIC_clust)))
names(cols_BIC) <- unique(ann$BIC_clust)

ann_colors = list(
  BIC_clust = cols_BIC,
  ER_status = c(positive = 'darkgreen', negative = 'yellow'),
  BRCA_defect = c(BRCA1 = 'blue', BRCA2 = 'red'),
  CHORD = c(cannot_be_determined = 'grey', HR_deficient = 'black', HR_proficient = 'white'),
  CHORD_type = c(cannot_be_determined = 'grey', BRCA1_type = 'blue', BRCA2_type = 'red', none = 'white')
)

# Use white -> navy scale
cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)

pheatmap(sigs_order, cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col = ann, annotation_colors = ann_colors,
         color = cols_scale)

# Naming clusters as signature phenotypes (by sight based on above heatmap)
pheno <- c('SBS5_1','HRD_ID8','HRD_ID6high','HRD_APOBEC','SBS5_SBS18',
           'SBS5_2','SBS5_3','SBS5_4','APOBEC_ID9','SBS5_5',
           'HRD_SBS8','APOBEC_SBS2','APOBEC_SBS13','SBS5_ID5','SBS5_6',
           'HRD_ID9','HRD_ID4','HRD_ID6mid','MMRD','SBS5_SBS39',
           'clustSmall1','clustSmall2')
# pheno <- c('SBS5_1', 'HRD_ID8high', 'HRD_ID9high', 'HRD_APOBEC', 'HRD_ID6high',
#            'HRD_ID8mid', 'SBS5_2', 'SBS5_SBS18', 'SBS5_3', 'SBS5_4',
#            'APOBEC_SBS2', 'SBS5_5', 'APOBEC_ID9', 'SBS5_6', 'SBS5_7',
#            'APOBEC_SBS13', 'SBS5_8', 'ID4', 'HRD_ID9mid', 'SBS5_ID1',
#            'HRD_SBS8', 'SBS5_SBS39', 'HRD_ID6mid', 'MMRD')

ann$Phenotype <- factor(pheno[as.numeric(ann$BIC_clust)],
                        levels = sort(pheno))

# Remove uninformative clusters:
ann <- ann[!(ann$Phenotype %in% c('clustSmall1','clustSmall2')), ]

pheno <- pheno[!(pheno %in% c('clustSmall1','clustSmall2'))]
ann$Phenotype <- factor(ann$Phenotype, levels = sort(pheno))

save(ann, file = '~/Projects/HRD_MutationalSignature/Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.Rdata')

ann <- ann[order(ann$Phenotype), ]

sigs_order <- sigs_order[, rownames(ann)]

# Redo heatmap with phenotype labels
cols_Pheno <- cols(length(unique(ann$Phenotype)))
names(cols_Pheno) <- unique(ann$Phenotype)

ann_colors = list(
  Phenotype = cols_Pheno,
  ER_status = c(positive = 'darkgreen', negative = 'yellow'),
  BRCA_defect = c(BRCA1 = 'blue', BRCA2 = 'red'),
  CHORD = c(cannot_be_determined = 'grey', HR_deficient = 'black', HR_proficient = 'white'),
  CHORD_type = c(cannot_be_determined = 'grey', BRCA1_type = 'blue', BRCA2_type = 'red', none = 'white')
)

pheatmap(sigs_order, cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col = ann[,c('Phenotype','BRCA_defect', 'ER_status',
                                 'HRDetect', 'CHORD', 'CHORD_type')], 
         annotation_colors = ann_colors,
         color = cols_scale, fontsize = 8, fontsize_row = 10
         ,filename = '../Figures/Supp_ICGCSignaturesHeatmap.pdf'
         )

## Investigate BRCA-defect distribution across HRD classifications
##   using the 'ann' data frame

ann.hrd <- ann
ann.hrd$HRD_cluster <- sapply(ann.hrd$Phenotype,
                              function(x) ifelse(grepl(pattern = 'HRD', x), 
                                                 'HRD', 'HR-proficient'))

# HRD clusters: sensitivity = 98.7%, specificity = 43.8%
table(ann.hrd$BRCA_defect, ann.hrd$HRD_cluster, 
      dnn = c('BRCA_defect', 'HRD cluster'), useNA = 'always')

ann.hrd$BRCA_defective <- ifelse(is.na(ann.hrd$BRCA_defect), 'BRCA+', 'BRCA_defective')
ann.hrd$HRD <- ifelse(grepl(pattern = 'HRD', ann.hrd$Phenotype), 'HRD', 'HR-proficient')

# HRDetect > 0.7: sensitivity = 98.7%, specificity = 61.7%
table(ann.hrd$BRCA_defect, ann.hrd$HRDetect > 0.7,
      dnn = c('BRCA_defect', 'HRDetect > 0.7'), useNA = 'always')

# CHORD: sensitivity = 85.3%, specificity = 66.7%
table(ann.hrd$BRCA_defect, ann.hrd$CHORD,
      dnn = c('BRCA_defect', 'CHORD'), useNA = 'always')

ann.hrd_summary <- ann.hrd %>%
  group_by(HRD, BRCA_defective) %>%
  summarise(n = n())
ann.hrd_summary$HRD <- factor(ann.hrd_summary$HRD, levels = c('HR-proficient','HRD'))
ann.hrd_summary$BRCA_defective <- factor(ann.hrd_summary$BRCA_defective, levels = c('BRCA_defective','BRCA+'))

g_annHRDsummary <- ggplot(ann.hrd_summary, aes(x = BRCA_defective, y = n, fill = HRD)) +
  geom_bar(stat = 'identity', position = 'fill') + theme_minimal() +
  scale_fill_brewer(palette = 'Paired') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())
ggsave(filename = '../Figures/Supp_ICGCHRDclassify.pdf',
       plot = g_annHRDsummary, width = 5, height = 6)


## BRCA-type specific HRD classification

hrd_brca1type <- c('HRD_APOBEC', 'HRD_ID6mid', 'HRD_ID8', 'HRD_SBS8')
hrd_brca2type <- c('HRD_ID6high')
hrd_undefined <- c('HRD_ID4','HRD_ID9')

ann.hrd$BRCAtype_HRD <- 'HR-proficient'
ann.hrd$BRCAtype_HRD[ann.hrd$Phenotype %in% hrd_brca1type] <- 'BRCA1-type HRD'
ann.hrd$BRCAtype_HRD[ann.hrd$Phenotype %in% hrd_brca2type] <- 'BRCA2-type HRD'
ann.hrd$BRCAtype_HRD[ann.hrd$Phenotype %in% hrd_undefined] <- 'HRD unassigned'

ann.hrd$BRCA_defect_label <- ann.hrd$BRCA_defect
ann.hrd$BRCA_defect_label[is.na(ann.hrd$BRCA_defect_label)] <- 'BRCA+'

ann.brca_summary <- ann.hrd %>%
  group_by(BRCA_defect_label, BRCAtype_HRD) %>%
  summarise(n = n())
ann.brca_summary$BRCA_defect_label <- factor(ann.brca_summary$BRCA_defect_label,
                                             levels = c('BRCA1','BRCA2','BRCA+'))
ann.brca_summary$BRCAtype_HRD <- factor(ann.brca_summary$BRCAtype_HRD,
                                        levels = c('HR-proficient', 'HRD unassigned', 'BRCA2-type HRD', 'BRCA1-type HRD'))

g_annBRCAsummary <- ggplot(ann.brca_summary, aes(x = BRCA_defect_label, y = n, fill = BRCAtype_HRD)) +
  geom_bar(stat = 'identity', position = 'fill') + theme_minimal() +
  scale_fill_manual(values = c('grey90', 'grey50','red','blue')) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank()) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))
ggsave(filename = '../Figures/Supp_ICGCBRCAclassify.pdf',
       plot = g_annBRCAsummary, width = 6, height = 7.5)

# CHORD:
#   BRCA1: sensitivity = 73.3%, specificity = 55.0%
#   BRCA2: sensitivity = 93.3%, specificity = 77.8%
table(ann.hrd$BRCA_defect, ann.hrd$CHORD_type,
      dnn = c('BRCA_defect', 'CHORD'), useNA = 'always')


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/ExomeClassifier/ICGC_PhenotypeClusteringAndSimulationAnalysis/ICGC_indelCounts_ExomeVsGenome.R`:

```````R
setwd('~/Data/ICGC/')
#####
## Processing of ICGC WGS breast cancer samples to infer SBS and indel mutation type contributions
####

# Load libraries
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)

## Load UK and EU BRCA data (US is exome, FR contains no indels)
# Extract WGS strategies, remove duplicated mutations

brca_uk <- read.table('~/Data/ICGC/simple_somatic_mutation.open.BRCA-UK.tsv.gz',
                      sep = '\t', header = TRUE)
brca_uk_wgs <- brca_uk[brca_uk$sequencing_strategy == 'WGS', ]
brca_uk_wgs <- brca_uk_wgs[!duplicated(brca_uk_wgs$icgc_mutation_id), ]

brca_eu <- read.table('~/Data/ICGC/simple_somatic_mutation.open.BRCA-EU.tsv.gz',
                      sep = '\t', header = TRUE)
brca_eu_wgs <- brca_eu[brca_eu$sequencing_strategy == 'WGS', ]
brca_eu_wgs <- brca_eu_wgs[!duplicated(brca_eu_wgs$icgc_mutation_id), ]

# Collate ICGC data and extract mutations in intron or intergenic regions

brca_wgs <- rbind(brca_uk_wgs, brca_eu_wgs)
brca_wgs_input <- data.frame(
  Tumor_Sample_Barcode = brca_wgs$icgc_donor_id,
  Hugo_Symbol = NA,
  Chromosome = brca_wgs$chromosome,
  Start_position = brca_wgs$chromosome_start,
  End_position = brca_wgs$chromosome_end,
  Variant_Classification = brca_wgs$consequence_type,
  Variant_Type = sapply(brca_wgs$mutation_type,
                        function(x) ifelse(x == 'single base substitution', 'SNP',
                                           ifelse(x == 'insertion of <=200bp', 'INS',
                                                  ifelse(x == 'deletion of <=200bp', 'DEL',NA)))),
  Reference_Allele = brca_wgs$reference_genome_allele,
  Tumor_Seq_Allele2 = brca_wgs$mutated_to_allele
)

# brca_exome_input <- brca_wgs_input[!(brca_wgs_input$Variant_Classification %in%
#                                        c('intergenic_region', 'intron_variant')), ]
brca_exome_input <- brca_wgs_input[brca_wgs_input$Variant_Classification != 'intergenic_region',]

brca_wgs_maf <- read.maf(maf = brca_wgs_input,
                         vc_nonSyn = names(table(brca_wgs_input$Variant_Classification)))

brca_exome_maf <- read.maf(maf = brca_exome_input,
                           vc_nonSyn = names(table(brca_exome_input$Variant_Classification)))

# Calculate ID-83 counts for each MAF file
mt_tally.brca_wgs <- sig_tally(
  object = brca_wgs_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ID',
  useSyn = TRUE
)

mt_tally.brca_exome <- sig_tally(
  object = brca_exome_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ID',
  useSyn = TRUE
)

# Count total mutations of each ID-83 type
icgc_indelCounts <- data.frame(
  wgs = apply(mt_tally.brca_wgs$all_matrices$ID_83, 2, sum),
  exome = apply(mt_tally.brca_exome$all_matrices$ID_83, 2, sum)
)

write.table(icgc_indelCounts, file = 'ICGC_BRCA_indelCounts.txt')

icgc_indelCounts$ex_gen <- icgc_indelCounts$exome/icgc_indelCounts$wgs
hist(icgc_indelCounts$ex_gen)

library(ggplot2)
ggplot(icgc_indelCounts, aes(x = log10(wgs), y = log10(exome))) +
  geom_point() + geom_smooth(method = 'lm')



```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/README.md`:

```````md
# Multiscale HRD classification

### Author: Daniel H Jacobson, UCL

This repository contains the workflow for classification of homologous recombination deficient breast cancers via mutational signature classification and transcriptomic signature approaches. Scripts are separated into these two groups, and are presented here in the order in which they must be executed.

![alt text](TCGA_HRDclassificationHeatmap.jpg)

# Scripts

## Exome Classifier

- **ICGC_BRCA_UKandEU_MatrixGeneration.R:** Collates mutation data from the BRCA-EU and BRCA-UK projects and uses the sigminer R package to collate the mutational spectra of each sample according to the SBS96 and ID83 signatures.
- **ICGC_deconstructSigs_genome2exome.R:** Calculates contributions of breast cancer-associated SBS and ID mutational signatures to each of the ICGC samples, establishing the signature profiles of each sample, alongside correction for genome-to-exome normalisation.
- **ICGC_PhenotypeClusterDevelopment_normalised.R:** Applies finite mixture modelling to cluster the ICGC samples depending on their signature profiles, and enables name assignment to each cluster generated based on their most prevalent signature contributions.
- **LikelihoodCluster_Generation.R:** Collates the signature profiles for all samples within a cluster and generates a mean mutational spectrum representative of them, therefore generating the probability distributions for each cluster.
- **TCGA_HRDclassification.R:** Using the prior probabilities and likelihoods generated using ICGC, we collate the mutational profiles of 986 exome sequenced breast cancers from TCGA and calculate the posterior probabilities of assignment to each of the ICGC-generated clusters. The sum of probabilities of assignment to the HRD-associated clusters equals to the probability of HRD.
- **TCGA_HRDhallmarks.R:** Comparison of HRD and HR-proficiency assignments in TCGA to HRD-associated features (Myriad HRD score, CX3 copy number signature contribution, POLQ expression, proliferation capacity)

## TranscriptionalSignature

- **TCGA_BRCA.RNAseq_prep.R:** Pre-processing of TCGA-BRCA expression data including removal of lowly expressed genes, cancer cell expression deconvolution using BayesPrism, and separation into training (~2/3) and testing (~1/3) cohorts. This includes both HRD/HR-proficiency assignment according to the exome classifier, as well as BRCA1/BRCA2/HRD_BRCA+/HR-proficiency classifications which are used for signature development.
- **MultinomialElasticNet_alpha0.25_1to100.R:** Performs 100 iterations of 10-fold cross-validated multinomial elastic net regression. On each iteration, the gene parameter coefficients are saved. To run all 1000 iterations, an addition nine scripts were run in parallel with the seed set to the iteration value. Additionally, these analyses were repeated with alpha = 0.5 to generate an alternative signature.
- **CentroidModelFormation.R:** Collates the coefficients from the 1000 iterations of elastic net regression and extracts the 228 genes which appear as non-zero in all of them. The median expression of each gene is calculated across te HRD/HR-proficiency and HRD/BRCA-defect groups to generate the templates.
- **TCGA_testScoring.R:** Correlate the TCGA-BRCA testing cohort against each of the templates, saving the Pearson's correlation coefficient which represents the associated 'score'. The 'HRD score' is calculated by subtracting the correlation with the HR-proficiency template against the correlation with the HRD template.
- **GSEA_pathfindR.R:** Runs gene set enrichment analysis using the pathfindR tool. To enable this, for each gene an ANOVA is run of correlation against the HRD/BRCA-defect group, with the significance saved and adjusted. 
- **CCLE_jointComparisons:** Analysis of associations between HRD scores and PARP inhibitor sensitivity in breast cancer cell lines obtained from the Cancer Cell Line Encyclopaedia, and comparison against alternative HRD signatures.
- **ISPY2_HRDscoring.R:** Analysis of HRD scores across breast cancer patients treated with olaparib and durvalumab as part of the I-SPY2 trial, and comparison with the PARPi7 score
- **Chung2017_analysis.R:** Separates the bulk and single cell expression profiles from the Chung 2017 single cell breast cancer atlas. The tumour cells are extracted from the single cell data. HRD scores are calculated for each sample and tumour cell, and the sample-wide HRD scores are compared with the mean scores generated in the individual tumour cells.
- **Qian2020_preprocessing.R:** Preprocessing of the Qian 2020 breast cancer cohort, including removal of cells with high cell stress and unreasonable gene counts, and normalisation of expression scores.
- **Qian2020_analysis.R:** Analyses of distributions of HRD scores across the Qian 2020 cohort, and file preparation for CellphoneDB analysis
- **Qian2020_HRDprofiling.R** Analyses of distributions of HRD scores across cancer cells and the tumour microenvironment across the Qian 2020 cohort and UMAP plotting.


# Copyright

This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/TranscriptionalSignatureDevelopment/RNAseqTransformationScript.R`:

```````R
setwd('/home/zcqsdhj/TranscriptionalSignatures/Revisions/Data')

load('TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')
z2 <- Z.tumor_training[,apply(Z.tumor_training, 2, sum) > 0]
z3 <- log2(z2+1)
z4 <- apply(z3, 2, scale)
rownames(z4) <- rownames(z3)

Z.tumor_training <- z4
save(Z.tumor_training, file = 'TCGA_BRCA.ExprDeconvolution_270723_p0.79_training.Rdata')

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/TranscriptionalSignatureDevelopment/MultinomialWEN_alpha0.25_1to250.R`:

```````R
#####
## Multinomial Elastic Net Regression to classify BRCA1-/-, BRCA2-/-, HRD, HR-proficient TCGA-BRCA samples
#####

setwd('~/TranscriptionalSignatures/Revisions')

# Load libraries
library(glmnet)

# Load in training data:
#   TCGA-BRCA, expression deconvolution by BayesPrism, should be followed by low-gene removal and log2-normalisation
#   Genes remain if they have >1 count in >70% samples
#   HRD/BRCA status
load('Data/TCGA_BRCA.ExprDeconvolution_270723_p0.79_training.Rdata')
rownames(Z.tumor_training) <- substr(rownames(Z.tumor_training), 1, 12)
#index.keep <- apply(Z.tumor_training, 2, function(x) sum(x > 1)/length(x) > 0.7)
#Z.tumor_training <- Z.tumor_training[,index.keep]
#Z.tumor_training <- log2(Z.tumor_training + 1)

# Organise output variable
load('Data/TCGA_HRDclassification_BRCAannotation_HRD0.79.Rdata')
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

# Match output to training data
samples.intersect <- intersect(rownames(Z.tumor_training), rownames(ann_tcga))

ann_tcga <- ann_tcga[samples.intersect, ]
ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

Z.tumor_training <- Z.tumor_training[samples.intersect, ]


# Separate out expression data from BRCA_defect status
input.x <- Z.tumor_training
input.y <- ann_tcga$group

# Calculate weightings
weight_function <- function(index, vector) {
        return(1/sum(vector == vector[index]))
}
weights.input <- sapply(1:length(input.y), function(x) weight_function(x,input.y))

# Conduct 500 iterations of 10-fold cross validation:
#   (500 iteration ~ 2 days)
#   For each iteration, calculate coefficients using cv$lambda.min
#     (lambda coefficient which gives smallest error)
#   set seed inside loop

iterations <- 1:250
coefs <- list()

for (i in iterations) {
  
  # Set seed with each iteration
  set.seed(123*i)
  
  print(paste0('Iteration ', i, ' of ', iterations[length(iterations)], '...', Sys.time()))
  
  # Conduct grouped multinomial 10-fold cross validation
  cv <- cv.glmnet(x = input.x, y = input.y, family = 'multinomial',
                  alpha = 0.25, type.multinomial = 'grouped', weights = weights.input)
  coefs[[i]] <- coef(cv$glmnet.fit, s = cv$lambda.min)
  
}

save(coefs, file = 'Results/CVcoefficients_MultiWEN_alpha0.25_p0.79_iter1to250.Rdata')

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/TranscriptionalSignatureDevelopment/ReduxSig_adjacencyMatrixFormation.R`:

```````R
expr <- read.csv('../Downloads/TcgaBrca_exprHRDsig.csv', row.names = 1)
expr <- as.data.frame(t(expr))

# Load important genes
imp.df <- read.csv('Data/imp_score_avg.csv')

genes.imp_hrd <- imp.df$gene_names[imp.df$HRD > 0.7]
genes.imp_hrprof <- imp.df$gene_names[imp.df$HR.proficient > 0.7]

genes.imp <- c(genes.imp_hrd, genes.imp_hrprof)
genes.imp_group <- c(rep('HRD', length(genes.imp_hrd)),
                     rep('HR-prof', length(genes.imp_hrprof)))

expr.redux <- expr[,genes.imp]

library(WGCNA)
library(tidyr)

adj <- as.data.frame(adjacency(expr.redux))
diag(adj) <- NA
adj$source <- rownames(adj)
adj$group <- genes.imp_group
adj[adj < 0.0002] <- NA

adj <- adj %>%
  pivot_longer(cols = -c(source, group), names_to = 'target', values_to = 'int')
adj <- adj[!is.na(adj$int), ]

adj$node1 <- apply(adj[,c('source','target')], 1, function(x) sort(x)[1])
adj$node2 <- apply(adj[,c('source','target')], 1, function(x) sort(x)[2])
adj$link <- apply(adj[,c('node1','node2')], 1, function(x) paste(x,collapse='_'))
adj <- adj[!duplicated(adj$link), ]

write.csv(adj[,c(1,3,4)], file = '~/Projects/HRD_TranscriptionalSignature/Results/adjacency_reduxSig.csv')

adj.df <- data.frame(
  Gene = genes.imp,
  Group = genes.imp_group
)
# adj.df <- adj.df[adj.df$name %in% c(adj$source,adj$target), ]
write.csv(adj.df, file = '~/Projects/HRD_TranscriptionalSignature/Results/adjInfo_reduxSig.csv')

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/TranscriptionalSignatureDevelopment/BayesPrism_TumorPurityComparison.R`:

```````R
library(TCGAbiolinks)
library(ggpubr)

load('~/Data/TCGA/TCGA_BRCA.BayesPrism.theta.Rdata')
df.theta <- data.frame(
  Sample.ID = sapply(rownames(theta), function(x) substr(x, 1, 16)),
  CancerCellFraction = theta[,'Cancer']
)

df.theta <- merge(x = df.theta, y = Tumor.purity[,c('Sample.ID','ESTIMATE')])

# Function for converting ESTIMATE column in Tumor.purity data.frame
comma.function <- function(x) {
  x.tmp = as.character(x)
  x.tmp = strsplit(x.tmp,split=',')[[1]]
  x.final = as.numeric(paste(x.tmp,collapse='.'))
  return(x.final)
}

df.theta$TumorPurity <- sapply(df.theta$ESTIMATE, comma.function)

# Compare estimates
g_tumorPurity <- ggplot(df.theta, aes(x = TumorPurity, y = CancerCellFraction)) +
  geom_point(alpha = 0.3) + geom_smooth(method = 'lm') + stat_cor() +
  theme_minimal()
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_BayesPrismEstimates.pdf',
       plot = g_tumorPurity)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/TranscriptionalSignatureDevelopment/CollateAlternativeSignatures.R`:

```````R
#####
## Script to generate centroid templates for alternative transcriptional signatures
#####

# Libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(readxl)

# Load TCGA transcriptional data
#   Templates will be formed from FPKM-normalised training data
#   Therefore, we must also load the training data to obtain sample IDs

load('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')
barcodes.training <- rownames(Z.tumor_training)

setwd('~/Data/TCGA')
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  barcode = barcodes.training
)
# GDCdownload(query)
expr.train <- GDCprepare(query = query)
# expr.train <- expr.train[,expr.train$sample_type == 'Primary Tumor']
expr.train <- expr.train[!duplicated(rowData(expr.train)$gene_name) &
                         !is.na(rowData(expr.train)$gene_name), ]
genes.tcga <- rowData(expr.train)$gene_name

# Get signatures
signatures_alternative <- list()

parpi7 <- c('BRCA1', 'MRE11', 'NBN', 'TDG', 'XPA', 'CHEK2', 'MAPKAPK2')
cin70 <- c('TPX2','PRC1','FOXM1','CDK1','TGIF2','MCM2','H2AZ1','TOP2A','PCNA','UBE2C',
           'MELK','TRIP13','CEP250','MCM7','RNASEH2A','RAD51AP1','KIF20A','CDC45','MAD2L1','ESPL1',
           'CCNB2','FEN1','TTK','CCT5','RFC4','ATAD2','CKAP5','NUP205','CDC20','CKS2',
           'RRM2','ELAVL1','CCNB1','RRM1','AURKB','MSH6','EZH2','CTPS1','DKC1','OIP5',
           'CDCA8','PTTG1','CEP55','H2AX','CMAS','NCAPH','MCM10','LSM4','NCAPG2','ASF1B',
           'ZWINT','PBK','ZWILCH','CDCA3','ECT2','CDC6','UNG','MTCH2','RAD21','ACTL6A',
           'PDCD2L','SRSF2','HDGF','NXT1','NEK2','DHCR7','AURKA','NDUFAB1','NEMP1','KIF4A')

signatures_alternative[['PARPi7']] <- parpi7
signatures_alternative[['CIN70']] <- cin70

# Severson
sev_init <- read_excel('../../../Downloads/BRCA1nessSignature_Severson2017.xlsx')
sev_init <- as.data.frame(sev_init)
sev <- sev_init[,1]
sev[which(!(sev %in% genes.tcga))] <- c('JAML','FAM241B','PIMREG',
                                        'HIF1A','PLAAT1','IRAG2')
signatures_alternative[['Severson']] <- sev

# Peng 2014
peng <- read_excel('../../../Downloads/HRDSignature_Peng2014.xlsx', skip=1)
peng <- as.data.frame(peng)
peng <- peng$`Gene Symbol`
peng[which(!(peng %in% genes.tcga))] <- c('MCMBP','FAM170B','DDIAS','SKA3','CEP128','TICRR',
                                          'TEDC2','METTL22','ATAD5','KIZ','ISM1','SMIM14',
                                          'SNHG32','DSCC1','DEFB1','DDX39A','HJURP','DLGAP5',
                                          'DNA2','RETREG1','H1-2','H2BC5','H2AC18','H2BC21',
                                          'HSP90AA2P','CREBRF','LOC554223','LOC649679','LOC729843','LOC91431',
                                          'VWA5A','ETFRF1','CENPU','MTARC1','BEX3','LRR1',
                                          'SRSF2','EPB41L4A-AS1','SLC35G1','TUBB4B','TUBB7P','WRAP53')
peng <- peng[peng %in% genes.tcga]
signatures_alternative[['Peng']] <- peng

save(signatures_alternative, file = '~/Projects/HRD_TranscriptionalSignature/AlternativeSignatures.Rdata')

# Collate centroids for alternative signatures
signature_alternative.centroid.list <- list()

# # Add Severson signature immediately
# rownames(sev_init) <- sev
# sev_init <- sev_init[,-1]
# colnames(sev_init) <- c('HRD','HR_proficient')
# signature_alternative.centroid.list[['Severson']] <- sev_init

# Prepare FPKM-normalised TCGA training data
expr.train_fpkm <- assay(expr.train, 'fpkm_unstrand')
rownames(expr.train_fpkm) <- rowData(expr.train)$gene_name
expr.train_fpkm <- log2(expr.train_fpkm + 1)
expr.train_fpkm <- apply(expr.train_fpkm, 1, scale)
rownames(expr.train_fpkm) <- sapply(expr.train$barcode, function(x) substr(x,1,12))
expr.train_fpkm <- as.data.frame(expr.train_fpkm)

# Add mutational signature/BRCA defect classification
load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('HRD','BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'
# ann_tcga$group[ann_tcga$HRD == 'HR-proficient'] <- 'HR-proficient'
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

samples.intersect <- intersect(rownames(expr.train_fpkm), rownames(ann_tcga))

ann_tcga <- ann_tcga[samples.intersect, ]
expr.train_fpkm <- expr.train_fpkm[samples.intersect, ]

# For the remaining three signatures:
#   Generate average centroids and add to signature list
for (signature in c('Severson','PARPi7','CIN70','Peng')) {
  sig.genes <- signatures_alternative[[signature]]
  
  centroid.signature <- data.frame(
    HRD = apply(expr.train_fpkm[ann_tcga$HRD == 'HRD',sig.genes], 2, mean),
    HR_proficient = apply(expr.train_fpkm[ann_tcga$HRD == 'HR-proficient',sig.genes],2,mean),
    BRCA1 = apply(expr.train_fpkm[ann_tcga$group == 'BRCA1',sig.genes],2,mean),
    BRCA2 = apply(expr.train_fpkm[ann_tcga$group == 'BRCA2',sig.genes],2,mean),
    HRD_BRCApos = apply(expr.train_fpkm[ann_tcga$group == 'HRD_BRCA+',sig.genes],2,mean),
    HR_BRCA_proficient = apply(expr.train_fpkm[ann_tcga$group == 'HR-proficient',sig.genes],2,mean)
  )
  
  if (signature == 'Severson') {
    centroid.signature$HRD = sev_init$`BRCA1ness template Pearson correlations`
    centroid.signature$HR_proficient = sev_init$`non-BRCAness template Pearson correlations`
  }
  
  signature_alternative.centroid.list[[signature]] <- centroid.signature
  
}

save(signature_alternative.centroid.list, file = '~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.Rdata')

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/TranscriptionalSignatureDevelopment/Qian2020_signatureDropOuts.R`:

```````R
setwd('~/Data/scRNASeq/Qian2020/')

load('exprData_Qian2020.Rdata')

meta.qian <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.qian <- meta.qian[match(colnames(expr.data_qian2020), meta.qian$Cell), ]
expr.cancer <- expr.data_qian2020[,meta.qian$CellType == 'Cancer']

# Collate genes lists
genes.g0 <- read.csv('~/Data/QuiescenceBiomarkers.csv')
genes.g0 <- intersect(genes.g0$Genes, rownames(expr.cancer))

load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
genes.hrd <- rownames(signature.centroid.list$ElasticNet_alpha0.25)

library(readxl)
genes.ddr_table <- readxl::read_excel('../../../../Downloads/Pearl_DDRgenes.xlsx', skip=1)
genes.ddr <- genes.ddr_table$`Gene ID`
genes.ddr <- genes.ddr[!is.na(genes.ddr)]
genes.ddr <- genes.ddr[!duplicated(genes.ddr)]
genes.ddr <- intersect(genes.ddr, rownames(expr.cancer))

# Summarise
summary(apply(expr.cancer[genes.g0,], 1, function(x) mean(x>0)))
hist(apply(expr.cancer[genes.g0,], 1, function(x) mean(x>0)), breaks = 50, 
     xlim = c(0,1), main = 'G0 genes expression')

summary(apply(expr.cancer[genes.ddr, ], 1, function(x) mean(x>0)))
hist(apply(expr.cancer[genes.ddr, ], 1, function(x) mean(x>0)), breaks = 50, 
     xlim = c(0,1), main = 'DDR genes expression')

summary(apply(expr.cancer[genes.hrd,], 1, function(x) mean(x>0)))
hist(apply(expr.cancer[genes.hrd,], 1, function(x) mean(x>0)), breaks = 50, 
     xlim = c(0,1), main = 'HRD genes expression')

summary(apply(expr.cancer[genes.hrd,], 2, function(x) sum(x>0)))

# 
# summary(apply(expr.cancer, 1, function(x) mean(x>0)))
# hist(apply(expr.cancer, 1, function(x) mean(x>0)), breaks = 50, 
#      xlim = c(0,1), main = 'All genes expression')


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/TranscriptionalSignatureDevelopment/CentroidModelFormation.R`:

```````R
# Collate centroid models for Runs

# Load deconvoluted training data for reference
load('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')
rownames(Z.tumor_training) <- sapply(rownames(Z.tumor_training), function(x) substr(x,1,12))
Z.tumor_training <- log2(Z.tumor_training + 1)
samples <- rownames(Z.tumor_training)
Z.tumor_training <- apply(Z.tumor_training, 2, scale)
rownames(Z.tumor_training) <- samples

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
# ann_tcga$HRD <- sapply(ann_tcga$HRD_prob, function(x) ifelse(x >= 0.32, 'HRD', 'HR-proficient'))
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'
# ann_tcga$group[ann_tcga$HRD == 'HR-proficient'] <- 'HR-proficient'
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

samples.intersect <- intersect(rownames(Z.tumor_training), rownames(ann_tcga))

ann_tcga_train <- ann_tcga[samples.intersect, ]
Z.tumor_training <- Z.tumor_training[samples.intersect, ]

# Z.tumor_training <- apply(Z.tumor_training, 2, scale)
# rownames(Z.tumor_training) <- rownames(ann_tcga_train)

# # Apply to non-deconvoluted samples
# library(TCGAbiolinks)
# library(SummarizedExperiment)
# 
# setwd('~/Data/TCGA')
# query <- GDCquery(
#   project = 'TCGA-BRCA',
#   data.category = 'Transcriptome Profiling',
#   data.type = 'Gene Expression Quantification',
#   workflow.type = 'STAR - Counts',
#   barcode = rownames(Z.tumor_training)
# )
# # GDCdownload(query)
# expr.train <- GDCprepare(query = query)
# expr.train <- expr.train[,expr.train$sample_type == 'Primary Tumor']
# expr.train <- expr.train[!duplicated(rowData(expr.train)$gene_name) &
#                          !is.na(rowData(expr.train)$gene_name), ]
# 
# library(SummarizedExperiment)
# expr.tumor_training <- assay(expr.train, 'fpkm_uq_unstrand')
# rownames(expr.tumor_training) <- rowData(expr.train)$gene_name
# colnames(expr.tumor_training) <- sapply(colnames(expr.tumor_training),
#                                        function(x) substr(x,1,12))
# expr.tumor_training <- log2(expr.tumor_training+1)
# 
# expr.tumor_training <- t(expr.tumor_training)
# expr.tumor_training <- expr.tumor_training[rownames(ann_tcga_train), ]

# Move to CV coefficients and generate centroids
setwd('~/Projects/HRD_TranscriptionalSignature/CV_Coefficients/Run6/')
signature.centroid.list <- list()

models <- list.files()
models <- sapply(models, function(x) paste0(strsplit(x,split='_')[[1]][2:3],collapse='_'))
models <- sapply(models, function(x) substr(x,6,nchar(x)))
models <- unique(models)

for (model in models) {
  
  print(model)
  
  files.coefs <- list.files(pattern = model)
  coefs_join <- list()
  for (file in files.coefs) {
    load(file)
    coefs_join <- c(coefs_join, coefs)
    rm(coefs)
  }
  
  coef.mat <- matrix(NA, nrow = nrow(coefs_join[[1]]$`HR-proficient`), ncol = length(coefs_join))
  rownames(coef.mat) <- rownames(coefs_join[[1]]$`HR-proficient`)
  for (i in 1:length(coefs_join)) {
    coef.mat[,i] <- coefs_join[[i]]$`HR-proficient`[,1] != 0
  }
  coef.mat <- coef.mat[-1,]
  genes.include <- rownames(coef.mat)[apply(coef.mat, 1, sum) == 1000]
  
  print(paste0('Number of genes in model ', model, ': ', length(genes.include)))
  
  ## Create centroids
  centroid.model <- data.frame(
    HRD = apply(Z.tumor_training[ann_tcga_train$HRD == 'HRD',genes.include], 2, mean),
    HR_proficient = apply(Z.tumor_training[ann_tcga_train$HRD == 'HR-proficient',genes.include],2,mean),
    BRCA1 = apply(Z.tumor_training[ann_tcga_train$group == 'BRCA1',genes.include],2,mean),
    BRCA2 = apply(Z.tumor_training[ann_tcga_train$group == 'BRCA2',genes.include],2,mean),
    HRD_BRCApos = apply(Z.tumor_training[ann_tcga_train$group == 'HRD_BRCA+',genes.include],2,mean),
    HR_BRCA_proficient = apply(Z.tumor_training[ann_tcga_train$group == 'HR-proficient',genes.include],2,mean)
  )
  
  signature.centroid.list[[model]] <- centroid.model
  
}

save(signature.centroid.list, file = '~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')





```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/TranscriptionalSignatureDevelopment/MultinomialElasticNet_alpha0.25_1to100.R`:

```````R
#####
## Multinomial Elastic Net Regression to classify BRCA1-/-, BRCA2-/-, HRD, HR-proficient TCGA-BRCA samples
#####

setwd('~/TranscriptionalSignatures/Revisions')

# Load libraries
library(glmnet)

# Load in training data:
#   TCGA-BRCA, expression deconvolution by BayesPrism, should be followed by low-gene removal and log2-normalisation
#   Genes remain if they have >1 count in >70% samples
#   HRD/BRCA status
load('Data/TCGA_BRCA.ExprDeconvolution_270723_p0.79_training.Rdata')
rownames(Z.tumor_training) <- substr(rownames(Z.tumor_training), 1, 12)
#index.keep <- apply(Z.tumor_training, 2, function(x) sum(x > 1)/length(x) > 0.7)
#Z.tumor_training <- Z.tumor_training[,index.keep]
#Z.tumor_training <- log2(Z.tumor_training + 1)

# Organise output variable
load('Data/TCGA_HRDclassification_BRCAannotation_HRD0.79.Rdata')
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

# Match output to training data
samples.intersect <- intersect(rownames(Z.tumor_training), rownames(ann_tcga))

ann_tcga <- ann_tcga[samples.intersect, ]
ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

Z.tumor_training <- Z.tumor_training[samples.intersect, ]


# Separate out expression data from BRCA_defect status
input.x <- Z.tumor_training
input.y <- ann_tcga$group

# Conduct 500 iterations of 10-fold cross validation:
#   (500 iteration ~ 2 days)
#   For each iteration, calculate coefficients using cv$lambda.min
#     (lambda coefficient which gives smallest error)
#   set seed inside loop

iterations <- 1:100
coefs <- list()

for (i in iterations) {
  
  # Set seed with each iteration
  set.seed(123*i)
  
  print(paste0('Iteration ', i, ' of ', iterations[length(iterations)], '...', Sys.time()))
  
  # Conduct grouped multinomial 10-fold cross validation
  cv <- cv.glmnet(x = input.x, y = input.y, family = 'multinomial',
                  alpha = 0.25, type.multinomial = 'grouped')
  coefs[[i]] <- coef(cv$glmnet.fit, s = cv$lambda.min)
  
}

save(coefs, file = 'Results/CVcoefficients_MultiElasticNet_alpha0.25_p0.79_iter1to100.Rdata')

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/TranscriptionalSignatureDevelopment/TCGA_BRCA.RNASeq_prep.R`:

```````R
#####
## Establishing training data and running BayesPrism for Expression Deconvolution of TCGA-BRCA
#####

setwd('~/Daniel/TCGAExpressionDeconvolution/')

# Load libraries
library(Seurat)
library(dplyr)
library(SummarizedExperiment)
library(caret)
library(BayesPrism)

# 1. Load Qian et al. 2020 data and create matrix

bc.data <- Read10X(data.dir = 'Data/Qian2020/BC_counts/')
bc_sc <- as.matrix(bc.data)
bc_sc <- t(bc_sc)

# 2. Match cell.type.labels from Qian et al. metadata
bc_met <- read.csv('Data/Qian2020/2103-Breastcancer_metadata.csv.gz')
bc_cellMeta <- bc_met$CellType

# 3a. Load TCGA bulk data obtained from TCGAbiolinks

load('Data/TCGA_BRCA.counts.SE.Rdata')
tcga_brca.rnaseq <- tcga_brca.rnaseq[!duplicated(rowData(tcga_brca.rnaseq)$gene_name), ]
tcga_brca.rnaseq <- tcga_brca.rnaseq[,!duplicated(tcga_brca.rnaseq$patient)]

# 3b. Match to HRD/BRCA-status output and remove PALB2/RAD51C defects
load('Data/TCGA_HRDclassification_BRCAannotation.Rdata')
# ann_tcga$patient <- rownames(ann_tcga)
# ann_tcga$HRD <- sapply(ann_tcga$HRD_prob, function(x) ifelse(x >= 0.79, 'HRD', 'HR-proficient'))
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga$HRD_BRCAstatus <- ann_tcga$BRCA_status
ann_tcga$HRD_BRCAstatus[ann_tcga$HRD == 'HRD' &
                          ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$HRD_BRCAstatus[ann_tcga$HRD_BRCAstatus == 'none'] <- 'HR-proficient'

ann_tcga <- ann_tcga[!(ann_tcga$HRD_BRCAstatus %in% c('PALB2', 'RAD51C')), ]

# 3c. Separate into training and testing (preserving class proportions), and save both
patients.intersect <- intersect(ann_tcga$Patient, tcga_brca.rnaseq$patient)
ann_tcga <- ann_tcga[match(patients.intersect, ann_tcga$patient), ]
tcga_brca.rnaseq <- tcga_brca.rnaseq[, match(patients.intersect, tcga_brca.rnaseq$patient)]

set.seed(1234)
inTrain <- createDataPartition(y = ann_tcga$HRD_BRCAstatus, p = 2/3, list = FALSE)
inTrain <- 1:ncol(tcga_brca.rnaseq) %in% inTrain
tcga_brca.rnaseq_train <- tcga_brca.rnaseq[, inTrain]
tcga_brca.rnaseq_test <- tcga_brca.rnaseq[, !inTrain]

tcga_brca.rnaseq_train$HRD_BRCAstatus <- ann_tcga$HRD_BRCAstatus[inTrain]
tcga_brca.rnaseq_test$HRD_BRCAstatus <- ann_tcga$HRD_BRCAstatus[!inTrain]

save(tcga_brca.rnaseq_train, file = 'Data/TCGA_BRCA.counts.SE_050823_training.Rdata')
save(tcga_brca.rnaseq_test, file = 'Data/TCGA_BRCA.counts.SE_050823_testing.Rdata')

bc_bulk <- t(assay(tcga_brca.rnaseq_train))
colnames(bc_bulk) <- rowData(tcga_brca.rnaseq_train)$gene_name

# 4. Preprocessing of Qian et al.

# Plot single cell outliers
sc.stat <- plot.scRNA.outlier(
  input = bc_sc,
  cell.type.labels = bc_met$CellType,
  species = 'hs',
  return.raw = FALSE,
  pdf.prefix = 'Qian2020_outlierPlot'
)

# Filter genes expressed in <2% cancer cells
bc_sc_cancer <- bc_sc[bc_cellMeta == 'Cancer', ]
index.gene_propCancer <- apply(bc_sc_cancer, 2, function(x) mean(x>0))
bc_sc <- bc_sc[,index.gene_propCancer > 0.02]

# Filter outlier genes from scRNA-seq data
bc_sc.filtered <- cleanup.genes(
  input = bc_sc,
  input.type = 'count.matrix',
  species = 'hs',
  gene.group = c('Rb', 'Mrp', 'other_Rb', 'chrM', 'MALAT1', 'chrX', 'chrY')
)
dim(bc_sc.filtered)

# Check concordance of varying gene types between scRNA-seq and bulk
plot.bulk.vs.sc(sc.input = bc_sc.filtered,
                bulk.input = bc_bulk)

# Subset protein coding genes, since these are most concordant and 
#   computation can be sped up
bc_sc.filtered.pc <- select.gene.type(
  input = bc_sc.filtered,
  gene.type = 'protein_coding'
)

# Construct a prism object
myPrism <- new.prism(
  reference = bc_sc.filtered.pc,
  mixture = bc_bulk,
  input.type = 'count.matrix',
  cell.type.labels = bc_cellMeta,
  cell.state.labels = bc_cellMeta,
  key = 'Cancer',
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

# Run BayesPrism
Sys.time(); tcga_brca.bayesPrism <- run.prism(prism = myPrism, n.cores = 20); Sys.time()

save(tcga_brca.bayesPrism, file = 'TCGA_BRCA.BayesPrism_050823_p0.79_training.Rdata')

# Extract cancer-specific expression
Z.tumor_training <- get.exp(
  bp = tcga_brca.bayesPrism, 
  state.or.type = 'type', cell.name = 'Cancer'
)
save(Z.tumor_training, file = 'TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/Bassez2021_preprocessing.R`:

```````R
#####
## Preprocessing Qian et al 2020 BC dataset
#####

setwd('~/Data/scRNASeq/Bassez2021/')

# Load libraries
library(Seurat)
library(dplyr)
library(cowplot)

# Load Qian et al. 2020 data
bc.data <- readRDS('1863-counts_cells_cohort1.rds')

# Extract pre-treatment cells
anno <- read.csv('1872-BIOKEY_metaData_cohort1_web.csv', header = TRUE)
anno <- anno[anno$timepoint == 'Pre', ]

bc.data <- bc.data[,anno$Cell]

# Initialize the Seurat object with the raw (non-normalized) data
#   init: 33694 genes, 44024 cells

## 8.3 Filtering low-quality cells
counts_per_cell <- Matrix::colSums(bc.data)
counts_per_gene <- Matrix::rowSums(bc.data)
genes_per_cell <- Matrix::colSums(bc.data>0)
cells_per_gene <- Matrix::rowSums(bc.data>0)

# 8.3.1 Summary counts for genes and cells
hist(log10(counts_per_cell+1), main='counts per cell', col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')
hist(log10(cells_per_gene+1), main='cells per gene', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat'); title('counts vs genes per cell')

# 8.3.2 Plot cells ranked by their number of detected genes
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')


## 8.4 Beginning with Seurat

# 8.4.1 Creating a seurat object

# Keep genes expressed in >= 3 cells (~.1% of the data).
# Keep all cells with at least 200 detected genes
#   now: 26040 genes, 38735 cells
seurat <- CreateSeuratObject(counts = bc.data, 
                             min.cells = 3, min.features = 200,
                             project = '10x_bc', assay = 'RNA')

## Load metadata
# anno <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.data <- seurat@meta.data
all(rownames(meta.data) == anno$Cell)
seurat$CellFromTumor <- TRUE
seurat$PatientNumber <- anno$patient_id
seurat$CellType <- anno$cellType

# # Keep only cells from tumours
# seurat <- subset(x = seurat, subset = CellFromTumor == 'TRUE')

## 8.5 Preprocessing Step 1: Filter out low-quality cells (in addition to 1. 200 minimum features)

# Common metric for judging damaged cells: relative expression of mitochondrially derived genes
#   Tell-tale sign of cell stress: widespread RNA degradation after apoptosis

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = '^MT-', x = rownames(x = seurat@assays$RNA@data), value=TRUE)
percent.mito <- Matrix::colSums(seurat@assays$RNA@data[mito.genes, ])/Matrix::colSums(seurat@assays$RNA@data)
seurat <- AddMetaData(object = seurat, metadata = percent.mito, col.name = 'percent.mito')
VlnPlot(object = seurat, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito'))

# Plot correlations RNA counts and other features
par(mfrow=c(1,2))
FeatureScatter(object = seurat, feature1 = 'nCount_RNA', feature2 = 'percent.mito')
FeatureScatter(object = seurat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

# Filter out cells with unique gene counts > 6,000 or mitochondrial content > 15%
seurat <- subset(x = seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                   percent.mito > -Inf & percent.mito < .15)

## 8.6.1 Preprocessing Step 2: Expression normalization
seurat <- NormalizeData(object = seurat, normalization.method = 'LogNormalize',
                        scale.factor = 10000)

# Save scRNA-seq data
expr.data_bassez2021 <- as.matrix(GetAssayData(seurat, slot = 'data'))
save(expr.data_bassez2021, file = 'exprData_Bassez2021.Rdata')

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/Qian2020_jointComparison.R`:

```````R
#####
## Check gene dropout in Qian et al 2020 BC dataset
#####

setwd('~/Data/scRNASeq/Qian2020/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)

# Load expression and metadata and subset for cancer cells
load('exprData_Qian2020.Rdata')

meta.qian <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.qian <- meta.qian[match(colnames(expr.data_qian2020), meta.qian$Cell), ]
expr.cancer <- expr.data_qian2020[,meta.qian$CellType == 'Cancer']

# Load signature centroids
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')

# load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives.Rdata')
# names(signature_alternative.centroid.list) <- paste0('Alternative_Sig_',names(signature_alternative.centroid.list))
# 
# signature.centroid.list <- c(signature.centroid.list, signature_alternative.centroid.list)
# rm(signature_alternative.centroid.list)

stats.df <- data.frame(
  Model = names(signature.centroid.list),
  prop_nonZeroCells = NA, 
  median_GenesExpressedPerCell = NA, mean_GenesExpressedPerCell = NA
)

for (i in 1:length(signature.centroid.list)) {
  
  print(names(signature.centroid.list)[i])
  
  sig.i <- signature.centroid.list[[i]]
  
  genes.intersect <- intersect(rownames(sig.i), rownames(expr.cancer))
  
  sig.i <- sig.i[genes.intersect, ]
  expr.cancer.i <- expr.cancer[genes.intersect, ]
  
  # Collate relevant results
  stats.df$prop_nonZeroCells[i] <- mean(apply(expr.cancer.i, 2, sum) > 0)
  stats.df$median_GenesExpressedPerCell[i] <- median(apply(expr.cancer.i, 2, function(x) sum(x>0)))
  stats.df$mean_GenesExpressedPerCell[i] <- mean(apply(expr.cancer.i, 2, function(x) sum(x>0)))
  
  
}

# stats.df$group <- sapply(stats.df$Model, function(x) strsplit(x,split='_')[[1]][1])

library(ggplot2)
library(tidyr)

# g1 <- ggplot(stats.df, aes(x = Model, y = prop_nonZeroCells, fill = Model)) + geom_bar(stat ='identity') + theme(axis.text.x = element_blank())
# g2 <- ggplot(stats.df, aes(x = Model, y = median_GenesExpressedPerCell, fill = Model)) + geom_bar(stat ='identity') + theme(axis.text.x = element_blank())
# g3 <- ggplot(stats.df, aes(x = Model, y = mean_GenesExpressedPerCell, fill = Model)) + geom_bar(stat ='identity') + theme(axis.text.x = element_blank())
# 
# ggarrange(plotlist = list(g1,g2,g3))

stats.df_plot <- stats.df %>%
  pivot_longer(cols = -Model, names_to = 'Measure')
g_dropout <- ggplot(stats.df_plot, aes(x = Model, y = value, fill = Model)) + 
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(legend.position = 'top',
        axis.text.x = element_blank()) +
  facet_wrap(~Measure, scales = 'free')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_QianSignatureComparison.pdf',
       plot = g_dropout, width = 8, height = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/Bassez2021_analysis.R`:

```````R
#####
## Apply signatures to Bassez et al 2021 BC dataset
#####

setwd('~/Data/scRNASeq/Bassez2021/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)

# Load expression and metadata and subset for cancer cells
load('exprData_Bassez2021.Rdata')

meta.bassez <- read.csv('1872-BIOKEY_metaData_cohort1_web.csv', header = TRUE)
meta.bassez <- meta.bassez[match(colnames(expr.data_bassez2021), meta.bassez$Cell), ]
expr.cancer <- expr.data_bassez2021[,meta.bassez$cellType == 'Cancer_cell']

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25

genes.intersect <- intersect(rownames(centroid.hrd), rownames(expr.cancer))
centroid.hrd <- centroid.hrd[genes.intersect, ]
expr.cancer <- expr.cancer[genes.intersect, ]

nonZero_genes <- apply(expr.cancer, 2, function(x) sum(x>0))
hist(nonZero_genes, breaks = 50)

# Calculate HRD scores across cells
hrdScores_bassez2021 <- data.frame(
  Cell = colnames(expr.cancer),
  Sample = sapply(colnames(expr.cancer), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(expr.cancer, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.cancer, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrdScores_bassez2021 <- hrdScores_bassez2021[!is.na(hrdScores_bassez2021$HRD), ]
hrdScores_bassez2021$HRD_score <- hrdScores_bassez2021$HRD - hrdScores_bassez2021$HR_proficient

# Prep for CellphoneDB analysis
cpdb.meta <- data.frame(
  Cell = meta.bassez$Cell,
  CellType = meta.bassez$cellType
)
hrdScores_bassez2021$HRD_group <- ifelse(hrdScores_bassez2021$HRD_score > 0, 'HRD', 'HR-proficient')
cpdb.meta <- merge(x = cpdb.meta, y = hrdScores_bassez2021[,c('Cell','HRD_group')], all.x = TRUE)
cpdb.meta$HRD_group[is.na(cpdb.meta$HRD_group)] <- ''
cpdb.meta$cell_type <- apply(cpdb.meta, 1, function(x) paste0(x[2],x[3],collapse = '_'))
cpdb.meta <- cpdb.meta[cpdb.meta$cell_type != 'Cancer_cell', ]

write.table(cpdb.meta[,c('Cell','cell_type')], file = '~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Bassez2021/Bassez2021_meta.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')

cpdb.expr <- as.data.frame(expr.data_bassez2021)
cpdb.expr <- cpdb.expr[,cpdb.meta$Cell]
cpdb.expr$Gene <- rownames(cpdb.expr)
cpdb.expr <- cpdb.expr[,c(ncol(cpdb.expr),1:(ncol(cpdb.expr)-1))]

write.table(cpdb.expr, file = '~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Bassez2021/Bassez2021_counts.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/Bassez2021_HRDprofiling.R`:

```````R
#####
## Apply signatures to bassez et al 2021 BC dataset
#####

setwd('~/Data/scRNASeq/Bassez2021/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)
library(wesanderson)
library(Seurat)

# Load expression and metadata and subset for cancer cells
load('exprData_Bassez2021.Rdata')

meta.bassez <- read.csv('1872-BIOKEY_metaData_cohort1_web.csv', header = TRUE)
meta.bassez <- meta.bassez[match(colnames(expr.data_bassez2021), meta.bassez$Cell), ]

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25
expr.hrd <- expr.data_bassez2021[rownames(centroid.hrd), meta.bassez$Cell]

# Gene inclusion plotting
expr.nonZero <- data.frame(
  Cell = colnames(expr.hrd), 
  Sample = sapply(colnames(expr.hrd), function(x) paste0(strsplit(x,split='_')[[1]][1:2], collapse='_')),
  prop_GenesExpressed = apply(expr.hrd, 2, function(x) 100*mean(x>0))
)
g_nonZero <- ggplot(expr.nonZero, aes(x = prop_GenesExpressed, fill = Sample)) +
  geom_density(alpha = 0.4) +
  theme_minimal() + theme(legend.position = 'none') +
  geom_vline(xintercept = mean(expr.nonZero$prop_GenesExpressed), col = 'red', linetype = 'dashed') +
  xlab('% Genes Expressed / Cell')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_BassezNonZeroGenes.pdf',
       plot = g_nonZero, width = 5.35, height = 4.79)

expr.nonZero_summary <- expr.nonZero %>%
  group_by(Sample) %>% summarise(meanExpression = mean(prop_GenesExpressed))
g_nonZeroSummary <- ggplot(expr.nonZero_summary, aes(x = meanExpression)) +
  geom_histogram(fill = 'lightblue', color = 'darkblue') +
  theme_minimal() + xlab('Mean % Genes Expressed Across Cells / Sample')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_BassezNonZeroGeneSummary.pdf',
       plot = g_nonZeroSummary, width = 5.35, height = 4.79)

# Calculate HRD scores across cells
hrdScores_bassez2021 <- data.frame(
  Cell = colnames(expr.hrd),
  CellType = meta.bassez$cellType,
  Sample = sapply(colnames(expr.hrd), function(x) paste0(strsplit(x,split='_')[[1]][1:2],collapse = '_')),
  HRD = apply(expr.hrd, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.hrd, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrdScores_bassez2021 <- hrdScores_bassez2021[!is.na(hrdScores_bassez2021$HRD), ]
hrdScores_bassez2021$HRD_score <- hrdScores_bassez2021$HRD - hrdScores_bassez2021$HR_proficient

# Factor results in decreasing average HRD score
bassez2021_hrdSummary <- hrdScores_bassez2021[hrdScores_bassez2021$CellType == 'Cancer_cell',] %>%
  group_by(Sample) %>% summarise(mean_HRD = mean(HRD_score)) %>%
  arrange(desc(mean_HRD))
hrdScores_bassez2021$Sample <- factor(hrdScores_bassez2021$Sample,
                                    levels = bassez2021_hrdSummary$Sample)
hrdScores_bassez2021$Cancer <- hrdScores_bassez2021$CellType == 'Cancer_cell'

mu <- hrdScores_bassez2021 %>%
  group_by(Sample, Cancer) %>% summarise(medianHRD = median(HRD_score))

g_tme <-ggplot(mu, aes(Cancer, medianHRD, fill=Cancer)) +
  geom_boxplot() +
  geom_point() + geom_line(aes(group = Sample)) +
  theme_minimal() +
  theme(legend.position = 'none') + ylab('median HRD score') +
  scale_fill_manual(values = wes_palette('Darjeeling2')) +
  ggtitle('Bassez et al.')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDinTME_Bassez.pdf',
       plot = g_tme, width = 4, height = 4)

# Map clinical data and proportion of cells with HRD > 0
bassez2021_hrd_props <- hrdScores_bassez2021[hrdScores_bassez2021$Cancer,] %>%
  group_by(as.character(Sample)) %>% summarise(prop_HRD = 100*mean(HRD_score > 0))
names(bassez2021_hrd_props)[1] <- 'Sample'
bassez2021_hrd_props$id <- sapply(bassez2021_hrd_props$Sample,
                                  function(x) as.numeric(strsplit(x,split='_')[[1]][2]))
bassez2021_hrd_props <- bassez2021_hrd_props[order(bassez2021_hrd_props$id), ]

bassez2021_hrd_props$BC_subtype = factor(c(
  'ER-HER2+','TNBC','TNBC','TNBC','TNBC',
  'ER-HER2+','TNBC','ER-HER2+','ER+HER2+','TNBC',
  'TNBC','ER+HER2-','ER+HER2-','ER+HER2-','TNBC',
  'TNBC','ER+HER2-','ER+HER2-','TNBC','ER+HER2-',
  'ER+HER2-','ER+HER2-','TNBC','ER+HER2-','TNBC',
  'ER+HER2+','TNBC','ER+HER2-','ER+HER2-','ER+HER2-','ER+HER2-'
), levels = c('ER+HER2+','ER+HER2-','ER-HER2+','TNBC'))
bassez2021_hrd_props <- bassez2021_hrd_props[order(bassez2021_hrd_props$prop_HRD), ]

bassez2021_hrd_props$Sample <- factor(bassez2021_hrd_props$Sample,
                                    levels = bassez2021_hrd_props$Sample)

g_bassez2021_HRDprops <- ggplot(data = bassez2021_hrd_props, aes(x = Sample, y = prop_HRD, fill = BC_subtype)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab('% HRD cells') +
  scale_fill_manual(values = c(wes_palette('Moonrise3')[c(5,2:4)]))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Bassez2021_propHRDscores.pdf',
       plot = g_bassez2021_HRDprops, width = 4, height = 4)

# Plot density plots
hrdScores_bassez2021 <- merge(x = hrdScores_bassez2021, y = bassez2021_hrd_props[,c('Sample','BC_subtype')])
g_densities <- ggplot(hrdScores_bassez2021[hrdScores_bassez2021$Cancer, ], aes(x = HRD_score, fill = BC_subtype)) +
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c(wes_palette('Moonrise3')[c(5,2:4)])) +
  theme(legend.position = 'top') +
  facet_wrap(~Sample, ncol = 4)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDdensities_Bassez.pdf',
       plot = g_densities, height = 8, width = 5.35)

## Bassez2021 UMAP plotting

# Create new Seurat object with normalised data
bassez_seurat <- CreateSeuratObject(counts = expr.data_bassez2021[, meta.bassez$Cell[meta.bassez$cellType == 'Cancer_cell']],
                                    project = 'bassez2021', min.cells = 3, min.features = 200)

# Identify highly variable features and scale data
bassez_seurat <- FindVariableFeatures(bassez_seurat, selection.method = 'vst',
                                    nfeatures = 2000)
bassez_seurat <- ScaleData(bassez_seurat)

# Perform linear dimensionality reduction, clustering, and UMAP
bassez_seurat <- RunPCA(bassez_seurat, features = VariableFeatures(object = bassez_seurat))
bassez_seurat <- FindNeighbors(bassez_seurat, dims = 1:10)
bassez_seurat <- FindClusters(bassez_seurat, resolution = 0.5)
bassez_seurat <- RunUMAP(bassez_seurat, dims = 1:10)

# UMAP plotting
umap.data <- as.data.frame(bassez_seurat@reductions$umap@cell.embeddings)
rownames(hrdScores_bassez2021) <- hrdScores_bassez2021$Cell
umap.data <- merge(x = umap.data, y = hrdScores_bassez2021[,c('Sample', 'HRD_score')], by=0)
umap.data <- merge(x = umap.data, y = bassez2021_hrd_props[,c('Sample','BC_subtype')], by = 'Sample')

# Plot UMAP coordinates coloured by breast cancer subtype
g_bassez2021_BC_subtype <- ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = BC_subtype)) +
  geom_point(size = .2) + theme_minimal() +
  scale_color_manual(values = c(wes_palette('Moonrise3')[c(5,2:4)])) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top') +
  guides(color = guide_legend(override.aes = list(size = 5),
                              nrow = 2))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/bassez2021_UMAP_BC_subtype.pdf',
       plot = g_bassez2021_BC_subtype, width = 4, height = 4)

# Plot UMAP coordinates coloured by HRD score
g_bassez2021_hrdScores <- ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = HRD_score)) +
  geom_point(size = .2) + theme_minimal() +
  scale_color_gradientn(colors = c(rep(wes_palette('GrandBudapest1')[1],2),
                                   'gray90',
                                   rep(wes_palette('GrandBudapest1')[2],2)),
                        values = c(0,.3,
                                   abs(min(umap.data$HRD_score))/(abs(min(umap.data$HRD_score)) + max(umap.data$HRD_score)),
                                   .7,1)) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/bassez2021_UMAP_HRDscores.pdf',
       plot = g_bassez2021_hrdScores, width = 4, height = 4)


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/Qian2020_analysis.R`:

```````R
#####
## Apply signatures to Qian et al 2020 BC dataset
#####

setwd('~/Data/scRNASeq/Qian2020/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)
library(wesanderson)

# Load expression and metadata and subset for cancer cells
load('exprData_Qian2020.Rdata')

meta.qian <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.qian <- meta.qian[match(colnames(expr.data_qian2020), meta.qian$Cell), ]
expr.cancer <- expr.data_qian2020[,meta.qian$CellType == 'Cancer']

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25

genes.intersect <- intersect(rownames(centroid.hrd), rownames(expr.cancer))
centroid.hrd <- centroid.hrd[genes.intersect, ]
expr.cancer <- expr.cancer[genes.intersect, ]

# nonZero_genes <- apply(expr.cancer, 2, function(x) sum(x>0))
# hist(nonZero_genes, breaks = 50)

# Calculate HRD scores across cells
hrdScores_qian2020 <- data.frame(
  Cell = colnames(expr.cancer),
  Sample = sapply(colnames(expr.cancer), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(expr.cancer, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.cancer, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrdScores_qian2020 <- hrdScores_qian2020[!is.na(hrdScores_qian2020$HRD), ]
hrdScores_qian2020$HRD_score <- hrdScores_qian2020$HRD - hrdScores_qian2020$HR_proficient

# Factor results in decreasing average HRD score
qian2020_hrdSummary <- hrdScores_qian2020 %>%
  group_by(Sample) %>% summarise(mean_HRD = mean(HRD_score)) %>%
  arrange(desc(mean_HRD))
hrdScores_qian2020$Sample <- factor(hrdScores_qian2020$Sample,
                                    levels = qian2020_hrdSummary$Sample)

# Plot density plots
ggplot(hrdScores_qian2020, aes(x = HRD_score)) +
  geom_density(fill = 'lightblue') + facet_wrap(~Sample)

# Compare BRCA1-/- vs Luminal A samples
hrdScores_plot <- hrdScores_qian2020
hrdScores_plot$label <- NA
hrdScores_plot$label[hrdScores_plot$Sample == 'sc5rJUQ033'] <- 'BRCA1-/- TNBC'
hrdScores_plot$label[hrdScores_plot$Sample == 'sc5rJUQ064'] <- 'Lum A-like'
hrdScores_plot <- hrdScores_plot[!is.na(hrdScores_plot$label), ]

g_plotTwo <- ggplot(hrdScores_plot, aes(x = HRD_score, fill = label)) +
  geom_density(alpha = 0.4) + theme_minimal() +
  theme(legend.position = 'top',
        legend.title = element_blank()) +
  ylab('density(cells)') +
  scale_fill_manual(values = c(wes_palette('Moonrise3')[c(1,5)]))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Qian2020_BRCA1_LumA.pdf',
       plot = g_plotTwo, width = 4, height = 3)

# Prep for CellphoneDB analysis
cpdb.meta <- data.frame(
  Cell = meta.qian$Cell,
  CellType = meta.qian$CellType
)
hrdScores_qian2020$HRD_group <- ifelse(hrdScores_qian2020$HRD_score > 0, 'HRD', 'HR-proficient')
cpdb.meta <- merge(x = cpdb.meta, y = hrdScores_qian2020[,c('Cell','HRD_group')], all.x = TRUE)
cpdb.meta$HRD_group[is.na(cpdb.meta$HRD_group)] <- ''
cpdb.meta$cell_type <- apply(cpdb.meta, 1, function(x) paste0(x[2],x[3],collapse = '_'))
cpdb.meta <- cpdb.meta[cpdb.meta$cell_type != 'Cancer', ]

write.table(cpdb.meta[,c('Cell','cell_type')], file = '~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Qian2020/Qian2020_meta.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')

cpdb.expr <- as.data.frame(expr.data_qian2020)
cpdb.expr <- cpdb.expr[,cpdb.meta$Cell]
cpdb.expr$Gene <- rownames(cpdb.expr)
cpdb.expr <- cpdb.expr[,c(ncol(cpdb.expr),1:(ncol(cpdb.expr)-1))]

write.table(cpdb.expr, file = '~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Qian2020/Qian2020_counts.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/Chung2017_analysis.R`:

```````R
#####
## Apply transcriptional signature to Chung et al. 2017 bulk and scRNAseq
#####

# Load libraries
library(dplyr)
library(ggpubr)

# Load signature
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25

# Load expression data, match with signature, and log2-normalise
expr.chung <- read.table('~/Data/scRNASeq/Chung2017/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt',h=T)

genes.intersect <- intersect(rownames(centroid.hrd), expr.chung$gene_name)

centroid.hrd <- centroid.hrd[genes.intersect, ]
expr.chung <- expr.chung[match(genes.intersect, expr.chung$gene_name), ]
rownames(expr.chung) <- expr.chung$gene_name
expr.chung <- expr.chung[,-c(1:3)]

expr.chung <- log2(expr.chung + 1)

# Separate into bulk and scRNAseq
expr.chung_bulk <- expr.chung[,grepl(pattern = 'Pooled', colnames(expr.chung))]
expr.chung_sc <- expr.chung[,15:ncol(expr.chung)]

# Extract only tumour cells from single cell data
chung_info <- read.table('~/Data/scRNASeq/Chung2017/GSE75688_final_sample_information.txt',h=T)
cells.tumor <- chung_info$sample[chung_info$type == 'SC' & chung_info$index == 'Tumor']

expr.chung_sc <- expr.chung_sc[,cells.tumor]

# Calculate HRD scores
hrd_scores.bk <- data.frame(
  sample = sapply(colnames(expr.chung_bulk), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(expr.chung_bulk, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.chung_bulk, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrd_scores.bk$HRD_score_bulk <- hrd_scores.bk$HRD - hrd_scores.bk$HR_proficient

hrd_scores.sc <- data.frame(
  cell = colnames(expr.chung_sc),
  HRD = apply(expr.chung_sc, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.chung_sc, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrd_scores.sc$HRD_score_sc <- hrd_scores.sc$HRD - hrd_scores.sc$HR_proficient
hrd_scores.sc$sample <- sapply(hrd_scores.sc$cell, function(x) strsplit(x,split='_')[[1]][1])

# Match bulk and single-cell HRD scores
hrd_scores.sc_summary <- hrd_scores.sc %>%
  group_by(sample) %>% summarise(mean_HRD_score_sc = mean(HRD_score_sc))

hrd_scores.df <- merge(x = hrd_scores.bk[,c('sample','HRD_score_bulk')],
                       y = hrd_scores.sc_summary)

# Plot results
g_chung <- ggplot(hrd_scores.df, aes(x = HRD_score_bulk, y = mean_HRD_score_sc)) +
  geom_point() + geom_smooth(method = 'lm', col = 'darkred') + stat_cor() +
  theme_minimal() +
  xlab('Bulk HRD') + ylab('Mean single-cell HRD')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Chung_BulkSingleCell.pdf',
       plot = g_chung, width = 4, height = 3)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/Qian2020_preprocessing.R`:

```````R
#####
## Preprocessing Qian et al 2020 BC dataset
#####

setwd('~/Data/scRNASeq/Qian2020/')

# Load libraries
library(Seurat)
library(dplyr)
library(cowplot)

# Load Qian et al. 2020 data
bc.data <- Read10X(data.dir = 
                     '~/Data/scRNASeq/Qian2020/export/BC_counts/')

# Initialize the Seurat object with the raw (non-normalized) data
#   init: 33694 genes, 44024 cells

## 8.3 Filtering low-quality cells
counts_per_cell <- Matrix::colSums(bc.data)
counts_per_gene <- Matrix::rowSums(bc.data)
genes_per_cell <- Matrix::colSums(bc.data>0)
cells_per_gene <- Matrix::rowSums(bc.data>0)

# 8.3.1 Summary counts for genes and cells
hist(log10(counts_per_cell+1), main='counts per cell', col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')
hist(log10(cells_per_gene+1), main='cells per gene', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat'); title('counts vs genes per cell')

# 8.3.2 Plot cells ranked by their number of detected genes
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')


## 8.4 Beginning with Seurat

# 8.4.1 Creating a seurat object

# Keep genes expressed in >= 3 cells (~.1% of the data).
# Keep all cells with at least 200 detected genes
#   now: 26040 genes, 38735 cells
seurat <- CreateSeuratObject(counts = bc.data, 
                             min.cells = 3, min.features = 200,
                             project = '10x_bc', assay = 'RNA')

## Load metadata
anno <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.data <- seurat@meta.data
all(rownames(meta.data) == anno$Cell)
seurat$CellFromTumor <- anno$CellFromTumor
seurat$PatientNumber <- anno$PatientNumber
seurat$CellType <- anno$CellType

# # Keep only cells from tumours
# seurat <- subset(x = seurat, subset = CellFromTumor == 'TRUE')

## 8.5 Preprocessing Step 1: Filter out low-quality cells (in addition to 1. 200 minimum features)

# Common metric for judging damaged cells: relative expression of mitochondrially derived genes
#   Tell-tale sign of cell stress: widespread RNA degradation after apoptosis

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = '^MT-', x = rownames(x = seurat@assays$RNA@data), value=TRUE)
percent.mito <- Matrix::colSums(seurat@assays$RNA@data[mito.genes, ])/Matrix::colSums(seurat@assays$RNA@data)
seurat <- AddMetaData(object = seurat, metadata = percent.mito, col.name = 'percent.mito')
VlnPlot(object = seurat, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito'))

# Plot correlations RNA counts and other features
par(mfrow=c(1,2))
FeatureScatter(object = seurat, feature1 = 'nCount_RNA', feature2 = 'percent.mito')
FeatureScatter(object = seurat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

# Filter out cells with unique gene counts > 6,000 or mitochondrial content > 15%
seurat <- subset(x = seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                   percent.mito > -Inf & percent.mito < .15)

## 8.6.1 Preprocessing Step 2: Expression normalization
seurat <- NormalizeData(object = seurat, normalization.method = 'LogNormalize',
                        scale.factor = 10000)

# Save scRNA-seq data
expr.data_qian2020 <- as.matrix(GetAssayData(seurat, slot = 'data'))
save(expr.data_qian2020, file = 'exprData_Qian2020.Rdata')

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/Qian2020_HRDprofiling.R`:

```````R
#####
## Apply signatures to Qian et al 2020 BC dataset
#####

setwd('~/Data/scRNASeq/Qian2020/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)
library(wesanderson)
library(Seurat)

# Load expression and metadata and subset for cancer cells
load('exprData_Qian2020.Rdata')

meta.qian <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.qian <- meta.qian[match(colnames(expr.data_qian2020), meta.qian$Cell), ]

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25
expr.hrd <- expr.data_qian2020[rownames(centroid.hrd), meta.qian$Cell]

# Gene inclusion plotting
expr.nonZero <- data.frame(
  Cell = colnames(expr.hrd), 
  Sample = sapply(colnames(expr.hrd), function(x) strsplit(x,split='_')[[1]][1]),
  prop_GenesExpressed = apply(expr.hrd, 2, function(x) 100*mean(x>0))
)
g_nonZero <- ggplot(expr.nonZero, aes(x = prop_GenesExpressed, fill = Sample)) +
  geom_density(alpha = 0.4) +
  theme_minimal() + theme(legend.position = 'none') +
  geom_vline(xintercept = mean(expr.nonZero$prop_GenesExpressed), col = 'red', linetype = 'dashed') +
  xlab('% Genes Expressed / Cell')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_QianNonZeroGenes.pdf',
       plot = g_nonZero)

expr.nonZero_summary <- expr.nonZero %>%
  group_by(Sample) %>% summarise(meanExpression = mean(prop_GenesExpressed))
g_nonZeroSummary <- ggplot(expr.nonZero_summary, aes(x = meanExpression)) +
  geom_histogram(fill = 'lightblue', color = 'darkblue') +
  theme_minimal() + xlab('Mean % Genes Expressed Across Cells / Sample')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_QianNonZeroGeneSummary.pdf',
       plot = g_nonZeroSummary)

# Calculate HRD scores across cells
hrdScores_qian2020 <- data.frame(
  Cell = colnames(expr.hrd),
  CellType = meta.qian$CellType,
  Sample = sapply(colnames(expr.hrd), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(expr.hrd, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.hrd, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrdScores_qian2020 <- hrdScores_qian2020[!is.na(hrdScores_qian2020$HRD), ]
hrdScores_qian2020$HRD_score <- hrdScores_qian2020$HRD - hrdScores_qian2020$HR_proficient

# Factor results in decreasing average HRD score
qian2020_hrdSummary <- hrdScores_qian2020[hrdScores_qian2020$CellType == 'Cancer',] %>%
  group_by(Sample) %>% summarise(mean_HRD = mean(HRD_score)) %>%
  arrange(desc(mean_HRD))
hrdScores_qian2020$Sample <- factor(hrdScores_qian2020$Sample,
                                    levels = qian2020_hrdSummary$Sample)
hrdScores_qian2020$Cancer <- hrdScores_qian2020$CellType == 'Cancer'

mu <- hrdScores_qian2020 %>%
  group_by(Sample, Cancer) %>% summarise(medianHRD = median(HRD_score))

g_tme <-ggplot(mu, aes(Cancer, medianHRD, fill=Cancer)) +
  geom_boxplot() +
  geom_point() + geom_line(aes(group = Sample)) +
  theme_minimal() +
  theme(legend.position = 'none') + ylab('median HRD score') +
  scale_fill_manual(values = wes_palette('Darjeeling2')) +
  ggtitle('Qian et al.')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDinTME_Qian.pdf',
       plot = g_tme, width = 4, height = 4)

# Map clinical data and proportion of cells with HRD > 0
qian2020_hrd_props <- hrdScores_qian2020[hrdScores_qian2020$Cancer,] %>%
  group_by(as.character(Sample)) %>% summarise(prop_HRD = 100*mean(HRD_score > 0))
qian2020_hrd_props$BC_subtype = factor(c(
  'HER2', 'TN', 'B1_TN', 'TN', 'TN', 'HER2',
  'TN', 'HER2', 'Lum_HER2', 'TN', 'TN', 'TN', 'LumB', 'LumA'
), levels = c('B1_TN','Lum_HER2','HER2',
              'TN','LumA','LumB'))
qian2020_hrd_props <- qian2020_hrd_props[order(qian2020_hrd_props$prop_HRD), ]

names(qian2020_hrd_props)[1] <- 'Sample'
qian2020_hrd_props$Sample <- factor(qian2020_hrd_props$Sample,
                                    levels = qian2020_hrd_props$Sample)

g_qian2020_HRDprops <- ggplot(data = qian2020_hrd_props, aes(x = Sample, y = prop_HRD, fill = BC_subtype)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab('% HRD cells') +
  scale_fill_manual(values = c(wes_palette('Moonrise3'),
                               wes_palette('Moonrise1')[1]))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Qian2020_propHRDscores.pdf',
       plot = g_qian2020_HRDprops, width = 4, height = 4)

# Plot density plots
hrdScores_qian2020 <- merge(x = hrdScores_qian2020, y = qian2020_hrd_props[,c('Sample','BC_subtype')])
g_densities <- ggplot(hrdScores_qian2020[hrdScores_qian2020$Cancer, ], aes(x = HRD_score, fill = BC_subtype)) +
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c(wes_palette('Moonrise3'),
                               wes_palette('Moonrise1')[1])) +
  theme(legend.position = 'top') +
  facet_wrap(~Sample, ncol = 4)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDdensities_Qian.pdf',
       plot = g_densities)

## Qian2020 UMAP plotting

# Create new Seurat object with normalised data
qian_seurat <- CreateSeuratObject(counts = expr.data_qian2020[, meta.qian$Cell[meta.qian$CellType == 'Cancer']],
                                  project = 'qian2020', min.cells = 3, min.features = 200)

# Identify highly variable features and scale data
qian_seurat <- FindVariableFeatures(qian_seurat, selection.method = 'vst',
                                    nfeatures = 2000)
qian_seurat <- ScaleData(qian_seurat)

# Perform linear dimensionality reduction, clustering, and UMAP
qian_seurat <- RunPCA(qian_seurat, features = VariableFeatures(object = qian_seurat))
qian_seurat <- FindNeighbors(qian_seurat, dims = 1:10)
qian_seurat <- FindClusters(qian_seurat, resolution = 0.5)
qian_seurat <- RunUMAP(qian_seurat, dims = 1:10)

# UMAP plotting
umap.data <- as.data.frame(qian_seurat@reductions$umap@cell.embeddings)
umap.data <- merge(x = umap.data, y = hrdScores_qian2020[,c('Sample', 'HRD_score')], by=0)
umap.data <- merge(x = umap.data, y = qian2020_hrd_props[,c('Sample','BC_subtype')], by = 'Sample')

# Plot UMAP coordinates coloured by breast cancer subtype
g_qian2020_BC_subtype <- ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = BC_subtype)) +
  geom_point(size = .2) + theme_minimal() +
  scale_color_manual(values = c(wes_palette('Moonrise3'),
                                wes_palette('Moonrise1')[1])) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top') +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Qian2020_UMAP_BC_subtype.pdf',
       plot = g_qian2020_BC_subtype, width = 4, height = 4)

# Plot UMAP coordinates coloured by HRD score
g_qian2020_hrdScores <- ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = HRD_score)) +
  geom_point(size = .2) + theme_minimal() +
  scale_color_gradientn(colors = c(rep(wes_palette('GrandBudapest1')[1],2),
                                   'gray90',
                                   rep(wes_palette('GrandBudapest1')[2],2)),
                        values = c(0,.3,
                                   abs(min(umap.data$HRD_score))/(abs(min(umap.data$HRD_score)) + max(umap.data$HRD_score)),
                                   .7,1)) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Qian2020_UMAP_HRDscores.pdf',
       plot = g_qian2020_hrdScores, width = 4, height = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/cpdb_interactionCounts.R`:

```````R
setwd('~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Results/')

library(tidyr)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(networkD3)
library(wesanderson)

qian <- read.delim('Bassez2021/significant_means.txt')
names(qian) <- gsub('HR.proficient','HR_proficient',names(qian))

# target_cancer.index <- which(sapply(names(qian), function(x) strsplit(x,split='[.]')[[1]][2])
#       %in% c('CancerHR_proficient','CancerHRD'))
qian <- qian[,c(2,13:ncol(qian))]

qian <- qian %>%
  pivot_longer(cols = -interacting_pair, names_to = 'cell_int', values_to = 'interaction')
qian$SOURCE <- sapply(qian$cell_int, function(x) strsplit(x,split='[.]')[[1]][1])
qian$TARGET <- sapply(qian$cell_int, function(x) strsplit(x,split='[.]')[[1]][2])

qian <- qian %>%
  group_by(SOURCE,TARGET) %>%
  summarise(interactions = sum(!is.na(interaction)))

qian_heatmap <- qian %>%
  pivot_wider(names_from = TARGET, values_from = interactions)
qian_heatmap <- as.data.frame(qian_heatmap)
rownames(qian_heatmap) <- qian_heatmap$SOURCE; qian_heatmap <- qian_heatmap[,-1]

col1 = 'dodgerblue4'; col2 = 'peachpuff'; col3 = 'deeppink4'
col.heatmap <- colorRampPalette(c(col1, col2, col3))(1000)
pheatmap(t(qian_heatmap), color = col.heatmap)

# cellType_order <- c('Fibroblast','EC','Myeloid','DC','Mast',
#                     'T_cell','B_cell','CancerHR_proficient','CancerHRD')
cellType_order <- c('Fibroblast','Endothelial_cell','Myeloid_cell','pDC','Mast_cell',
                    'T_cell','B_cell','Cancer_cellHR_proficient','Cancer_cellHRD')

qian_heatmap <- qian_heatmap[cellType_order, cellType_order]
pheatmap(t(qian_heatmap), color = col.heatmap,
         cluster_rows = FALSE, cluster_cols = FALSE,
         filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Bassez_Heatmap.pdf',
         width = 5, height = 5)

qian_sankey <- qian[!grepl(pattern = 'Cancer',qian$SOURCE) &
                      grepl(pattern = 'Cancer',qian$TARGET), ]
qian_sankey.Nodes <- data.frame(
  name = c(as.character(qian_sankey$SOURCE),
           as.character(qian_sankey$TARGET)) %>% unique()
)
qian_sankey$IDsource <- match(qian_sankey$SOURCE, qian_sankey.Nodes$name)-1
qian_sankey$IDtarget <- match(qian_sankey$TARGET, qian_sankey.Nodes$name)-1

# sankeyNetwork(Links = qian_sankey, Nodes = qian_sankey.Nodes,
#               Source = 'IDsource', Target = 'IDtarget',
#               Value = 'interactions', NodeID = 'name', fontSize = 20)

g_counts <- ggplot(qian_sankey, aes(x = TARGET, y = interactions, fill = TARGET)) + 
  geom_bar(stat = 'identity') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  facet_wrap(~SOURCE, scales = 'free', nrow=2)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Bassez_countBarplot.pdf',
       plot = g_counts, width = 6, height = 3.5)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SingleCellAnalysis/cpdb_individualInteractions.R`:

```````R
setwd('~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Results/')

library(dplyr)
library(tidyr)
library(ggpubr)
library(wesanderson)

cellType.of.interest <- 'Myeloid_cell'

# Process Qian
qian <- read.delim('Qian2020/significant_means.txt')
qian <- qian[,c('interacting_pair',
                paste(cellType.of.interest,'CancerHR.proficient',sep='.'),
                paste(cellType.of.interest,'CancerHRD',sep='.'))]
qian$status <- NA
qian$status[!is.na(qian[,2])] <- 'ONLY HR-proficient'
qian$status[!is.na(qian[,3])] <- 'ONLY HRD'
qian$status[!is.na(qian[,2]) &
              !is.na(qian[,3])] <- 'Both'
qian <- qian[!is.na(qian$status), ]
qian$Dataset <- 'Qian2020'

# Process Bassez
bassez <- read.delim('Bassez2021/significant_means.txt')
bassez <- bassez[,c('interacting_pair',
                    paste(cellType.of.interest,'Cancer_cellHR.proficient',sep='.'),
                    paste(cellType.of.interest,'Cancer_cellHRD',sep='.'))]
bassez$status <- NA
bassez$status[!is.na(bassez[,2])] <- 'ONLY HR-proficient'
bassez$status[!is.na(bassez[,3])] <- 'ONLY HRD'
bassez$status[!is.na(bassez[,2]) &
              !is.na(bassez[,3])] <- 'Both'
bassez <- bassez[!is.na(bassez$status), ]
bassez$Dataset <- 'Bassez2021'

names(bassez) <- names(qian)

# Merge
df <- rbind(qian,bassez)
names(df) <- c('interacting_pair','int_HRproficient','int_HRD','status','Dataset')
df 
df <- df %>%
  pivot_longer(cols = c(int_HRproficient, int_HRD),
               names_to = 'CancerStatus', values_to = 'Interaction')
df$group <- paste(df$Dataset, df$CancerStatus, sep = '_')
df$group <- factor(df$group,
                   levels = c('Qian2020_int_HRD','Bassez2021_int_HRD',
                              'Qian2020_int_HRproficient','Bassez2021_int_HRproficient'))

g_ints <- ggballoonplot(df, x = 'interacting_pair', y = 'group',
              size = 'Interaction', fill = 'status') +
  scale_fill_manual(values = c('grey90',wes_palette('GrandBudapest1')[1:2])) +
  ggtitle(paste0('Cell type: ',cellType.of.interest)) +
  geom_hline(yintercept = 2.5, color = 'red', linetype = 'dashed') +
  theme(legend.position = 'top')
ggsave(filename = paste0('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_CPDB_Interactions_',cellType.of.interest,'.pdf'),
       plot = g_ints, height = 4, width = 10)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/ISPY2_HRDscoring.R`:

```````R
#####
## Correlate HRD transcriptional score with pCR status in I-SPY2 trial
#####

# Load libraries
library(ggpubr)

## Organise Data

# Load I-SPY2 expression and clinical data from Puzstai et al. 2021
ispy2_expr <- read.table('~/Data/ClinicalData/ISPY2_Puzstai2021_expression.txt',h=T,row.names = 1)

ispy2_response <- read.csv('~/Data/ClinicalData/GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv')
ispy2_response <- ispy2_response[ispy2_response$Arm == 'durvalumab/olaparib', ]
ispy2_response$pCR.status[ispy2_response$pCR.status == -1] <- 0 # present in control arm

# Extract clinical arm from expression data
ispy2_response$ResearchID <- paste0('X',ispy2_response$ResearchID)
ispy2_expr <- ispy2_expr[,ispy2_response$ResearchID]

# Load signature centroids, extract ElasticNet_alpha0.25, and organise it with expression matrix
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
sig <- signature.centroid.list$ElasticNet_alpha0.25

genes.intersect <- intersect(rownames(sig), rownames(ispy2_expr))
sig <- sig[genes.intersect, ]
ispy2_expr <- ispy2_expr[genes.intersect, ]

## Calculate HRD scores and match with relevant clinical data
ispy2_hrd <- data.frame(
  ResearchID = colnames(ispy2_expr),
  HRD_score = apply(ispy2_expr, 2, function(x) cor(x,sig$HRD)) -
    apply(ispy2_expr, 2, function(x) cor(x,sig$HR_proficient))
)

ispy2_hrd <- merge(x = ispy2_hrd, y = ispy2_response[,c('ResearchID','pCR.status', 'PARPi7_sig.')])

# Format pCR status into responders vs non-responders
ispy2_hrd$Response <- sapply(ispy2_hrd$pCR.status, function(x)
  ifelse(x == 1, 'Responder', 'Non-responder'))
ispy2_hrd$Response <- factor(ispy2_hrd$Response,
                             levels = c('Non-responder','Responder'))

# Plot HRD score against pCR status
g_ispy2 <- ggboxplot(data = ispy2_hrd, x = 'Response', y = 'HRD_score',
          add = 'jitter', color = 'Response') +
  stat_compare_means() + scale_color_brewer(palette = 'Paired') +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) +
  ylab('HRD score')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure5/ISPY2_Response.pdf',
       plot = g_ispy2, width = 5, height = 5)
# ggsave(filename = '~/Projects/Thesis/Chapter 5/ISPY2_controlArm.pdf',
#        plot = g_ispy2, width = 5, height = 5)

# Plot HRD score against PARPi7 score
g_hrdvsparpi7 <- ggplot(data = ispy2_hrd, aes(x = PARPi7_sig., y = HRD_score)) +
  geom_point(aes(color = Response)) + 
  geom_smooth(method = 'lm', color = 'gray60') + stat_cor() +
  xlab('PARPi7 Signature Score') + ylab('HRD score') +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = 'top') +
  scale_color_brewer(palette = 'Paired')
# ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure5/ISPY2_HRDvsPARPi7.pdf',
#        plot = g_hrdvsparpi7, width = 5, height = 5)

# Plot PARPi7 score against pCR status
g_parpi7 <- ggboxplot(data = ispy2_hrd, x = 'Response', y = 'PARPi7_sig.',
                     add = 'jitter', color = 'Response') +
  stat_compare_means() + scale_color_brewer(palette = 'Paired') +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) +
  ylab('PARPi7 Signature Score')
# ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_ISPY2_PARPi7.pdf',
#        plot = g_parpi7, width = 5, height = 5)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/CCLE_jointComparisons.R`:

```````R
#####
## Scripts to compare centroid models against PARPi sensitivity in CCLE
#####

# Load libraries
library(data.table)
library(ggplot2)
library(tidyr)
library(ggpubr)

# Load CCLE expression data
setwd('~/Data/CCLE')

expr <- read.csv('Expression_Public_23Q2_subsetted.csv', row.names = 1)

# Load signature centroids (subset for ElasticNet_alpha0.25)
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
signature.centroid.list <- list(signature.centroid.list$ElasticNet_alpha0.25)
names(signature.centroid.list) <- c('ElasticNet_alpha0.25')

load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.Rdata')
names(signature_alternative.centroid.list) <- paste0('Alternative_',names(signature_alternative.centroid.list))

signature.centroid.list <- c(signature.centroid.list, signature_alternative.centroid.list)
rm(signature_alternative.centroid.list)

# Load drug sensitivity data
drugs <- read.csv('CCLE_breast_PRISM_drugSensitivity.csv', row.names = 1)
drugs <- drugs[,c(which(grepl(pattern = 'OLAPARIB', names(drugs))),
                  which(grepl(pattern = 'TALAZOPARIB', names(drugs))),
                  which(grepl(pattern = 'NIRAPARIB', names(drugs))),
                  which(grepl(pattern = 'RUCAPARIB', names(drugs))))]
names(drugs) <- sapply(names(drugs), function(x) strsplit(x,split='[..]')[[1]][1])
drugs$CellLine <- rownames(drugs)

# Initialise results matrix
res_olaparib <- data.frame(Signature = names(signature.centroid.list), Drug = 'olaparib', cor = NA, pVal = NA)
res_talazoparib <- data.frame(Signature = names(signature.centroid.list), Drug = 'talazoparib', cor = NA, pVal = NA)
res_niraparib <- data.frame(Signature = names(signature.centroid.list), Drug = 'niraparib', cor = NA, pVal = NA)
res_rucaparib <- data.frame(Signature = names(signature.centroid.list), Drug = 'rucaparib', cor = NA, pVal = NA)

# Correlate each signature with sensivitiy to each PARP inhibitor
for (i in 1:length(signature.centroid.list)) {
  
  print(names(signature.centroid.list)[i])
  
  sig.i <- signature.centroid.list[[i]] # Extract ith signature
  genes.intersect <- intersect(rownames(sig.i), colnames(expr))
  
  # if (i == 12) genes.intersect <- genes.intersect[genes.intersect != 'FAM170B']
  
  sig.i <- sig.i[genes.intersect, ]
  expr.i <- expr[, genes.intersect]
  
  # Calculate HRD scores
  expr.i.hrdScores <- data.frame(
    CellLine = rownames(expr.i),
    HRD = apply(expr.i, 1, function(x) cor(x,sig.i$HRD)),
    HR_proficient = apply(expr.i, 1, function(x) cor(x,sig.i$HR_proficient))
  )
  expr.i.hrdScores$HRD_score <- expr.i.hrdScores$HRD - expr.i.hrdScores$HR_proficient
  
  # Merge with drug data
  df.i <- merge(x = expr.i.hrdScores, y = drugs)
  
  # Save correlations in relevant results tables
  cor.ol <- cor.test(df.i$HRD_score, df.i$OLAPARIB)
  res_olaparib$cor[i] <- cor.ol$estimate
  res_olaparib$pVal[i] <- cor.ol$p.value
  
  cor.tal <- cor.test(df.i$HRD_score, df.i$TALAZOPARIB)
  res_talazoparib$cor[i] <- cor.tal$estimate
  res_talazoparib$pVal[i] <- cor.tal$p.value
  
  cor.nir <- cor.test(df.i$HRD_score, df.i$NIRAPARIB)
  res_niraparib$cor[i] <- cor.nir$estimate
  res_niraparib$pVal[i] <- cor.nir$p.value
  
  cor.ruc <- cor.test(df.i$HRD_score, df.i$RUCAPARIB)
  res_rucaparib$cor[i] <- cor.ruc$estimate
  res_rucaparib$pVal[i] <- cor.ruc$p.value
  
}

# Organise data table and plot relative significance values
library(data.table)
res <- rbindlist(list(res_olaparib,res_talazoparib,res_niraparib,res_rucaparib))
res$logP <- -log10(res$pVal)
res$group <- sapply(res$Signature, function(x) strsplit(x,split='_')[[1]][1])

library(ggplot2)
# g_pval <- ggplot(res, aes(x = Signature, y = logP, fill = group)) +
#   geom_bar(stat = 'identity') + geom_hline(yintercept = -log10(0.05), col = 'red', linetype = 'dashed') +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = 'top', legend.title = element_blank()) +
#   facet_wrap(~Drug, scales = 'free', nrow = 1)
# ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_CCLEsignificance.pdf',
#        plot = g_pval, height = 4, width = 7)

res$Method <- sapply(res$Signature, function(x) strsplit(x,split='_')[[1]][2])
res$Method[res$Method == 'alpha0.25'] <- 'ElasticNet'
res$Method <- factor(res$Method,
                     levels = c('ElasticNet','Severson','PARPi7','CIN70','Peng'))

g_pval <- ggplot(res, aes(x = Method, y = logP, fill = group)) +
  geom_bar(stat = 'identity') + geom_hline(yintercept = -log10(0.05), col = 'red', linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'top', legend.title = element_blank(), axis.title.x = element_blank()) +
  facet_wrap(~Drug, scales = 'free', nrow = 1) +
  scale_fill_brewer(palette = 'Pastel1')
ggsave(filename = '~/Projects/Thesis/Chapter 5/CCLE_methodComparisons.pdf',
       plot = g_pval, height = 4, width = 8)

# Plot correlate of ElasticNet_alpha0.25 score against PARPi response

sig.interest <- signature.centroid.list$ElasticNet_alpha0.25

genes.intersect <- intersect(rownames(sig.interest), colnames(expr))
sig.interest <- sig.interest[genes.intersect, ]
expr.hrd <- expr[,genes.intersect]

# Calculate HRD scores and match drug sensitivity data
df.hrd <- data.frame(
  CellLine = rownames(expr.hrd),
  HRD_score = apply(expr.hrd, 1, function(x) cor(x,sig.interest$HRD)) -
    apply(expr.hrd, 1, function(x) cor(x,sig.interest$HR_proficient))
)
df.hrd <- merge(x = df.hrd, y = drugs)

# Reorder df.hrd and plot correlations

df.hrd <- df.hrd %>%
  pivot_longer(cols = -c(CellLine, HRD_score),
               names_to = 'Drug', values_to = 'PRISM')
df.hrd$Drug <- factor(df.hrd$Drug,
                      levels = c('RUCAPARIB','NIRAPARIB',
                                 'TALAZOPARIB','OLAPARIB'))
g_ccle <- ggplot(df.hrd, aes(x = HRD_score, y = PRISM)) + 
  geom_point() + geom_smooth(method = 'lm', color = 'darkred') +
  stat_cor() + theme_minimal() +
  facet_wrap(~Drug, nrow = 1)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure5/CCLE_PARPiResponse.pdf',
       plot = g_ccle, width = 8, height = 4)


```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/GSEA_enrichR.R`:

```````R
library(enrichR)

dbs <- listEnrichrDbs()
dbs <- dbs$libraryName[c(212:214,173)]

load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
enriched <- enrichr(rownames(signature.centroid.list$ElasticNet_alpha0.25),
                    dbs)

pdf('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_GSEAenrichR.pdf')
plotEnrich(enriched$KEGG_2021_Human, showTerms = 15)
dev.off()

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/TCGA_HRDscoreByERstatus.R`:

```````R
#####
## Compare HRD transcriptional scores across ER-status
#####

# Load libraries
library(pROC)
library(ggplot2)
library(ggpubr)
library(wesanderson)

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
ann_tcga <- ann_tcga[ann_tcga$ER_status %in% c('Negative','Positive'),]

# Load reference testing data
load('~/Data/TCGA/TCGA_BRCA.ExprDeconvolution_050823_p0.79_testing.Rdata')
rownames(ann_tcga) <- ann_tcga$Patient

samples.intersect <- intersect(substr(rownames(Z.tumor_testing),1,12), rownames(ann_tcga))

ann_tcga_test <- ann_tcga[samples.intersect, ]

# Apply to non-deconvoluted samples
library(TCGAbiolinks)
library(SummarizedExperiment)

setwd('~/Data/TCGA')
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  barcode = rownames(Z.tumor_testing)
)
# GDCdownload(query)
expr.test <- GDCprepare(query = query)
expr.test <- expr.test[,expr.test$sample_type == 'Primary Tumor']
expr.test <- expr.test[!duplicated(rowData(expr.test)$gene_name) &
                         !is.na(rowData(expr.test)$gene_name), ]

library(SummarizedExperiment)
expr.tumor_testing <- assay(expr.test, 'fpkm_uq_unstrand')
rownames(expr.tumor_testing) <- rowData(expr.test)$gene_name
colnames(expr.tumor_testing) <- sapply(colnames(expr.tumor_testing),
                                       function(x) substr(x,1,12))
expr.tumor_testing <- log2(expr.tumor_testing+1)

expr.tumor_testing <- t(expr.tumor_testing)
expr.tumor_testing <- expr.tumor_testing[rownames(ann_tcga_test), ]

# Load signature centroids
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
signature.of.interest <- signature.centroid.list$ElasticNet_alpha0.25

# Calculate HRD scores for testing data and match with ann_tcga_testing
hrdScore_func <- function(expr) {
  expr.hrd = expr[,rownames(signature.of.interest)]
  cor_hrd = apply(expr.hrd, 1, function(x) cor(x,signature.of.interest$HRD))
  cor_hrproficient = apply(expr.hrd, 1, function(x) cor(x,signature.of.interest$HR_proficient))
  hrdscores = cor_hrd - cor_hrproficient
  return(hrdscores)
}

res.df <- data.frame(
  Patient = rownames(expr.tumor_testing),
  HRD_score = hrdScore_func(expr.tumor_testing)
)

res.df <- merge(x = res.df, y = ann_tcga_test[,c('Patient','HRD','ER_status')])

# Calculate AUC values for positive and negative ER status
auc.erNeg <- roc(HRD ~ HRD_score, data = res.df[res.df$ER_status == 'Negative',])
auc.erPos <- roc(HRD ~ HRD_score, data = res.df[res.df$ER_status == 'Positive',])

# Plot results by ER status
g_erPos <- ggboxplot(data = res.df[res.df$ER_status == 'Positive',],
                     x = 'HRD', y = 'HRD_score', add = 'jitter', color = 'HRD')+
  scale_color_manual(values = wes_palette('GrandBudapest1')) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) +
  stat_compare_means() +
  ggtitle(paste0('ER-positive: AUC = ',round(auc.erPos$auc,2)))

g_erNeg <- ggboxplot(data = res.df[res.df$ER_status == 'Negative',],
                     x = 'HRD', y = 'HRD_score', add = 'jitter', color = 'HRD')+
  scale_color_manual(values = wes_palette('GrandBudapest1')) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) +
  stat_compare_means() +
  ggtitle(paste0('ER-negative: AUC = ',round(auc.erNeg$auc,2)))

g_join <- ggarrange(plotlist = list(g_erPos, g_erNeg))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/SupplementaryFigures/Supp_TCGA_HRDscoreByERstatus.pdf',
       plot = g_join, width = 7, height = 4)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/GSEA_pathfindR.R`:

```````R
#####
## Run Gene Set Enrichment Analysis on the 130-gene signature using the pathfindR package
#####

library(pathfindR)

load('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')
rownames(Z.tumor_training) <- substr(rownames(Z.tumor_training), 1, 12)
Z.tumor_training <- log2(Z.tumor_training + 1)

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
# ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'
ann_tcga$group[ann_tcga$HRD == 'HR-proficient'] <- 'HR-proficient'
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

samples.intersect <- intersect(rownames(Z.tumor_training), rownames(ann_tcga))

ann_tcga <- ann_tcga[samples.intersect, ]
Z.tumor_training <- Z.tumor_training[samples.intersect, ]

# Extract relevant signature
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.brca <- signature.centroid.list$ElasticNet_alpha0.25

input.x = Z.tumor_training[,rownames(centroid.brca)]

# For each gene, run an ANOVA against the four HRD/BRCA-defect groups
#   Save the p-values and adjust accordingly
df.res <- data.frame(
  Gene.symbol = colnames(input.x),
  pVal = apply(input.x, 2, function(x)
    summary(aov(x ~ ann_tcga$BRCA_status))[[1]]$`Pr(>F)`[1])
)
df.res$adj.P.Val = p.adjust(df.res$pVal, method = 'BH')
df.res <- df.res[,c(1,3)]

# Run pathfindR and save results
df.res_output <- run_pathfindR(df.res, p_val_threshold = 1,
                               gene_sets = 'GO-BP',
                               min_gset_size = 10,
                               output_dir = '~/Projects/HRD_TranscriptionalSignature/Results/pathfindR')

# Hierarchical clustering of enriched terms and extract representative clusters
res_clustered <- cluster_enriched_terms(
  df.res_output, plot_dend = FALSE, plot_clusters_graph = FALSE)
res_clustered <- res_clustered[res_clustered$Status == 'Representative', ]

pdf('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_GSEApathfindR.pdf')
enrichment_chart(res_clustered)
dev.off()

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/TCGA_testScoring.R`:

```````R
#####
## Apply final HRD transcriptional signature to TCGA-BRCA test cohort
#####

# Load libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggpubr)
library(wesanderson)
library(tidyr)
library(pROC)

# Load data
load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata') # Group Annotation
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('Patient', 'HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'

ann_tcga$HRD <- factor(ann_tcga$HRD, levels = c('HR-proficient','HRD'))
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

load('~/Data/TCGA/TCGA_BRCA.ExprDeconvolution_050823_p0.79_testing.Rdata') # reference dataset to define test cohort

# Define test annotations
samples.intersect <- intersect(substr(rownames(Z.tumor_testing),1,12), rownames(ann_tcga))
ann_tcga_test <- ann_tcga[samples.intersect, ]

# Obain non-deconvoluted test cohort via TCGAbiolinks

setwd('~/Data/TCGA')

query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  barcode = rownames(Z.tumor_testing)
)
# GDCdownload(query)
expr.test <- GDCprepare(query = query)
expr.test <- expr.test[,expr.test$sample_type == 'Primary Tumor']
expr.test <- expr.test[!duplicated(rowData(expr.test)$gene_name) &
                         !is.na(rowData(expr.test)$gene_name), ]

# Extract FPKM-normalised cohort, log2 normalise and match to annotation
expr.tumor_testing <- assay(expr.test, 'fpkm_uq_unstrand')
rownames(expr.tumor_testing) <- rowData(expr.test)$gene_name
colnames(expr.tumor_testing) <- sapply(colnames(expr.tumor_testing),
                                       function(x) substr(x,1,12))
expr.tumor_testing <- log2(expr.tumor_testing+1)

expr.tumor_testing <- t(expr.tumor_testing)
expr.tumor_testing <- expr.tumor_testing[rownames(ann_tcga_test), ]

# Load and extract the ElasticNet_alpha0.25 signature centroid
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
sig <- signature.centroid.list$ElasticNet_alpha0.25

# Calculate and plot HRD scores (including AUC values)
results_hrd <- data.frame(
  Patient = rownames(expr.tumor_testing),
  HRD_score = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$HRD)) -
    apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$HR_proficient))
)

results_hrd <- merge(x = results_hrd, y = ann_tcga_test)

roc(HRD ~ HRD_score, data = results_hrd) # prints AUC estimate

g_hrdvshrd <- ggboxplot(data = results_hrd, x = 'HRD', y = 'HRD_score', fill = 'HRD') +
  stat_compare_means() + ylab('HRD score') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_HRDvsHRD.pdf',
       plot = g_hrdvshrd, width = 4, height = 4)

group_comparisons <- list(c('HR-proficient','BRCA2'),
                          c('HR-proficient','BRCA1'),
                          c('HR-proficient','HRD_BRCA+'))
g_hrdvsbrca <- ggboxplot(data = results_hrd, x = 'group', y = 'HRD_score', fill = 'group') +
  stat_compare_means(comparisons = group_comparisons) + ylab('HRD score') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_HRDvsBRCA.pdf',
       plot = g_hrdvsbrca, width = 5, height = 4)

# Calculate and plot BRCA defect-specific scores

results_brca <- data.frame(
  Patient = rownames(expr.tumor_testing),
  BRCA1 = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$BRCA1)),
  BRCA2 = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$BRCA2)),
  HRD_BRCApos = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$HRD_BRCApos)),
  HR_proficient = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$HR_BRCA_proficient))
)

results_brca <- merge(x = results_brca, y = ann_tcga_test)

results_brca.plot <- results_brca[,c(1:5,8)] %>%
  pivot_longer(cols = -c(Patient, group), names_to = 'Signature', values_to = 'Score')
results_brca.plot$Signature <- factor(results_brca.plot$Signature,
                                      levels = c('BRCA1','BRCA2',
                                                 'HRD_BRCApos','HR_proficient'))

g_brcavsbrca <- ggboxplot(data = results_brca.plot, x = 'group', y = 'Score', fill = 'group') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  facet_wrap(~Signature, scales = 'free')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_BRCAvsBRCA.pdf',
       plot = g_brcavsbrca, width = 8, height = 5)


## Analysis of reduced signature
genes.importance <- read.csv('~/Data/imp_score_avg.csv')
thres.imp <- 0.7
genes.important <- c(genes.importance$gene_names[genes.importance$HR.proficient > thres.imp],
                     genes.importance$gene_names[genes.importance$HRD > thres.imp])

sig_redux <- sig[genes.important,]

results_HRDredux <- data.frame(
  Patient = rownames(expr.tumor_testing),
  HRD_score = apply(expr.tumor_testing[,rownames(sig_redux)], 1, function(x) cor(x,sig_redux$HRD)) -
    apply(expr.tumor_testing[,rownames(sig_redux)], 1, function(x) cor(x,sig_redux$HR_proficient))
)

results_HRDredux <- merge(x = results_HRDredux, y = ann_tcga_test)
g_redux_hrdvshrd <- ggboxplot(data = results_HRDredux, x = 'HRD', y = 'HRD_score', fill = 'HRD') +
  stat_compare_means() + ylab('HRD score') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_TCGA_HRDvsHRD_redux.pdf',
       plot = g_redux_hrdvshrd, width = 4, height = 4)

group_comparisons <- list(c('HR-proficient','BRCA2'),
                          c('HR-proficient','BRCA1'),
                          c('HR-proficient','HRD_BRCA+'))
g_redux_hrdvsbrca <- ggboxplot(data = results_HRDredux, x = 'group', y = 'HRD_score', fill = 'group') +
  stat_compare_means(comparisons = group_comparisons) + ylab('HRD score') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_TCGA_HRDvsBRCA_redux.pdf',
       plot = g_redux_hrdvsbrca, width = 5, height = 4)

# Barplot of relevant importance values
genes.imp <- genes.importance[genes.importance$gene_names %in% genes.important, ]
genes.imp <- genes.imp[,-1]

genes.imp$Enriched <- sapply(genes.imp$HRD, function(x)
  ifelse(x > 0.7, 'HRD', 'HR-proficient'))
genes.imp$HR.proficient[genes.imp$Enriched == 'HRD'] <- 0
genes.imp$HRD[genes.imp$Enriched == 'HR-proficient'] <- 0
genes.imp$enrich_score <- genes.imp$HRD - genes.imp$HR.proficient
genes.imp <- genes.imp[order(genes.imp$enrich_score), ]
genes.imp$gene_names <- factor(genes.imp$gene_names,
                               levels = genes.imp$gene_names)

g_importance <- ggplot(genes.imp, aes(x = gene_names, y = enrich_score, fill = Enriched)) +
  geom_bar(stat = 'identity') +
  theme_minimal() + theme(axis.title.y = element_blank(),
                          axis.title.x = element_blank(),
                          legend.position = 'top',
                          legend.title = element_blank(),
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure4/Gene_ImportanceRank.pdf',
       plot = g_importance, height = 3)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/SMC_validation.R`:

```````R
library(readxl)
library(ggpubr)
library(wesanderson)
library(tidyr)

load('~/Projects/HRD_MutationalSignature/Results/SMC_HRD_resultsSummary.Rdata')
data.brca <- read_excel('~/Data/SMC_BRCA/SMC_BRCA.BrcaStatus.xlsx', skip=2)

results.smc_df$sample_id <- sapply(results.smc_df$Patient, 
                                   function(x) substr(x, 15, nchar(x)))
results.smc_df <- merge(x = results.smc_df, y = data.brca[,c('sample_id', 'gene_symbol')],
                        all.x = TRUE)
results.smc_df$group <- results.smc_df$gene_symbol
results.smc_df$group[results.smc_df$HRD_prob > 0.79 & 
                       is.na(results.smc_df$gene_symbol)] <- 'HRD_BRCA+'
results.smc_df$group[is.na(results.smc_df$group)] <- 'HR-proficient'
results.smc_df$group <- factor(results.smc_df$group,
                               levels = c('HR-proficient','HRD_BRCA+',
                                          'BRCA1','BRCA2'))

# Load validation data and match with HRD classifications
smc_rnaseq <- read.delim('~/Data/SMC_BRCA/data_mrna_seq_tpm.txt', sep='\t')
smc_rnaseq <- smc_rnaseq[!duplicated(smc_rnaseq$Hugo_Symbol), ]
rownames(smc_rnaseq) <- smc_rnaseq$Hugo_Symbol
smc_rnaseq <- smc_rnaseq[,3:ncol(smc_rnaseq)]

samples.intersect <- intersect(colnames(smc_rnaseq), results.smc_df$Patient)
smc_rnaseq <- smc_rnaseq[,samples.intersect]
results.smc_df <- results.smc_df[match(samples.intersect, results.smc_df$Patient), ]

# Load signature centroids
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
sig.interest <- signature.centroid.list$ElasticNet_alpha0.25

# Calculate HRD and BRCA-defect scores
genes.intersect <- intersect(rownames(sig.interest), rownames(smc_rnaseq))
sig.interest <- sig.interest[genes.intersect, ]
smc_rnaseq_hrd <- smc_rnaseq[genes.intersect, ]

df.hrd <- data.frame(
  Patient = colnames(smc_rnaseq_hrd),
  HRD_score = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$HRD)) -
    apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$HR_proficient)),
  BRCA1 = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$BRCA1)),
  BRCA2 = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$BRCA2)),
  HRD_BRCApos = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$HRD_BRCApos)),
  HR_BRCA_proficient = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$HR_BRCA_proficient))
)

# Match with HRD classification
df.hrd <- merge(x = df.hrd, y = results.smc_df[,c(2,6,8)])
g_HRDbyHRD <- ggboxplot(df.hrd, x = 'HRD', y = 'HRD_score', fill = 'HRD') +
  stat_compare_means() + scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_HRDvsHRD.pdf',
       plot = g_HRDbyHRD, width = 4, height = 4)

brca_comparisons <- list(c('HR-proficient','HRD_BRCA+'),
                         c('HR-proficient','BRCA1'),
                         c('HR-proficient','BRCA2'))
g_HRDvsBRCA <- ggboxplot(df.hrd, x = 'group', y = 'HRD_score', fill = 'group') +
  stat_compare_means(comparisons = brca_comparisons) + 
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_HRDvsBRCA.pdf',
       plot = g_HRDvsBRCA, width = 6, height = 4)

# Match with BRCA classifications
df.hrd_plot <- df.hrd[,c(1,3:6,8)] %>%
  pivot_longer(cols = -c(Patient,group), names_to = 'Signature', values_to = 'Correlation')
g_BRCAvsBRCA <- ggboxplot(df.hrd_plot, x = 'group', y = 'Correlation', fill = 'group') +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~Signature, scales = 'free')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_BRCAvsBRCA.pdf',
       plot = g_BRCAvsBRCA, width = 8, height = 5)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/SMC_signatureComparisons.R`:

```````R
library(pROC)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(readxl)

load('~/Projects/HRD_MutationalSignature/Results/SMC_HRD_resultsSummary.Rdata')
data.brca <- read_excel('~/Data/SMC_BRCA/SMC_BRCA.BrcaStatus.xlsx', skip=2)

results.smc_df$sample_id <- sapply(results.smc_df$Patient, 
                                   function(x) substr(x, 15, nchar(x)))
results.smc_df <- merge(x = results.smc_df, y = data.brca[,c('sample_id', 'gene_symbol')],
                        all.x = TRUE)
results.smc_df$group <- results.smc_df$gene_symbol
results.smc_df$group[results.smc_df$HRD_prob > 0.79 & 
                       is.na(results.smc_df$gene_symbol)] <- 'HRD_BRCA+'
results.smc_df$group[is.na(results.smc_df$group)] <- 'HR-proficient'
results.smc_df$group <- factor(results.smc_df$group,
                               levels = c('HR-proficient','HRD_BRCA+',
                                          'BRCA1','BRCA2'))

# Load validation data and match with HRD classifications
smc_rnaseq <- read.delim('~/Data/SMC_BRCA/data_mrna_seq_tpm.txt', sep='\t')
smc_rnaseq <- smc_rnaseq[!duplicated(smc_rnaseq$Hugo_Symbol), ]
rownames(smc_rnaseq) <- smc_rnaseq$Hugo_Symbol
smc_rnaseq <- smc_rnaseq[,3:ncol(smc_rnaseq)]

samples.intersect <- intersect(colnames(smc_rnaseq), results.smc_df$Patient)
smc_rnaseq <- smc_rnaseq[,samples.intersect]
smc_rnaseq <- as.data.frame(t(smc_rnaseq))
results.smc_df <- results.smc_df[match(samples.intersect, results.smc_df$Patient), ]

load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.Rdata')
names(signature_alternative.centroid.list) <- paste0('Alternative_',names(signature_alternative.centroid.list))

signature.centroid.list <- c(list(signature.centroid.list[[1]]), signature_alternative.centroid.list)
rm(signature_alternative.centroid.list)
names(signature.centroid.list)[1] <- 'ElasticNet'

auc.df <- data.frame(
  Model = names(signature.centroid.list),
  BRCA1 = NA, BRCA2 = NA, HRD_BRCApos = NA, HR_BRCA_proficient = NA, 
  HRD = NA
)

for (i in 1:length(signature.centroid.list)) {
  
  print(names(signature.centroid.list)[i])
  
  sig.i <- signature.centroid.list[[i]]
  
  genes.include <- intersect(rownames(sig.i), colnames(smc_rnaseq))
  sig.i <- sig.i[genes.include, ]
  
  # Start with BRCA
  sig.group <- strsplit(names(signature.centroid.list)[i],split = '_')[[1]][1]
  
  results_brca.testing <- data.frame(
    group = results.smc_df$group,
    BRCA1 = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$BRCA1)),
    BRCA2 = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$BRCA2)),
    HRD_BRCApos = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$HRD_BRCApos)),
    HR_proficient = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$HR_BRCA_proficient))
  )
  auc.df$BRCA1[i] <- roc(group == 'BRCA1' ~ BRCA1, data = results_brca.testing)$auc
  auc.df$BRCA2[i] <- roc(group == 'BRCA2' ~ BRCA2, data = results_brca.testing)$auc
  auc.df$HRD_BRCApos[i] <- roc(group == 'HRD_BRCA+' ~ HRD_BRCApos, data = results_brca.testing)$auc
  auc.df$HR_BRCA_proficient[i] <- roc(group == 'HR-proficient' ~ HR_proficient, data = results_brca.testing)$auc
  
  # Then HRD
  results_hrd.testing <- data.frame(
    HRD_status = results.smc_df$HRD,
    HRD = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$HRD)),
    HR_proficient = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$HR_proficient))
  )
  results_hrd.testing$HRD_score <- results_hrd.testing$HRD - results_hrd.testing$HR_proficient
  auc.df$HRD[i] <- roc(HRD_status == 'HRD' ~ HRD_score, data = results_hrd.testing)$auc
  
}

gene.markers <- c('BRCA1','BRCA2','POLQ','PARP1')
auc_gene.df <- data.frame(
  Model = paste0('Gene_',gene.markers),
  BRCA1 = NA, BRCA2 = NA, HRD_BRCApos = NA, HR_BRCA_proficient = NA, 
  HRD = NA
)

for (i in 1:length(gene.markers)) {
  
  df.i <- results.smc_df
  df.i$expr_i <- smc_rnaseq[,gene.markers[i]]
  
  auc_gene.df$BRCA1[i] <- roc(group  == 'BRCA1' ~ expr_i, data = df.i)$auc
  auc_gene.df$BRCA2[i] <- roc(group  == 'BRCA2' ~ expr_i, data = df.i)$auc
  auc_gene.df$HRD_BRCApos[i] <- roc(group  == 'HRD_BRCA+' ~ expr_i, data = df.i)$auc
  auc_gene.df$HR_BRCA_proficient[i] <- roc(group  == 'HR-proficient' ~ expr_i, data = df.i)$auc
  
  auc_gene.df$HRD[i] <- roc(HRD == 'HRD' ~ expr_i, data = df.i)$auc
  
}

auc.df <- rbind(auc.df, auc_gene.df)

auc.df <- auc.df[order(auc.df$HRD, decreasing = TRUE), ]
auc.df$Model <- factor(auc.df$Model, levels = auc.df$Model)
# auc.df$group <- sapply(auc.df$Model, function(x) strsplit(x,split='_')[[1]][1])
ggplot(auc.df, aes(x = Model, y = HRD, fill = Model)) + geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Spectral') +
  coord_cartesian(ylim = c(0.45,1.05)) +
  theme(axis.text.x = element_blank())

auc.df_plot <- auc.df[,-6] %>%
  pivot_longer(cols = -c(Model), names_to = 'Signature', values_to = 'AUC')
ggplot(auc.df_plot, aes(x = Model, y = AUC, fill = Model)) + geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Spectral') +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(0.45,1.05)) +
  facet_wrap(~Signature, scales = 'free', nrow = 1)

# Plot for figures
auc.df <- auc.df[order(auc.df$HRD, decreasing = TRUE), ]
auc.df$Model <- factor(auc.df$Model, levels = auc.df$Model)

g1 <- ggplot(auc.df, aes(x = Model, y = HRD, fill = Model)) +
  geom_bar(stat = 'identity') + ylim(c(0,1)) + ylab('AUC') +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))

auc.dfb <- auc.df[,1:5] %>%
  pivot_longer(cols = -Model, names_to = 'Signature', values_to = 'AUC')
auc.dfb$Signature <- factor(auc.dfb$Signature,
                            levels = c('BRCA1','BRCA2',
                                       'HRD_BRCApos', 'HR_BRCA_proficient'))
g2 <- ggplot(auc.dfb, aes(x = Model, y = AUC, fill = Model)) +
  geom_bar(stat = 'identity') + ylim(c(0,1)) + 
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18)) +
  facet_wrap(~Signature, nrow = 1)

g_join <- ggarrange(plotlist = list(g1,g2))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_AUCresults.pdf',
       plot = g_join, width = 16, height = 6)

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/TCGA_trainHeatmap.R`:

```````R
#####
## Plot Heatmap of TCGA-BRCA training cohort for chosen signature
#####

# Load libraries
library(wesanderson)
library(pheatmap)

# Load and process data
load('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')

samples <- sapply(rownames(Z.tumor_training), function(x) substr(x, 1, 12))
Z.tumor_training <- log2(Z.tumor_training + 1)
Z.tumor_training <- apply(Z.tumor_training, 2, scale)
rownames(Z.tumor_training) <- samples

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata') # Group Annotation
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('BRCA_status', 'HRD')]
ann_tcga_train <- ann_tcga[rownames(Z.tumor_training), ]

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
sig <- signature.centroid.list$ElasticNet_alpha0.25

Z.tumor_training <- Z.tumor_training[,rownames(sig)]

ann_cols <- list(
  BRCA_status = c('BRCA1' = 'blue', 'BRCA2' = 'red', 'none' = 'white'),
  HRD = c('HRD' = wes_palette('GrandBudapest1')[2], 'HR-proficient' = wes_palette('GrandBudapest1')[1])
)

# Set colour range
paletteLength <- 100
myColor <- colorRampPalette(c('navy', 'darkblue', 'white', 'red', 'darkred'))(paletteLength)
myBreaks <- c(seq(min(Z.tumor_training), -3, length.out = ceiling(paletteLength/4)+1),
              seq(-3, 0, length.out = floor(paletteLength/4))[-1],
              seq(0, 3, length.out = floor(paletteLength/4))[-1],
              seq(3, max(Z.tumor_training), length.out = floor(paletteLength/4))[-1])

pheatmap(t(Z.tumor_training), show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = ann_tcga_train, annotation_colors = ann_cols,
         color = myColor, breaks = myBreaks, clustering_method = 'average'
         , filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_Heatmap.pdf'
         )

```````

`/Users/leojo/Developer/git/alexandrov_sd/MultiscaleHRD/Scripts/TranscriptionalSignature/SignatureValidation/TCGA_testScoring_signatureComparisons.R`:

```````R
#####
## Compare transcriptional signatures against TCGA-BRCA training cohort
#####

# Load libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggpubr)
library(wesanderson)
library(tidyr)
library(pROC)

# Load data
load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata') # Group Annotation
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('Patient', 'HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'

ann_tcga$HRD <- factor(ann_tcga$HRD, levels = c('HR-proficient','HRD'))
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

load('~/Data/TCGA/TCGA_BRCA.ExprDeconvolution_050823_p0.79_testing.Rdata') # reference dataset to define test cohort

# Define test annotations
samples.intersect <- intersect(substr(rownames(Z.tumor_testing),1,12), rownames(ann_tcga))
ann_tcga_test <- ann_tcga[samples.intersect, ]

# Obain non-deconvoluted test cohort via TCGAbiolinks

setwd('~/Data/TCGA')

query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  barcode = rownames(Z.tumor_testing)
)
# GDCdownload(query)
expr.test <- GDCprepare(query = query)
expr.test <- expr.test[,expr.test$sample_type == 'Primary Tumor']
expr.test <- expr.test[!duplicated(rowData(expr.test)$gene_name) &
                         !is.na(rowData(expr.test)$gene_name), ]

# Extract FPKM-normalised cohort, log2 normalise and match to annotation
expr.tumor_testing <- assay(expr.test, 'fpkm_uq_unstrand')
rownames(expr.tumor_testing) <- rowData(expr.test)$gene_name
colnames(expr.tumor_testing) <- sapply(colnames(expr.tumor_testing),
                                       function(x) substr(x,1,12))
expr.tumor_testing <- log2(expr.tumor_testing+1)

expr.tumor_testing <- t(expr.tumor_testing)
expr.tumor_testing <- expr.tumor_testing[rownames(ann_tcga_test), ]

# Load complete set of signatures
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.Rdata')
names(signature_alternative.centroid.list) <- paste0('Alternative_',names(signature_alternative.centroid.list))

signature.centroid.list <- c(signature.centroid.list, signature_alternative.centroid.list)
rm(signature_alternative.centroid.list)

# Initialise data frame to collate AUC values
auc.df <- data.frame(
  Model = names(signature.centroid.list),
  BRCA1 = NA, BRCA2 = NA, HRD_BRCApos = NA, HR_BRCA_proficient = NA, 
  HRD = NA
)

# For each model:
#   Calculate HRD and BRCA-specific signature scores
#   Calculate the relevant AUC for each statistic

for (i in 1:length(signature.centroid.list)) {
  
  print(names(signature.centroid.list)[i])
  
  sig.i <- signature.centroid.list[[i]]
  genes.include <- rownames(sig.i)
  
  # Start with BRCA
  sig.group <- strsplit(names(signature.centroid.list)[i],split = '_')[[1]][1]
  
  results_brca.testing <- data.frame(
    group = ann_tcga_test$group,
    BRCA1 = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$BRCA1)),
    BRCA2 = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$BRCA2)),
    HRD_BRCApos = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$HRD_BRCApos)),
    HR_proficient = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$HR_BRCA_proficient))
  )
  
  auc.df$BRCA1[i] <- roc(group == 'BRCA1' ~ BRCA1, data = results_brca.testing)$auc
  auc.df$BRCA2[i] <- roc(group == 'BRCA2' ~ BRCA2, data = results_brca.testing)$auc
  auc.df$HRD_BRCApos[i] <- roc(group == 'HRD_BRCA+' ~ HRD_BRCApos, data = results_brca.testing)$auc
  auc.df$HR_BRCA_proficient[i] <- roc(group == 'HR-proficient' ~ HR_proficient, data = results_brca.testing)$auc
  
  # Then HRD
  results_hrd.testing <- data.frame(
    HRD_status = ann_tcga_test$HRD,
    HRD = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$HRD)),
    HR_proficient = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$HR_proficient))
  )
  results_hrd.testing$HRD_score <- results_hrd.testing$HRD - results_hrd.testing$HR_proficient
  auc.df$HRD[i] <- roc(HRD_status == 'HRD' ~ HRD_score, data = results_hrd.testing)$auc
  
}

# Repeat for relevant gene markers and add to AUC results

gene.markers <- c('BRCA1','BRCA2','POLQ','PARP1')
auc_gene.df <- data.frame(
  Model = paste0('Gene_',gene.markers),
  BRCA1 = NA, BRCA2 = NA, HRD_BRCApos = NA, HR_BRCA_proficient = NA, 
  HRD = NA
)

for (i in 1:length(gene.markers)) {
  
  df.i <- ann_tcga_test
  df.i$expr_i <- expr.tumor_testing[,gene.markers[i]]
  
  auc_gene.df$BRCA1[i] <- roc(group  == 'BRCA1' ~ expr_i, data = df.i)$auc
  auc_gene.df$BRCA2[i] <- roc(group  == 'BRCA2' ~ expr_i, data = df.i)$auc
  auc_gene.df$HRD_BRCApos[i] <- roc(group  == 'HRD_BRCA+' ~ expr_i, data = df.i)$auc
  auc_gene.df$HR_BRCA_proficient[i] <- roc(group  == 'HR-proficient' ~ expr_i, data = df.i)$auc
  
  auc_gene.df$HRD[i] <- roc(HRD == 'HRD' ~ expr_i, data = df.i)$auc
  
}

auc.df <- rbind(auc.df, auc_gene.df)

# For plotting, remove other regression signatures, and plot AUC in descending order of HRD prediction
auc.df_plot <- auc.df[c(1,5:nrow(auc.df)), ]
auc.df_plot$Signature <- sapply(auc.df_plot$Model, function(x) strsplit(x,split='_')[[1]][2])
auc.df_plot$Signature[1] <- 'ElasticNet'
auc.df_plot <- auc.df_plot[order(auc.df_plot$HRD, decreasing = TRUE), ]
auc.df_plot$Signature <- factor(auc.df_plot$Signature, levels = auc.df_plot$Signature)

# Plot HRD AUCs
g_auc_hrd <- ggplot(auc.df_plot, aes(x = Signature, y = HRD, fill = Signature)) + 
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(0.5, 1)) + ylab('AUC') +
  scale_fill_brewer(palette = 'Spectral')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_AUCbyHRD.pdf',
       plot = g_auc_hrd, height = 4, width = 7)

# Plot BRCA-specific AUCs
auc.df_plot2 <- auc.df_plot[,c(7,2:5)] %>%
  pivot_longer(cols = -Signature, names_to = 'group', values_to = 'AUC')
auc.df_plot2$group <- factor(auc.df_plot2$group,
                             levels = c('BRCA1','BRCA2',
                                        'HRD_BRCApos','HR_BRCA_proficient'))
g_auc_brca <- ggplot(auc.df_plot2, aes(x = Signature, y = AUC, fill = Signature)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_fill_brewer(palette = 'Spectral') +
  facet_wrap(~ group, nrow = 1)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_TCGA_AUCbyBRCA.pdf',
       plot = g_auc_brca, width = 7, height = 4)

```````
