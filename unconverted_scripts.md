Scripts
├── ExomeClassifier
│ ├── ExomeHRDClassification
│ │ ├── SMC_HRDclassification.R (done)
│ │ ├── TCGA_geneCNAenrichment.R (done)
│ │ ├── decoupleR_full.R (done)
│ │ ├── TCGA_chromosomeArmAlterations_HRDonly.R (done)
│ │ ├── HRD_Hypoxia.R (done)
│ │ ├── TCGA_chromosomeArmAlterations.R (done)
│ │ ├── dNdS_Analysis.R (done)
│ │ ├── TCGA_HRDclassification.R (done)
│ │ ├── TCGA_HRDhallmarks.R (done)
│ │ ├── TCGA_geneCNAenrichment_HRDonly.R (done)
│ │ ├── TCGA_compareHRDclassifiers.R
│ │ └── TCGA_HRDclassification_thresholds.R
│ └── ICGC_PhenotypeClusteringAndSimulationAnalysis
│ ├── Likelihood_Cluster_Generation.R
│ ├── ICGC_deconstructSigs_genome2exome.R (done)
│ ├── ICGC_BRCA_sigProfilerContributions.R (done)
│ ├── ICGC_simulations_indelProportions.R (done)
│ ├── ICGC_BRCA_UKandEU_MatrixGeneration.R
│ ├── ICGC_simulations_indelLikelihoodAlterations.R (done)
│ ├── ICGC_simulations_clusterReassignment.R (done)
│ ├── ICGC_PhenotypeClusterDevelopment_normalised.R (done)
│ └── ICGC_indelCounts_ExomeVsGenome.R
├── TCGA_HRDclassificationHeatmap.jpg
├── README.md
└── TranscriptionalSignature
├── TranscriptionalSignatureDevelopment
│ ├── RNAseqTransformationScript.R
│ ├── MultinomialWEN_alpha0.25_1to250.R
│ ├── ReduxSig_adjacencyMatrixFormation.R
│ ├── BayesPrism_TumorPurityComparison.R
│ ├── CollateAlternativeSignatures.R
│ ├── Qian2020_signatureDropOuts.R
│ ├── CentroidModelFormation.R
│ ├── MultinomialElasticNet_alpha0.25_1to100.R
│ └── TCGA_BRCA.RNASeq_prep.R
├── SingleCellAnalysis
│ ├── Bassez2021_preprocessing.R
│ ├── Qian2020_jointComparison.R
│ ├── Bassez2021_analysis.R
│ ├── Bassez2021_HRDprofiling.R
│ ├── Qian2020_analysis.R
│ ├── Chung2017_analysis.R
│ ├── Qian2020_preprocessing.R
│ ├── Qian2020_HRDprofiling.R
│ ├── cpdb_interactionCounts.R
│ └── cpdb_individualInteractions.R
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