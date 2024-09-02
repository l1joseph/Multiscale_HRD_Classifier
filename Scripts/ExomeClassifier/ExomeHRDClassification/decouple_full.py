import os
import pandas as pd
import numpy as np
from TCGAbiolinks import GDCquery, GDCdownload, GDCprepare
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
import seaborn as sns

# Set working directory
os.chdir("~/Data/TCGA")

# Load required R packages
pandas2ri.activate()
deseq2 = importr("DESeq2")
decoupleR = importr("decoupleR")

# Query TCGA data
query = GDCquery(
    project="TCGA-BRCA",
    data_category="Transcriptome Profiling",
    data_type="Gene Expression Quantification",
    workflow_type="STAR - Counts",
)
# GDCdownload(query)
data = GDCprepare(query=query)

# Filter data
data = data[data.sample_type == "Primary Tumor"]
data = data[~data.patient.duplicated()]

# Load HRD classification
r(
    'load("~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata")'
)
ann_tcga = r("ann_tcga")
ann_tcga["HRD"] = ["HRD" if x >= 0.79 else "HR-proficient" for x in ann_tcga.HRD_prob]
ann_tcga = ann_tcga[~ann_tcga.ER_status.isna()]

# Intersect patients
patients_intersect = set(ann_tcga.Patient).intersection(set(data.patient))
ann_tcga = ann_tcga[ann_tcga.Patient.isin(patients_intersect)]
data = data[data.patient.isin(patients_intersect)]

data["HRD_status"] = pd.Categorical(ann_tcga.HRD, categories=["HR-proficient", "HRD"])

# Construct DESeqDataSet
dds = deseq2.DESeqDataSetFromMatrix(
    countData=data.counts, colData=data.colData, design=r("~ HRD_status")
)

# Gene count filtering
dds = dds[r("rowSums(counts(dds)) > 10"),]

# Normalization
dds = deseq2.estimateSizeFactors(dds)

# Differential gene expression analysis
dds_DGE = deseq2.DESeq(dds)
dds_DGE_results = deseq2.results(dds_DGE)

dds_statVals = pd.DataFrame(
    {"EnsemblID": dds_DGE_results.names, "stat": dds_DGE_results.stat}
)
dds_statVals["ID"] = data.rowData.gene_name[
    data.rowData.index.isin(dds_statVals.EnsemblID)
]
dds_statVals = dds_statVals[~dds_statVals.ID.duplicated()]
dds_statVals = dds_statVals.reset_index(drop=True)

deg = dds_statVals.set_index("ID")["stat"].to_frame()

counts = data.assay.loc[dds_statVals.EnsemblID]
counts.index = dds_statVals.ID
counts.columns = [x[:12] for x in counts.columns]

counts_logNorm = np.log2(counts + 1)

design = ann_tcga[["Patient", "HRD"]]
design.columns = ["sample", "condition"]

# Run progeny
net = decoupleR.get_progeny(organism="human", top=100)

# Run mlm
contrast_acts = decoupleR.run_mlm(
    mat=deg, net=net, source="source", target="target", mor="weight", minsize=5
)

# Plot
g_progeny = sns.barplot(
    x="source",
    y="score",
    data=contrast_acts,
    order=contrast_acts.sort_values("score", ascending=False).source,
)
g_progeny.set_xticklabels(g_progeny.get_xticklabels(), rotation=45, ha="right")
plt.title("Enrichment in HRD Samples")
plt.tight_layout()
plt.savefig("~/Projects/HRD_MutationalSignature/Figures/Figure2/HRD_decoupleR.pdf")
plt.close()
