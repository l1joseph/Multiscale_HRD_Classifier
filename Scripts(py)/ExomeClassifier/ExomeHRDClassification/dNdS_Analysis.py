import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import pyreadr  # for reading R data files
import mygene  # for gene ID conversion


# Load data and dN/dS references
ref_cds = pyreadr.read_r("~/Data/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda")["RefCDS"]

# Load TCGA data
# Note: We might need to implement a  function to replicate GDCquery, GDCdownload, and GDCprepare
# I'm pretty sure data is already downloaded and prepared in git repo though.
mutations = pd.read_csv("~/Data/TCGA/mutations.csv")

# Process mutations data
data = mutations[
    [
        "Tumor_Sample_Barcode",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
    ]
]
data.columns = ["sampleID", "chr", "pos", "ref", "mut"]

# Load HRD/BRCA groups and arrange data accordingly
ann_tcga = pd.read_pickle("Results/TCGA_HRDclassification_BRCAannotation.pkl")
ann_tcga = ann_tcga.dropna(subset=["BRCA_status"])
ann_tcga = ann_tcga[~ann_tcga["BRCA_status"].isin(["RAD51C", "PALB2"])]

data = pd.merge(
    data,
    ann_tcga[["Patient", "HRD", "BRCA_status"]],
    left_on="sampleID",
    right_on="Patient",
)

# HRD vs HR-proficient
data_HRD = data[data["HRD"] == "HRD"].iloc[:, :5]
data_HRproficient = data[data["HRD"] == "HR-proficient"].iloc[:, :5]

# BRCA-defect categories
data_BRCA1 = data[data["BRCA_status"] == "BRCA1"].iloc[:, :5]
data_BRCA2 = data[data["BRCA_status"] == "BRCA2"].iloc[:, :5]
data_RAD51C = data[data["BRCA_status"] == "RAD51C"].iloc[:, :5]
data_HRDBRCApos = data[(data["HRD"] == "HRD") & (data["BRCA_status"] == "none")].iloc[
    :, :5
]

# Run dN/dS analysis on each group separately
# Note: dndscv doesn't have a direct Python equivalent. We might need to implement this functionality
# or use a different tool. I used a placeholder function.


def dndscv(data, refdb):
    # Placeholder function. Implement dN/dS calculation here.
    return pd.DataFrame({"gene_name": ["GENE1", "GENE2"], "qind_cv": [0.01, 0.05]})


dndsoutBRCA1 = dndscv(data_BRCA1, ref_cds)
dndsoutBRCA2 = dndscv(data_BRCA2, ref_cds)
dndsoutRAD51C = dndscv(data_RAD51C, ref_cds)
dndsoutHRDBRCApos = dndscv(data_HRDBRCApos, ref_cds)
dndsoutHRD = dndscv(data_HRD, ref_cds)
dndsoutHRprof = dndscv(data_HRproficient, ref_cds)

# Find genes under positive selection in at least one group
sig_genes = set(dndsoutHRD[dndsoutHRD["qind_cv"] < 0.1]["gene_name"]) | set(
    dndsoutHRprof[dndsoutHRprof["qind_cv"] < 0.1]["gene_name"]
)

sel_cvSigGenes = pd.merge(
    dndsoutHRD[dndsoutHRD["gene_name"].isin(sig_genes)][
        ["gene_name", "wind_cv", "qind_cv"]
    ],
    dndsoutHRprof[dndsoutHRprof["gene_name"].isin(sig_genes)][
        ["gene_name", "wind_cv", "qind_cv"]
    ],
    on="gene_name",
    suffixes=("_HRD", "_HRprof"),
)

# Prepare labels for plotting
sel_cvSigGenes["Significant"] = "Both"
sel_cvSigGenes.loc[sel_cvSigGenes["qind_cv_HRD"] > 0.1, "Significant"] = (
    "HR-proficient ONLY"
)
sel_cvSigGenes.loc[sel_cvSigGenes["qind_cv_HRprof"] > 0.1, "Significant"] = "HRD ONLY"

sel_cvSigGenes["dNdS_HRD"] = np.log2(sel_cvSigGenes["wind_cv_HRD"] + 1)
sel_cvSigGenes["dNdS_HRprof"] = np.log2(sel_cvSigGenes["wind_cv_HRprof"] + 1)

# Plot dN/dS of indsense variants across all groups
sel_cvFull = pd.concat(
    [
        dndsoutHRprof[["gene_name", "wmis_cv", "qmis_cv"]],
        dndsoutHRD[["gene_name", "wmis_cv", "qmis_cv"]],
        dndsoutBRCA1[["gene_name", "wmis_cv", "qmis_cv"]],
        dndsoutBRCA2[["gene_name", "wmis_cv", "qmis_cv"]],
        dndsoutRAD51C[["gene_name", "wmis_cv", "qmis_cv"]],
        dndsoutHRDBRCApos[["gene_name", "wmis_cv", "qmis_cv"]],
    ]
)
sel_cvFull["group"] = np.repeat(
    ["HR-proficient", "HRD full", "BRCA1", "BRCA2", "RAD51C", "HRD BRCA+"],
    len(sel_cvFull) // 6,
)

sel_cvFull = sel_cvFull[sel_cvFull["qmis_cv"] < 0.05]

sel_cvFull["log_dNdS"] = np.log2(sel_cvFull["wmis_cv"])
sel_cvFull["log_FDR"] = -np.log10(sel_cvFull["qmis_cv"])
sel_cvFull.loc[sel_cvFull["log_FDR"] == np.inf, "log_FDR"] = 9
sel_cvFull.loc[sel_cvFull["log_FDR"] >= 9, "log_FDR"] = 9

sel_cvFull["group"] = pd.Categorical(
    sel_cvFull["group"],
    categories=["HR-proficient", "HRD full", "BRCA1", "BRCA2", "RAD51C", "HRD BRCA+"],
)

# Plot
plt.figure(figsize=(4, 7))
sns.scatterplot(
    data=sel_cvFull,
    x="gene_name",
    y="group",
    size="log_dNdS",
    hue="log_FDR",
    palette="Blues",
    sizes=(20, 200),
)
plt.axhline(y=2.5, color="red", linestyle="--")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.xlabel("Gene")
plt.ylabel("Group")
plt.title("dN/dS Analysis")
plt.tight_layout()
plt.savefig("Figures/Figure2/HRD_dNdSballoonPlot.pdf")
plt.close()


# notes

# need to implement dndscv function
# same for GDCquery, GDCdownload, and GDCprepare if data isn't as we expect.
# also pyreadr is to read r files, but I'm not sure if it is needed if there is a Non-R-data equivilant to the dnds data.
# I'm also not sure if mygene is the best package for geneID conversion. I know that there is similar api called Biothings.
