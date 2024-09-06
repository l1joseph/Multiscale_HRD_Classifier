import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import gseapy as gp

# Load deconvoluted training data for reference
Z_tumor_training = pd.read_pickle(
    "~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.pkl"
)
Z_tumor_training.index = Z_tumor_training.index.str[:12]
Z_tumor_training = np.log2(Z_tumor_training + 1)

# Load HRD/BRCA status annotation
ann_tcga = pd.read_pickle(
    "~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl"
)
ann_tcga = ann_tcga.set_index("Patient")
ann_tcga = ann_tcga[["HRD", "BRCA_status"]]
ann_tcga = ann_tcga[~ann_tcga["BRCA_status"].isna()]
ann_tcga = ann_tcga[~ann_tcga["BRCA_status"].isin(["PALB2", "RAD51C"])]

ann_tcga["group"] = ann_tcga["BRCA_status"]
ann_tcga.loc[
    (ann_tcga["HRD"] == "HRD") & (ann_tcga["BRCA_status"] == "none"), "group"
] = "HRD_BRCA+"
ann_tcga.loc[ann_tcga["group"] == "none", "group"] = "HR-proficient"
ann_tcga["group"] = pd.Categorical(
    ann_tcga["group"], categories=["HR-proficient", "HRD_BRCA+", "BRCA1", "BRCA2"]
)

samples_intersect = Z_tumor_training.index.intersection(ann_tcga.index)

ann_tcga = ann_tcga.loc[samples_intersect]
Z_tumor_training = Z_tumor_training.loc[samples_intersect]

# Extract relevant signature
signature_centroid_list = pd.read_pickle(
    "~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl"
)
centroid_brca = signature_centroid_list["ElasticNet_alpha0.25"]

input_x = Z_tumor_training[centroid_brca.index]

# For each gene, run an ANOVA against the four HRD/BRCA-defect groups
# Save the p-values and adjust accordingly
df_res = pd.DataFrame(
    {
        "Gene.symbol": input_x.columns,
        "pVal": [
            stats.f_oneway(
                *[
                    input_x[gene][ann_tcga["group"] == group]
                    for group in ann_tcga["group"].cat.categories
                ]
            )[1]
            for gene in input_x.columns
        ],
    }
)
df_res["adj.P.Val"] = multipletests(df_res["pVal"], method="fdr_bh")[1]
df_res = df_res[["Gene.symbol", "adj.P.Val"]]

# Run GSEA using gseapy
gmt = gp.get_library("GO_Biological_Process_2021")
res = gp.enrichr(
    gene_list=df_res["Gene.symbol"].tolist(),
    gene_sets=gmt,
    cutoff=1,  # set cutoff to 1 to include all genes
    outdir="~/Projects/HRD_TranscriptionalSignature/Results/gseapy_output",
    no_plot=True,
)

# Filter and sort results
res_df = res.results
res_df = res_df.sort_values("Adjusted P-value")

# Plot results
plt.figure(figsize=(12, 10))
sns.barplot(x="Adjusted P-value", y="Term", data=res_df.head(20), palette="YlOrRd")
plt.title("GO Biological Process Enrichment")
plt.tight_layout()
plt.savefig("~/Projects/HRD_TranscriptionalSignature/Figures/Supp_GSEApathfindR.pdf")
plt.close()

# Save full results
res_df.to_csv(
    "~/Projects/HRD_TranscriptionalSignature/Results/gseapy_results.csv", index=False
)


# notes:

# Instead of the R package pathfindR, I'm using the Python package gseapy, which provides similar functionality for GSEA.
