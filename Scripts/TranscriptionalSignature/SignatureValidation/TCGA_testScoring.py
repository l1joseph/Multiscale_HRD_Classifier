import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import roc_curve, auc
import anndata as ad

# Load data
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
ann_tcga["HRD"] = pd.Categorical(ann_tcga["HRD"], categories=["HR-proficient", "HRD"])

Z_tumor_testing = pd.read_pickle(
    "~/Data/TCGA/TCGA_BRCA.ExprDeconvolution_050823_p0.79_testing.pkl"
)

samples_intersect = [
    sample[:12] for sample in Z_tumor_testing.index if sample[:12] in ann_tcga.index
]
ann_tcga_test = ann_tcga.loc[samples_intersect]

# Load non-deconvoluted samples
expr_test = ad.read_h5ad("~/Data/TCGA/TCGA_BRCA.counts.SE_050823_testing.h5ad")
expr_test = expr_test[expr_test.obs["sample_type"] == "Primary Tumor"]
expr_test = expr_test[
    :, ~expr_test.var_names.duplicated() & expr_test.var["gene_name"].notna()
]

expr_tumor_testing = pd.DataFrame(
    expr_test.layers["fpkm_uq_unstrand"].T,
    index=expr_test.var["gene_name"],
    columns=[sample[:12] for sample in expr_test.obs_names],
)
expr_tumor_testing = np.log2(expr_tumor_testing + 1)
expr_tumor_testing = expr_tumor_testing.loc[:, ann_tcga_test.index]

# Load signature centroids
signature_centroid_list = pd.read_pickle(
    "~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl"
)
sig = signature_centroid_list["ElasticNet_alpha0.25"]

# Calculate and plot HRD scores (including AUC values)
results_hrd = pd.DataFrame(
    {
        "Patient": expr_tumor_testing.columns,
        "HRD_score": expr_tumor_testing.loc[sig.index].apply(
            lambda x: np.corrcoef(x, sig["HRD"])[0, 1]
            - np.corrcoef(x, sig["HR_proficient"])[0, 1]
        ),
    }
)

results_hrd = results_hrd.merge(ann_tcga_test, left_on="Patient", right_index=True)

fpr, tpr, _ = roc_curve(results_hrd["HRD"] == "HRD", results_hrd["HRD_score"])
roc_auc = auc(fpr, tpr)
print(f"AUC: {roc_auc:.3f}")

plt.figure(figsize=(4, 4))
sns.boxplot(data=results_hrd, x="HRD", y="HRD_score", palette="Set3")
sns.swarmplot(data=results_hrd, x="HRD", y="HRD_score", color=".2", size=3)
plt.title("HRD vs HRD score")
plt.xlabel("")
plt.ylabel("HRD score")
statistic, pvalue = stats.ttest_ind(
    results_hrd[results_hrd["HRD"] == "HRD"]["HRD_score"],
    results_hrd[results_hrd["HRD"] == "HR-proficient"]["HRD_score"],
)
plt.text(0.5, plt.ylim()[1], f"p-value: {pvalue:.4e}", horizontalalignment="center")
plt.tight_layout()
plt.savefig("~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_HRDvsHRD.pdf")
plt.close()

group_comparisons = [
    ("HR-proficient", "BRCA2"),
    ("HR-proficient", "BRCA1"),
    ("HR-proficient", "HRD_BRCA+"),
]
plt.figure(figsize=(5, 4))
sns.boxplot(data=results_hrd, x="group", y="HRD_score", palette="Set3")
sns.swarmplot(data=results_hrd, x="group", y="HRD_score", color=".2", size=3)
plt.title("HRD score vs BRCA status")
plt.xlabel("")
plt.ylabel("HRD score")
for i, comparison in enumerate(group_comparisons):
    group1 = results_hrd[results_hrd["group"] == comparison[0]]["HRD_score"]
    group2 = results_hrd[results_hrd["group"] == comparison[1]]["HRD_score"]
    statistic, pvalue = stats.ttest_ind(group1, group2)
    plt.text(
        0.5,
        plt.ylim()[1] - (i + 1) * 0.1 * (plt.ylim()[1] - plt.ylim()[0]),
        f"{comparison[0]} vs {comparison[1]}: p={pvalue:.4e}",
        horizontalalignment="center",
    )
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_HRDvsBRCA.pdf"
)
plt.close()

# Calculate and plot BRCA defect-specific scores
results_brca = pd.DataFrame(
    {
        "Patient": expr_tumor_testing.columns,
        "BRCA1": expr_tumor_testing.loc[sig.index].apply(
            lambda x: np.corrcoef(x, sig["BRCA1"])[0, 1]
        ),
        "BRCA2": expr_tumor_testing.loc[sig.index].apply(
            lambda x: np.corrcoef(x, sig["BRCA2"])[0, 1]
        ),
        "HRD_BRCApos": expr_tumor_testing.loc[sig.index].apply(
            lambda x: np.corrcoef(x, sig["HRD_BRCApos"])[0, 1]
        ),
        "HR_proficient": expr_tumor_testing.loc[sig.index].apply(
            lambda x: np.corrcoef(x, sig["HR_BRCA_proficient"])[0, 1]
        ),
    }
)

results_brca = results_brca.merge(
    ann_tcga_test[["group"]], left_on="Patient", right_index=True
)

results_brca_plot = results_brca.melt(
    id_vars=["Patient", "group"], var_name="Signature", value_name="Score"
)
results_brca_plot["Signature"] = pd.Categorical(
    results_brca_plot["Signature"],
    categories=["BRCA1", "BRCA2", "HRD_BRCApos", "HR_proficient"],
)

plt.figure(figsize=(8, 5))
sns.boxplot(
    data=results_brca_plot, x="group", y="Score", hue="Signature", palette="Set3"
)
plt.title("BRCA-specific scores")
plt.xlabel("")
plt.legend(title="")
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_BRCAvsBRCA.pdf"
)
plt.close()

# Analysis of reduced signature
genes_importance = pd.read_csv("~/Data/imp_score_avg.csv")
thres_imp = 0.7
genes_important = (
    genes_importance[genes_importance["HR.proficient"] > thres_imp][
        "gene_names"
    ].tolist()
    + genes_importance[genes_importance["HRD"] > thres_imp]["gene_names"].tolist()
)

sig_redux = sig.loc[genes_important]

results_HRDredux = pd.DataFrame(
    {
        "Patient": expr_tumor_testing.columns,
        "HRD_score": expr_tumor_testing.loc[sig_redux.index].apply(
            lambda x: np.corrcoef(x, sig_redux["HRD"])[0, 1]
            - np.corrcoef(x, sig_redux["HR_proficient"])[0, 1]
        ),
    }
)

results_HRDredux = results_HRDredux.merge(
    ann_tcga_test, left_on="Patient", right_index=True
)

plt.figure(figsize=(4, 4))
sns.boxplot(data=results_HRDredux, x="HRD", y="HRD_score", palette="Set3")
sns.swarmplot(data=results_HRDredux, x="HRD", y="HRD_score", color=".2", size=3)
plt.title("HRD vs HRD score (reduced signature)")
plt.xlabel("")
plt.ylabel("HRD score")
statistic, pvalue = stats.ttest_ind(
    results_HRDredux[results_HRDredux["HRD"] == "HRD"]["HRD_score"],
    results_HRDredux[results_HRDredux["HRD"] == "HR-proficient"]["HRD_score"],
)
plt.text(0.5, plt.ylim()[1], f"p-value: {pvalue:.4e}", horizontalalignment="center")
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/Supp_TCGA_HRDvsHRD_redux.pdf"
)
plt.close()

plt.figure(figsize=(5, 4))
sns.boxplot(data=results_HRDredux, x="group", y="HRD_score", palette="Set3")
sns.swarmplot(data=results_HRDredux, x="group", y="HRD_score", color=".2", size=3)
plt.title("HRD score vs BRCA status (reduced signature)")
plt.xlabel("")
plt.ylabel("HRD score")
for i, comparison in enumerate(group_comparisons):
    group1 = results_HRDredux[results_HRDredux["group"] == comparison[0]]["HRD_score"]
    group2 = results_HRDredux[results_HRDredux["group"] == comparison[1]]["HRD_score"]
    statistic, pvalue = stats.ttest_ind(group1, group2)
    plt.text(
        0.5,
        plt.ylim()[1] - (i + 1) * 0.1 * (plt.ylim()[1] - plt.ylim()[0]),
        f"{comparison[0]} vs {comparison[1]}: p={pvalue:.4e}",
        horizontalalignment="center",
    )
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/Supp_TCGA_HRDvsBRCA_redux.pdf"
)
plt.close()

# Barplot of relevant importance values
genes_imp = genes_importance[
    genes_importance["gene_names"].isin(genes_important)
].copy()
genes_imp["Enriched"] = np.where(genes_imp["HRD"] > 0.7, "HRD", "HR-proficient")
genes_imp.loc[genes_imp["Enriched"] == "HRD", "HR.proficient"] = 0
genes_imp.loc[genes_imp["Enriched"] == "HR-proficient", "HRD"] = 0
genes_imp["enrich_score"] = genes_imp["HRD"] - genes_imp["HR.proficient"]
genes_imp = genes_imp.sort_values("enrich_score")
genes_imp["gene_names"] = pd.Categorical(
    genes_imp["gene_names"], categories=genes_imp["gene_names"]
)

plt.figure(figsize=(10, 3))
sns.barplot(
    data=genes_imp, x="gene_names", y="enrich_score", hue="Enriched", palette="Set3"
)
plt.title("Gene Importance Ranking")
plt.xlabel("")
plt.ylabel("Enrichment Score")
plt.xticks(rotation=90)
plt.legend(title="")
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/Figure4/Gene_ImportanceRank.pdf"
)
plt.close()
