import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import roc_curve, auc

# Load TCGA annotation data
ann_tcga = pd.read_pickle(
    "~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl"
)
ann_tcga = ann_tcga[ann_tcga["ER_status"].isin(["Negative", "Positive"])]

# Load reference testing data
Z_tumor_testing = pd.read_pickle(
    "~/Data/TCGA/TCGA_BRCA.ExprDeconvolution_050823_p0.79_testing.pkl"
)
ann_tcga = ann_tcga.set_index("Patient")

samples_intersect = [
    sample[:12] for sample in Z_tumor_testing.index if sample[:12] in ann_tcga.index
]

ann_tcga_test = ann_tcga.loc[samples_intersect]

# Load non-deconvoluted samples
import anndata as ad

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
signature_of_interest = signature_centroid_list["ElasticNet_alpha0.25"]


# Calculate HRD scores for testing data and match with ann_tcga_testing
def hrdScore_func(expr):
    expr_hrd = expr.loc[signature_of_interest.index]
    cor_hrd = expr_hrd.apply(
        lambda x: np.corrcoef(x, signature_of_interest["HRD"])[0, 1]
    )
    cor_hrproficient = expr_hrd.apply(
        lambda x: np.corrcoef(x, signature_of_interest["HR_proficient"])[0, 1]
    )
    return cor_hrd - cor_hrproficient


res_df = pd.DataFrame(
    {
        "Patient": expr_tumor_testing.columns,
        "HRD_score": hrdScore_func(expr_tumor_testing),
    }
)

res_df = res_df.merge(
    ann_tcga_test[["HRD", "ER_status"]], left_on="Patient", right_index=True
)


# Calculate AUC values for positive and negative ER status
def calculate_auc(df, er_status):
    subset = df[df["ER_status"] == er_status]
    fpr, tpr, _ = roc_curve(subset["HRD"] == "HRD", subset["HRD_score"])
    return auc(fpr, tpr)


auc_erNeg = calculate_auc(res_df, "Negative")
auc_erPos = calculate_auc(res_df, "Positive")

# Plot results by ER status
plt.figure(figsize=(7, 4))
sns.boxplot(
    data=res_df[res_df["ER_status"] == "Positive"],
    x="HRD",
    y="HRD_score",
    palette="Set3",
)
sns.stripplot(
    data=res_df[res_df["ER_status"] == "Positive"],
    x="HRD",
    y="HRD_score",
    color=".3",
    size=3,
)
plt.title(f"ER-positive: AUC = {auc_erPos:.2f}")
plt.xlabel("")
statistic, pvalue = stats.ttest_ind(
    res_df[(res_df["ER_status"] == "Positive") & (res_df["HRD"] == "HRD")]["HRD_score"],
    res_df[(res_df["ER_status"] == "Positive") & (res_df["HRD"] == "HR-proficient")][
        "HRD_score"
    ],
)
plt.text(0.5, plt.ylim()[1], f"p-value: {pvalue:.4f}", horizontalalignment="center")
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/SupplementaryFigures/Supp_TCGA_HRDscoreByERstatus_Positive.pdf"
)
plt.close()

plt.figure(figsize=(7, 4))
sns.boxplot(
    data=res_df[res_df["ER_status"] == "Negative"],
    x="HRD",
    y="HRD_score",
    palette="Set3",
)
sns.stripplot(
    data=res_df[res_df["ER_status"] == "Negative"],
    x="HRD",
    y="HRD_score",
    color=".3",
    size=3,
)
plt.title(f"ER-negative: AUC = {auc_erNeg:.2f}")
plt.xlabel("")
statistic, pvalue = stats.ttest_ind(
    res_df[(res_df["ER_status"] == "Negative") & (res_df["HRD"] == "HRD")]["HRD_score"],
    res_df[(res_df["ER_status"] == "Negative") & (res_df["HRD"] == "HR-proficient")][
        "HRD_score"
    ],
)
plt.text(0.5, plt.ylim()[1], f"p-value: {pvalue:.4f}", horizontalalignment="center")
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/SupplementaryFigures/Supp_TCGA_HRDscoreByERstatus_Negative.pdf"
)
plt.close()
