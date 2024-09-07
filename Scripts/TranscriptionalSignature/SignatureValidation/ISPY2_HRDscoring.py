import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Load I-SPY2 expression and clinical data from Puzstai et al. 2021
ispy2_expr = pd.read_csv(
    "~/Data/ClinicalData/ISPY2_Puzstai2021_expression.txt", sep="\t", index_col=0
)

ispy2_response = pd.read_csv(
    "~/Data/ClinicalData/GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv"
)
ispy2_response = ispy2_response[ispy2_response["Arm"] == "durvalumab/olaparib"]
ispy2_response.loc[ispy2_response["pCR.status"] == -1, "pCR.status"] = (
    0  # present in control arm
)

# Extract clinical arm from expression data
ispy2_response["ResearchID"] = "X" + ispy2_response["ResearchID"]
ispy2_expr = ispy2_expr[ispy2_response["ResearchID"]]

# Load signature centroids, extract ElasticNet_alpha0.25, and organize it with expression matrix
signature_centroid_list = pd.read_pickle(
    "~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl"
)
sig = signature_centroid_list["ElasticNet_alpha0.25"]

genes_intersect = np.intersect1d(sig.index, ispy2_expr.index)
sig = sig.loc[genes_intersect]
ispy2_expr = ispy2_expr.loc[genes_intersect]

# Calculate HRD scores and match with relevant clinical data
ispy2_hrd = pd.DataFrame(
    {
        "ResearchID": ispy2_expr.columns,
        "HRD_score": ispy2_expr.apply(
            lambda x: np.corrcoef(x, sig["HRD"])[0, 1]
            - np.corrcoef(x, sig["HR_proficient"])[0, 1]
        ),
    }
)

ispy2_hrd = pd.merge(
    ispy2_hrd,
    ispy2_response[["ResearchID", "pCR.status", "PARPi7_sig."]],
    on="ResearchID",
)

# Format pCR status into responders vs non-responders
ispy2_hrd["Response"] = ispy2_hrd["pCR.status"].map(
    {1: "Responder", 0: "Non-responder"}
)
ispy2_hrd["Response"] = pd.Categorical(
    ispy2_hrd["Response"], categories=["Non-responder", "Responder"]
)

# Plot HRD score against pCR status
plt.figure(figsize=(5, 5))
sns.boxplot(data=ispy2_hrd, x="Response", y="HRD_score", palette="Paired")
sns.swarmplot(data=ispy2_hrd, x="Response", y="HRD_score", color=".25")

statistic, pvalue = stats.ttest_ind(
    ispy2_hrd[ispy2_hrd["Response"] == "Responder"]["HRD_score"],
    ispy2_hrd[ispy2_hrd["Response"] == "Non-responder"]["HRD_score"],
)
plt.text(0.5, plt.ylim()[1], f"p-value: {pvalue:.4f}", horizontalalignment="center")

plt.title("HRD Score vs Response")
plt.xlabel("")
plt.ylabel("HRD score")
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/Figure5/ISPY2_Response.pdf"
)
plt.close()

# Plot HRD score against PARPi7 score
plt.figure(figsize=(5, 5))
sns.scatterplot(
    data=ispy2_hrd, x="PARPi7_sig.", y="HRD_score", hue="Response", palette="Paired"
)
sns.regplot(data=ispy2_hrd, x="PARPi7_sig.", y="HRD_score", scatter=False, color="gray")

r, p = stats.pearsonr(ispy2_hrd["PARPi7_sig."], ispy2_hrd["HRD_score"])
plt.text(
    0.05,
    0.95,
    f"r = {r:.2f}\np = {p:.2e}",
    transform=plt.gca().transAxes,
    verticalalignment="top",
)

plt.title("HRD Score vs PARPi7 Score")
plt.xlabel("PARPi7 Signature Score")
plt.ylabel("HRD score")
plt.legend(title="")
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/Figure5/ISPY2_HRDvsPARPi7.pdf"
)
plt.close()

# Plot PARPi7 score against pCR status
plt.figure(figsize=(5, 5))
sns.boxplot(data=ispy2_hrd, x="Response", y="PARPi7_sig.", palette="Paired")
sns.swarmplot(data=ispy2_hrd, x="Response", y="PARPi7_sig.", color=".25")

statistic, pvalue = stats.ttest_ind(
    ispy2_hrd[ispy2_hrd["Response"] == "Responder"]["PARPi7_sig."],
    ispy2_hrd[ispy2_hrd["Response"] == "Non-responder"]["PARPi7_sig."],
)
plt.text(0.5, plt.ylim()[1], f"p-value: {pvalue:.4f}", horizontalalignment="center")

plt.title("PARPi7 Score vs Response")
plt.xlabel("")
plt.ylabel("PARPi7 Signature Score")
plt.tight_layout()
plt.savefig("~/Projects/HRD_TranscriptionalSignature/Figures/Supp_ISPY2_PARPi7.pdf")
plt.close()


# notes:

# I'm using numpy's corrcoef function, which is equivalent to R's cor function for correlation calculations.
# I'm  using scipy's ttest_ind function, which is equivalent to R's t.test function for independent t-tests.
