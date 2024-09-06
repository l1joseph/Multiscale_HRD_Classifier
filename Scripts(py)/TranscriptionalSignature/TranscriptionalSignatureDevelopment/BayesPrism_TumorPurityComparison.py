import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Load BayesPrism results
theta = pd.read_pickle("~/Data/TCGA/TCGA_BRCA.BayesPrism.theta.pkl")
df_theta = pd.DataFrame(
    {"Sample.ID": [x[:16] for x in theta.index], "CancerCellFraction": theta["Cancer"]}
)

# Load ESTIMATE scores
Tumor_purity = pd.read_csv("~/Data/TCGA/TCGA_purity_ESTIMATE.txt", sep="\t")

# Merge data
df_theta = pd.merge(df_theta, Tumor_purity[["Sample.ID", "ESTIMATE"]], on="Sample.ID")


# Function for converting ESTIMATE column
def comma_function(x):
    return float(str(x).replace(",", "."))


df_theta["TumorPurity"] = df_theta["ESTIMATE"].apply(comma_function)

# Compare estimates
plt.figure(figsize=(10, 8))
sns.scatterplot(data=df_theta, x="TumorPurity", y="CancerCellFraction", alpha=0.3)
sns.regplot(
    data=df_theta, x="TumorPurity", y="CancerCellFraction", scatter=False, color="red"
)

# Add correlation coefficient and p-value
r, p = stats.pearsonr(df_theta["TumorPurity"], df_theta["CancerCellFraction"])
plt.text(
    0.05,
    0.95,
    f"r = {r:.2f}\np = {p:.2e}",
    transform=plt.gca().transAxes,
    verticalalignment="top",
)

plt.title("BayesPrism vs ESTIMATE Tumor Purity")
plt.xlabel("ESTIMATE Tumor Purity")
plt.ylabel("BayesPrism Cancer Cell Fraction")
plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/Supp_BayesPrismEstimates.pdf"
)
plt.close()
