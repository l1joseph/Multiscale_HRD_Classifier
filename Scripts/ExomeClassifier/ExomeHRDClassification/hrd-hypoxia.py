import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

# Set working directory
os.chdir("~/Data/TCGA/brca_tcga_pan_can_atlas_2018/")

# Load data
ann_tcga = pd.read_pickle(
    "~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl"
)
hypoxia = pd.read_csv("data_clinical_supp_hypoxia.txt", sep="\t")

# Filter data
ann_tcga = ann_tcga[
    ~ann_tcga["BRCA_status"].isna() & ~ann_tcga["BRCA_status"].isin(["PALB2", "RAD51C"])
]
ann_tcga = ann_tcga[ann_tcga["ER_status"].isin(["Negative", "Positive"])]

# Merge data
df = pd.merge(
    ann_tcga[["Patient", "HRD", "BRCA_status", "ER_status"]],
    hypoxia[["PATIENT_ID", "BUFFA_HYPOXIA_SCORE"]],
    left_on="Patient",
    right_on="PATIENT_ID",
)
df = df.rename(columns={"BUFFA_HYPOXIA_SCORE": "Hypoxia_Buffa"})
df["ER_status"] = df["ER_status"].map({"Negative": "ER-", "Positive": "ER+"})
df["ER_status"] = pd.Categorical(df["ER_status"], categories=["ER+", "ER-"])

df["Group"] = df["BRCA_status"]
df.loc[(df["BRCA_status"] == "none") & (df["HRD"] == "HRD"), "Group"] = "HRD_BRCA+"
df.loc[df["Group"] == "none", "Group"] = "HR-proficient"
df["Group"] = pd.Categorical(
    df["Group"], categories=["HR-proficient", "HRD_BRCA+", "BRCA1", "BRCA2"]
)

# Define comparisons
comps = [("HRD_BRCA+", "BRCA2"), ("HR-proficient", "HRD_BRCA+"), ("HRD_BRCA+", "BRCA1")]

# Create split violin plot with boxplot
plt.figure(figsize=(10, 6))

# Violin plot
sns.violinplot(
    x="Group",
    y="Hypoxia_Buffa",
    hue="ER_status",
    data=df,
    split=True,
    inner=None,
    scale="width",
    alpha=0.4,
)

# Box plot
sns.boxplot(
    x="Group",
    y="Hypoxia_Buffa",
    hue="ER_status",
    data=df,
    width=0.3,
    color="white",
    saturation=0.5,
    showfliers=False,
)

# Add statistical comparisons
y_max = df["Hypoxia_Buffa"].max()
y_increment = (y_max - df["Hypoxia_Buffa"].min()) / 20

for i, comp in enumerate(comps):
    group1 = df[df["Group"] == comp[0]]["Hypoxia_Buffa"]
    group2 = df[df["Group"] == comp[1]]["Hypoxia_Buffa"]
    t, p = stats.ttest_ind(group1, group2)
    x1 = df["Group"].cat.categories.get_loc(comp[0])
    x2 = df["Group"].cat.categories.get_loc(comp[1])
    y = y_max + (i + 1) * y_increment
    plt.plot(
        [x1, x1, x2, x2],
        [y, y + y_increment / 2, y + y_increment / 2, y],
        color="black",
    )
    plt.text((x1 + x2) / 2, y + y_increment / 2, f"p={p:.3f}", ha="center", va="bottom")

plt.ylabel("Hypoxia Score (Buffa)")
plt.title("Hypoxia Score by Group and ER Status")
plt.legend(title="")
plt.tight_layout()
plt.savefig("~/Projects/Thesis/Chapter 4/HRD_Hypoxia.pdf")
plt.close()

# Plot for ER+ samples only
plt.figure(figsize=(8, 6))
sns.boxplot(x="Group", y="Hypoxia_Buffa", data=df[df["ER_status"] == "ER+"])
plt.title("Hypoxia Score (Buffa) for ER+ Samples")

# Add statistical comparisons for ER+ samples
y_max = df[df["ER_status"] == "ER+"]["Hypoxia_Buffa"].max()
y_increment = (y_max - df[df["ER_status"] == "ER+"]["Hypoxia_Buffa"].min()) / 20

for i, comp in enumerate(comps):
    group1 = df[(df["Group"] == comp[0]) & (df["ER_status"] == "ER+")]["Hypoxia_Buffa"]
    group2 = df[(df["Group"] == comp[1]) & (df["ER_status"] == "ER+")]["Hypoxia_Buffa"]
    t, p = stats.ttest_ind(group1, group2)
    x1 = df["Group"].cat.categories.get_loc(comp[0])
    x2 = df["Group"].cat.categories.get_loc(comp[1])
    y = y_max + (i + 1) * y_increment
    plt.plot(
        [x1, x1, x2, x2],
        [y, y + y_increment / 2, y + y_increment / 2, y],
        color="black",
    )
    plt.text((x1 + x2) / 2, y + y_increment / 2, f"p={p:.3f}", ha="center", va="bottom")

plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_MutationalSignature/Figures/Figure2/HRD_Hypoxia_ERpositive.pdf"
)
plt.close()

# Linear model and ANOVA
formula = "Hypoxia_Buffa ~ C(ER_status) + C(BRCA_status) * C(HRD)"
model = ols(formula, data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print("\nANOVA Results:")
print(anova_table)
