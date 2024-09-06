import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Load Gene-level CNA data, and organize to remove metadata
cna = pd.read_csv("~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_cna.txt", sep="\t")
cna = cna.drop_duplicates(subset="Hugo_Symbol")

# Extract cancer genes with alterations in >5% cases
cgc_data = pd.read_csv("~/Data/Census_allMon Jul  3 15_56_40 2023.csv")
cgc_data_genes = cgc_data["Gene.Symbol"].tolist()

cna = cna[cna["Hugo_Symbol"].isin(cgc_data_genes)]

# Sort dataframe
cna = cna.set_index("Hugo_Symbol")
cna = cna.drop("Entrez_Gene_Id", axis=1)
cna = cna.T
cna["Patient"] = cna.index.map(lambda x: "-".join(x.split("-")[:3]))

# Load HRD/HR-proficiency labels and merge with CNA data
ann_tcga = pd.read_pickle("Results/TCGA_HRDclassification_BRCAannotation.pkl")
ann_tcga = ann_tcga[ann_tcga["BRCA_status"].notna()]
ann_tcga = ann_tcga[ann_tcga["HRD"] == "HRD"]
ann_tcga["group"] = ann_tcga["BRCA_status"].apply(
    lambda x: "BRCA+" if x == "none" else "BRCA-defective"
)

df_ann = ann_tcga[["Patient", "group"]]

df = pd.merge(df_ann, cna, on="Patient")

# Initialize data frames to track gain/loss enrichments
fishers_gain = pd.DataFrame(
    {"Gene": df.columns[2:], "Estimate": np.nan, "pVal": np.nan}
)
fishers_loss = pd.DataFrame(
    {"Gene": df.columns[2:], "Estimate": np.nan, "pVal": np.nan}
)

# For each chromosome Gene:
#   Conduct a Fisher's exact test for HRD and HR-proficient samples against:
#     1) Gain vs not-Gain
#     2) Loss vs not-Loss
#   Save odds ratios and p-values in initialized datasets
for i, gene in enumerate(df.columns[2:]):
    gene_i = gene

    df_i = pd.DataFrame(
        {"group": df["group"], "Gain": df[gene_i] > 0, "Loss": df[gene_i] < 0}
    )
    df_i = df_i.dropna(subset=["Gain"])

    df_i["group"] = pd.Categorical(
        df_i["group"], categories=["BRCA-defective", "BRCA+"]
    )
    df_i["Gain"] = pd.Categorical(df_i["Gain"], categories=[False, True])
    df_i["Loss"] = pd.Categorical(df_i["Loss"], categories=[False, True])

    # Gains
    table_gain_i = pd.crosstab(df_i["group"], df_i["Gain"])
    fishers_gain_i = stats.fisher_exact(table_gain_i)

    fishers_gain.loc[i, "Estimate"] = fishers_gain_i[0]
    fishers_gain.loc[i, "pVal"] = fishers_gain_i[1]

    # Losses
    table_loss_i = pd.crosstab(df_i["group"], df_i["Loss"])
    fishers_loss_i = stats.fisher_exact(table_loss_i)

    fishers_loss.loc[i, "Estimate"] = fishers_loss_i[0]
    fishers_loss.loc[i, "pVal"] = fishers_loss_i[1]

# Add information to gain/loss results:
#   - Adjust p-values
#   - If Estimate == 'Inf', set to maximum estimate
#   - If Estimate == 0, set to minimum estimate
#   - log2-normalize estimates and -log10(p-adjust)
#   - Add labels for Genes with enrichment significance < .05
#   - 'Volcano plot' displaying results

# Gains
fishers_gain["padj"] = multipletests(fishers_gain["pVal"], method="fdr_bh")[1]
fishers_gain.loc[fishers_gain["Estimate"] == np.inf, "Estimate"] = (
    fishers_gain["Estimate"].replace(np.inf, np.nan).max()
)
fishers_gain.loc[fishers_gain["Estimate"] == 0, "Estimate"] = (
    fishers_gain["Estimate"].replace(0, np.nan).min()
)
fishers_gain["l2fc"] = np.log2(fishers_gain["Estimate"])
fishers_gain["logp"] = -np.log10(fishers_gain["padj"])
fishers_gain["label"] = fishers_gain["Gene"]
fishers_gain.loc[fishers_gain["padj"] > 0.05, "label"] = np.nan

plt.figure(figsize=(10, 6))
sns.scatterplot(data=fishers_gain, x="l2fc", y="logp")
for _, row in fishers_gain.dropna(subset=["label"]).iterrows():
    plt.annotate(row["label"], (row["l2fc"], row["logp"]))
plt.title("Gain vs !Gain")
plt.savefig("Figures/Supplementary/ChrGeneEnrich_GainVsNotGain.pdf")
plt.close()

# Losses
fishers_loss["padj"] = multipletests(fishers_loss["pVal"], method="fdr_bh")[1]
fishers_loss.loc[fishers_loss["Estimate"] == np.inf, "Estimate"] = (
    fishers_loss["Estimate"].replace(np.inf, np.nan).max()
)
fishers_loss.loc[fishers_loss["Estimate"] == 0, "Estimate"] = (
    fishers_loss["Estimate"].replace(0, np.nan).min()
)
fishers_loss["l2fc"] = np.log2(fishers_loss["Estimate"])
fishers_loss["logp"] = -np.log10(fishers_loss["padj"])
fishers_loss["label"] = fishers_loss["Gene"]
fishers_loss.loc[fishers_loss["padj"] > 0.05, "label"] = np.nan

plt.figure(figsize=(10, 6))
sns.scatterplot(data=fishers_loss, x="l2fc", y="logp")
for _, row in fishers_loss.dropna(subset=["label"]).iterrows():
    plt.annotate(row["label"], (row["l2fc"], row["logp"]))
plt.title("Loss vs !Loss")
plt.savefig("Figures/Supplementary/ChrGeneEnrich_LossVsNotLoss.pdf")
plt.close()

# Merge log-fold changes (+ labels, showing significant Genes for each test)
df_fishers = pd.merge(
    fishers_gain[["Gene", "l2fc", "label"]],
    fishers_loss[["Gene", "l2fc", "label"]],
    on="Gene",
    suffixes=("_GAIN", "_LOSS"),
)

# Add GAIN/LOSS/GAIN+LOSS labels for which tests are significant
df_fishers["Label"] = "none"
df_fishers.loc[df_fishers["label_GAIN"].notna(), "Label"] = "GAIN"
df_fishers.loc[df_fishers["label_LOSS"].notna(), "Label"] = "LOSS"
df_fishers.loc[
    (df_fishers["label_GAIN"].notna()) & (df_fishers["label_LOSS"].notna()), "Label"
] = "GAIN+LOSS"
df_fishers["Gene_Label"] = np.nan
df_fishers.loc[df_fishers["Label"] != "none", "Gene_Label"] = df_fishers.loc[
    df_fishers["Label"] != "none", "Gene"
]

# Plot
plt.figure(figsize=(4, 4))
sns.scatterplot(
    data=df_fishers, x="l2fc_GAIN", y="l2fc_LOSS", hue="Label", palette=["gray90"]
)
for _, row in df_fishers.dropna(subset=["Gene_Label"]).iterrows():
    plt.annotate(row["Gene_Label"], (row["l2fc_GAIN"], row["l2fc_LOSS"]))
plt.axhline(y=0, color="black", linestyle="--")
plt.axvline(x=0, color="black", linestyle="--")
plt.xlabel("log2(Fold Change) - Gain")
plt.ylabel("log2(Fold Change) - Loss")
plt.legend(title="", loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=4)
plt.savefig("Figures/Supp_GeneCnaEnrich_GainVsLoss_HRDonly.pdf", bbox_inches="tight")
plt.close()
