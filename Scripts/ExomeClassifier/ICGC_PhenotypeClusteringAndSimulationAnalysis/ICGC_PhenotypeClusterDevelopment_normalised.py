import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from scipy import stats

# Load ICGC deconstructSigs data
sigs_complete = pd.read_pickle(
    "Results/ICGC_BRCA_deconstructSigs_Cutoff0.01_SBSandIDnormalised.pkl"
)

# Run mixture modelling using GaussianMixture
# Note: This is an approximation of mclust. We might need to adjust parameters.
bic_scores = []
n_components_range = range(10, 26)
for n_components in n_components_range:
    gmm = GaussianMixture(n_components=n_components, random_state=42)
    gmm.fit(sigs_complete)
    bic_scores.append(gmm.bic(sigs_complete))

plt.figure(figsize=(10, 6))
plt.plot(n_components_range, bic_scores, marker="o")
plt.xlabel("Number of components")
plt.ylabel("BIC")
plt.title("BIC Scores")
plt.savefig("Figures/Supp_ICGCSignaturesMixtureModelling.pdf")
plt.close()

# Apply optimal clustering (let's say it's 20 components)
gmm = GaussianMixture(n_components=20, random_state=42)
mod_sigs_BIC = gmm.fit(sigs_complete)
classification = mod_sigs_BIC.predict(sigs_complete)

# Form annotation document for visualisation

# CHORD
chord = pd.read_excel("~/Data/ICGC/CHORD_output.xlsx", sheet_name="CHORD")
chord = chord[~chord["group"].str.contains("HMF")]
chord = chord[["sample", "response", "hr_status", "hrd_type"]]
chord.columns = ["sample", "response", "CHORD", "CHORD_type"]

# HRDetect
hrdetect_samples = pd.read_excel(
    "~/Data/ICGC/HRDetect_sampleInfo.xlsx", sheet_name="Data", skiprows=2
)
hrdetect_pred = pd.read_excel(
    "~/Data/ICGC/HRDetect_predictions.xlsx", sheet_name="b.Predictor", skiprows=2
)
hrdetect_pred = hrdetect_pred.loc[
    hrdetect_pred["sample"].isin(hrdetect_samples["Sample"])
]

hrdetect = pd.merge(
    hrdetect_samples[
        [
            "Sample",
            "ER status",
            "BRCA1",
            "BRCA2",
            "PALB2",
            "RAD51C",
            "RAD51D",
            "isBrcaMonoallelic",
        ]
    ],
    hrdetect_pred[["sample", "predictorProb"]],
    left_on="Sample",
    right_on="sample",
)
hrdetect.loc[hrdetect["isBrcaMonoallelic"], "Gene"] = np.nan
hrdetect = hrdetect[["Sample", "ER status", "Gene", "predictorProb"]]
hrdetect.columns = ["sample", "ER_status", "BRCA_defect", "HRDetect"]
hrdetect["sample"] = hrdetect["sample"].str[:-1]

# HRDetect validation: sensitivity = 76/77 = 98.7%
print(
    pd.crosstab(
        hrdetect["BRCA_defect"].notna(), hrdetect["HRDetect"] > 0.7, dropna=False
    )
)

# Load ICGC sample data (for sampleID matching)
samples_eu = pd.read_csv("~/Data/ICGC/sample.BRCA-EU.tsv.gz", sep="\t")
samples_uk = pd.read_csv("~/Data/ICGC/sample.BRCA-UK.tsv.gz", sep="\t")
samples_icgc = pd.concat([samples_eu, samples_uk])
samples_icgc = samples_icgc[
    ["project_code", "submitted_sample_id", "icgc_donor_id"]
].drop_duplicates(subset="icgc_donor_id")

samples_icgc["final_letter"] = samples_icgc["submitted_sample_id"].str[-1]
samples_icgc = samples_icgc[samples_icgc["final_letter"].isin(["a", "b"])]
samples_icgc["sample_id"] = samples_icgc["submitted_sample_id"].str[:-1]
samples_icgc = samples_icgc.drop_duplicates(subset="sample_id")

# Combine sample data with HRDetect and CHORD
ann_icgc = pd.merge(
    samples_icgc, hrdetect, left_on="sample_id", right_on="sample", how="left"
)
ann_icgc = pd.merge(ann_icgc, chord, left_on="sample_id", right_on="sample", how="left")
ann_icgc = ann_icgc.set_index("icgc_donor_id")
ann_icgc = ann_icgc[["ER_status", "BRCA_defect", "HRDetect", "CHORD", "CHORD_type"]]

# Create annotation with finite mixture model clusters
ann = pd.DataFrame({"BIC_clust": pd.Categorical(classification)})
ann = pd.merge(ann, ann_icgc, left_index=True, right_index=True, how="left")
ann = ann.sort_values("BIC_clust")

# Order samples by classification
sigs_order = sigs_complete.loc[ann.index].T

# Set colours
colors = sns.color_palette("Set1", n_colors=len(ann["BIC_clust"].unique()))
cols_BIC = dict(zip(sorted(ann["BIC_clust"].unique()), colors))

ann_colors = {
    "BIC_clust": cols_BIC,
    "ER_status": {"positive": "darkgreen", "negative": "yellow"},
    "BRCA_defect": {"BRCA1": "blue", "BRCA2": "red"},
    "CHORD": {
        "cannot_be_determined": "grey",
        "HR_deficient": "black",
        "HR_proficient": "white",
    },
    "CHORD_type": {
        "cannot_be_determined": "grey",
        "BRCA1_type": "blue",
        "BRCA2_type": "red",
        "none": "white",
    },
}

# Plot heatmap
plt.figure(figsize=(20, 10))
sns.heatmap(
    sigs_order, cmap="Blues", xticklabels=False, yticklabels=False, cbar=False, center=0
)
for i, col in enumerate(ann.columns):
    if col in ann_colors:
        colors = pd.Series(ann[col]).map(ann_colors[col])
        plt.colorbar(
            plt.cm.ScalarMappable(cmap=plt.ListedColormap(colors.unique())),
            ax=plt.gca(),
            orientation="horizontal",
            aspect=10,
            pad=0.05 + i * 0.05,
        )
plt.title("ICGC Signatures Heatmap")
plt.tight_layout()
plt.savefig("Figures/Supp_ICGCSignaturesHeatmap.pdf")
plt.close()

# Naming clusters as signature phenotypes
pheno = [
    "SBS5_1",
    "HRD_ID8",
    "HRD_ID6high",
    "HRD_APOBEC",
    "SBS5_SBS18",
    "SBS5_2",
    "SBS5_3",
    "SBS5_4",
    "APOBEC_ID9",
    "SBS5_5",
    "HRD_SBS8",
    "APOBEC_SBS2",
    "APOBEC_SBS13",
    "SBS5_ID5",
    "SBS5_6",
    "HRD_ID9",
    "HRD_ID4",
    "HRD_ID6mid",
    "MMRD",
    "SBS5_SBS39",
]

ann["Phenotype"] = pd.Categorical(
    [pheno[i] for i in ann["BIC_clust"]], categories=sorted(pheno)
)

# Remove uninformative clusters:
ann = ann[~ann["Phenotype"].isin(["clustSmall1", "clustSmall2"])]
pheno = [p for p in pheno if p not in ["clustSmall1", "clustSmall2"]]
ann["Phenotype"] = pd.Categorical(ann["Phenotype"], categories=sorted(pheno))

ann.to_pickle("Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.pkl")

ann = ann.sort_values("Phenotype")
sigs_order = sigs_order[ann.index]

# Redo heatmap with phenotype labels
cols_Pheno = dict(
    zip(
        sorted(ann["Phenotype"].unique()),
        sns.color_palette("Set1", n_colors=len(ann["Phenotype"].unique())),
    )
)

ann_colors["Phenotype"] = cols_Pheno

plt.figure(figsize=(20, 10))
sns.heatmap(
    sigs_order, cmap="Blues", xticklabels=False, yticklabels=False, cbar=False, center=0
)
for i, col in enumerate(
    ["Phenotype", "BRCA_defect", "ER_status", "HRDetect", "CHORD", "CHORD_type"]
):
    if col in ann_colors:
        colors = pd.Series(ann[col]).map(ann_colors[col])
        plt.colorbar(
            plt.cm.ScalarMappable(cmap=plt.ListedColormap(colors.unique())),
            ax=plt.gca(),
            orientation="horizontal",
            aspect=10,
            pad=0.05 + i * 0.05,
        )
plt.title("ICGC Signatures Heatmap with Phenotypes")
plt.tight_layout()
plt.savefig("Figures/Supp_ICGCSignaturesHeatmap_withPhenotypes.pdf")
plt.close()

# Investigate BRCA-defect distribution across HRD classifications
ann_hrd = ann.copy()
ann_hrd["HRD_cluster"] = ann_hrd["Phenotype"].apply(
    lambda x: "HRD" if "HRD" in x else "HR-proficient"
)

# HRD clusters: sensitivity = 98.7%, specificity = 43.8%
print(pd.crosstab(ann_hrd["BRCA_defect"].notna(), ann_hrd["HRD_cluster"], dropna=False))

ann_hrd["BRCA_defective"] = np.where(
    ann_hrd["BRCA_defect"].isna(), "BRCA+", "BRCA_defective"
)
ann_hrd["HRD"] = np.where(
    ann_hrd["Phenotype"].str.contains("HRD"), "HRD", "HR-proficient"
)

# HRDetect > 0.7: sensitivity = 98.7%, specificity = 61.7%
print(
    pd.crosstab(ann_hrd["BRCA_defect"].notna(), ann_hrd["HRDetect"] > 0.7, dropna=False)
)

# CHORD: sensitivity = 85.3%, specificity = 66.7%
print(pd.crosstab(ann_hrd["BRCA_defect"].notna(), ann_hrd["CHORD"], dropna=False))

ann_hrd_summary = (
    ann_hrd.groupby(["HRD", "BRCA_defective"]).size().reset_index(name="n")
)
ann_hrd_summary["HRD"] = pd.Categorical(
    ann_hrd_summary["HRD"], categories=["HR-proficient", "HRD"]
)
ann_hrd_summary["BRCA_defective"] = pd.Categorical(
    ann_hrd_summary["BRCA_defective"], categories=["BRCA_defective", "BRCA+"]
)

plt.figure(figsize=(5, 6))
sns.barplot(x="BRCA_defective", y="n", hue="HRD", data=ann_hrd_summary)
plt.xlabel("")
plt.ylabel("")
plt.legend(title="", loc="upper right")
plt.title("HRD Classification")
plt.tight_layout()
plt.savefig("Figures/Supp_ICGCHRDclassify.pdf")
plt.close()

# BRCA-type specific HRD classification
hrd_brca1type = ["HRD_APOBEC", "HRD_ID6mid", "HRD_ID8", "HRD_SBS8"]
hrd_brca2type = ["HRD_ID6high"]
hrd_undefined = ["HRD_ID4", "HRD_ID9"]

ann_hrd["BRCAtype_HRD"] = "HR-proficient"
ann_hrd.loc[ann_hrd["Phenotype"].isin(hrd_brca1type), "BRCAtype_HRD"] = "BRCA1-type HRD"
ann_hrd.loc[ann_hrd["Phenotype"].isin(hrd_brca2type), "BRCAtype_HRD"] = "BRCA2-type HRD"
ann_hrd.loc[ann_hrd["Phenotype"].isin(hrd_undefined), "BRCAtype_HRD"] = "HRD unassigned"

ann_hrd["BRCA_defect_label"] = ann_hrd["BRCA_defect"].fillna("BRCA+")

ann_brca_summary = (
    ann_hrd.groupby(["BRCA_defect_label", "BRCAtype_HRD"]).size().reset_index(name="n")
)
ann_brca_summary["BRCA_defect_label"] = pd.Categorical(
    ann_brca_summary["BRCA_defect_label"], categories=["BRCA1", "BRCA2", "BRCA+"]
)
ann_brca_summary["BRCAtype_HRD"] = pd.Categorical(
    ann_brca_summary["BRCAtype_HRD"],
    categories=["HR-proficient", "HRD unassigned", "BRCA2-type HRD", "BRCA1-type HRD"],
)

plt.figure(figsize=(6, 7.5))
sns.barplot(
    x="BRCA_defect_label",
    y="n",
    hue="BRCAtype_HRD",
    data=ann_brca_summary,
    palette=["grey90", "grey50", "red", "blue"],
)
plt.xlabel("")
plt.ylabel("")
plt.legend(title="", loc="upper right", ncol=2)
plt.title("BRCA-type HRD Classification")
plt.tight_layout()
plt.savefig("Figures/Supp_ICGCBRCAclassify.pdf")
plt.close()

# CHORD:
#   BRCA1: sensitivity = 73.3%, specificity = 55.0%
#   BRCA2: sensitivity = 93.3%, specificity = 77.8%
print(pd.crosstab(ann_hrd["BRCA_defect"], ann_hrd["CHORD_type"], dropna=False))

# notes:

# I'm assuming that  input data (like sigs_complete) is in the same format as in the R version. If this isn't the case, we will need to preprocess the data differently.
