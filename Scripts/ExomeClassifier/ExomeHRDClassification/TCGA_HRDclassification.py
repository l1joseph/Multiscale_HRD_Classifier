import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import pyreadr  # for reading R data files
from sklearn.mixture import GaussianMixture
from Bio import SeqIO  # for handling genomic sequences
import mygene  # for gene ID conversion

# Set working directory
os.chdir("~/Data/TCGA/")

# Load libraries
# Note: Some R libraries used in the original script don't have direct Python equivalents.
# We'll need to implement some functionality manually or use alternative libraries.

# Load TCGA data and tally mutation contributions
# Note: We might need to implement custom functions to replicate GDCquery, GDCdownload, and GDCprepare
# I assumed data is already downloaded and prepared
tcga_mutations = pd.read_csv("tcga_mutations.csv")

# Exclude mutations from non-primary tumours
tcga_mutations["sample_type_code"] = tcga_mutations["Tumor_Sample_Barcode"].str[13:15]
tcga_mutations = tcga_mutations[tcga_mutations["sample_type_code"] == "01"]

# Note: need to implement a Python equivalent of read.maf and sig_tally
# Used placeholder functions for now
mut_maf = read_maf(
    tcga_mutations,
    is_tcga=True,
    vc_non_syn=tcga_mutations["Variant_Classification"].unique(),
)

mt_tally_tcga = sig_tally(
    mut_maf, ref_genome="BSgenome.Hsapiens.UCSC.hg38", mode="ALL", use_syn=True
)

# Check proportion of samples with SBS/ID loads >= 50 and plot distributions
sample_loads = pd.DataFrame(
    {
        "Sample": mt_tally_tcga["SBS_96"].index,
        "SBS": np.log2(mt_tally_tcga["SBS_96"].sum(axis=1) + 1),
        "ID": np.log2(mt_tally_tcga["ID_83"].sum(axis=1) + 1),
    }
)

sample_loads_long = sample_loads.melt(
    id_vars=["Sample"], var_name="MutationType", value_name="log_count"
)
sample_loads_long["MutationType"] = pd.Categorical(
    sample_loads_long["MutationType"], categories=["SBS", "ID"]
)

g_mut_loads = sns.FacetGrid(sample_loads_long, col="MutationType", col_wrap=2)
g_mut_loads.map(sns.histplot, "log_count", bins=30, kde=True)
g_mut_loads.map(plt.axvline, x=np.log2(50 + 1), color="red", linestyle="dashed")
g_mut_loads.map(
    lambda x, **kwargs: plt.axvline(x=np.median(x), color="blue", linestyle="dashed")
)
g_mut_loads.set_axis_labels("log2(count + 1)", "Frequency")
g_mut_loads.savefig(
    "~/Projects/HRD_MutationalSignature/Figures/SupplementaryFigures/Supp_TCGAmutationLoads.pdf"
)
plt.close()

print((sample_loads["SBS"] >= np.log2(50 + 1)).value_counts())
print((sample_loads["ID"] >= np.log2(50 + 1)).value_counts())

# Collate SBS_96 and ID_83 contributions and save results
tcga_muts = pd.concat([mt_tally_tcga["SBS_96"], mt_tally_tcga["ID_83"]], axis=1)

tcga_muts.to_pickle("TCGA_BRCA_mutContributions.pkl")

# Load relevant data for classifier
os.chdir("~/Projects/HRD_MutationalSignature/")

# Prior cluster mean distributions
mut_dists_mean = pd.read_pickle(
    "Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.pkl"
)

# Signature Phenotype Assignment
ann = pd.read_pickle("Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.pkl")
pheno_assigned = ann["Phenotype"]


def likelihood_calc(input_data, cluster_distributions, cluster_assign):
    log_likelihoods = pd.DataFrame(
        index=input_data.index, columns=cluster_distributions.index
    )

    print("Calculating log-likelihoods...")
    for i, row in input_data.iterrows():
        print(
            f"Calculating log likelihoods for sample {i+1} of {len(input_data)}: {row.name}"
        )
        log_likelihoods.loc[i] = cluster_distributions.apply(
            lambda x: np.sum(np.log10(x) * row), axis=1
        )

    marginal_probs = cluster_assign.value_counts() / len(cluster_assign)
    marginal_probs = marginal_probs[log_likelihoods.columns]

    log_posteriors = log_likelihoods.add(np.log10(marginal_probs), axis=1)

    final_probs = log_posteriors.apply(lambda x: 10 ** (x - x.max()), axis=1)
    final_probs = final_probs.div(final_probs.sum(axis=1), axis=0)

    return final_probs


# Apply log-likelihood approach
results_tcga_loglik = likelihood_calc(
    input_data=tcga_muts,
    cluster_distributions=mut_dists_mean,
    cluster_assign=pheno_assigned,
)

results_tcga_df = pd.DataFrame(
    {
        "Patient": results_tcga_loglik.index,
        "Phenotype_Assigned": results_tcga_loglik.idxmax(axis=1),
        "Phenotype_Assigned_prob": results_tcga_loglik.max(axis=1),
        "HRD_prob": results_tcga_loglik.filter(regex="HRD").sum(axis=1),
    }
)

results_tcga_df["HRD"] = results_tcga_df["HRD_prob"].apply(
    lambda x: "HRD" if x > 0.79 else "HR-proficient"
)
results_tcga_df.to_pickle("Results/TCGA_HRD_resultsSummary.pkl")

# Compare with BRCA defects and clinical features
valieris = pd.read_excel(
    "~/Data/TCGA/cancers-12-03687-s001/BRCA.xlsx", sheet_name="class-original"
)
valieris = valieris[valieris["BRCA1_somatic_null"] != "NA"]
valieris = valieris[
    ["sample", "event.BRCA1", "event.BRCA2", "event.RAD51C", "event.PALB2"]
]
valieris["BRCA1"] = valieris["event.BRCA1"] != 0
valieris["BRCA2"] = ~valieris["event.BRCA2"].isin([0, "Mono-allelic-inactivation"])
valieris["RAD51C"] = valieris["event.RAD51C"] != 0
valieris["PALB2"] = valieris["event.PALB2"] != 0

valieris["BRCA_status"] = "none"
valieris.loc[valieris["PALB2"], "BRCA_status"] = "PALB2"
valieris.loc[valieris["RAD51C"], "BRCA_status"] = "RAD51C"
valieris.loc[valieris["BRCA2"], "BRCA_status"] = "BRCA2"
valieris.loc[valieris["BRCA1"], "BRCA_status"] = "BRCA1"

# Clinical subtypes
clin = pd.read_csv("~/Data/TCGA/TCGA_clinicalStatus.txt", sep="\t")
tcga_clinical = pd.merge(
    valieris[["sample", "BRCA_status"]],
    clin[["bcr_patient_barcode", "er_status_by_ihc"]],
    left_on="sample",
    right_on="bcr_patient_barcode",
)
tcga_clinical = tcga_clinical.rename(columns={"er_status_by_ihc": "ER_status"})

# Create annotation
ann_tcga = pd.merge(
    results_tcga_df[["Patient", "Phenotype_Assigned", "HRD", "HRD_prob"]],
    tcga_clinical,
    left_on="Patient",
    right_on="sample",
    how="left",
)
ann_tcga = ann_tcga.drop_duplicates(subset="Patient")
ann_tcga = ann_tcga.sort_values("Phenotype_Assigned")
ann_tcga = ann_tcga.set_index("Patient")

# Plot heatmap of results
results_tcga_plot = results_tcga_loglik.loc[
    ann_tcga.index, sorted(results_tcga_loglik.columns)
]

# Sort colours
cols = sns.color_palette("Set1", n_colors=len(ann_tcga["Phenotype_Assigned"].unique()))
cols_pheno = dict(zip(sorted(ann_tcga["Phenotype_Assigned"].unique()), cols))

ann_tcga_colours = {
    "Phenotype_Assigned": cols_pheno,
    "HRD": {"HRD": "black", "HR-proficient": "white"},
    "BRCA_status": {
        "BRCA1": "blue",
        "BRCA2": "red",
        "PALB2": "gold",
        "RAD51C": "darkgreen",
        "none": "white",
    },
    "ER_status": {
        "[Not Evaluated]": "white",
        "Indeterminate": "navy",
        "Negative": "gold",
        "Positive": "darkgreen",
    },
}

# Use white -> navy scale
cols_scale = sns.color_palette("Blues", n_colors=1000)

plt.figure(figsize=(10, 10))
sns.heatmap(
    results_tcga_plot.T,
    cmap=cols_scale,
    cbar=False,
    xticklabels=False,
    yticklabels=False,
)

for i, col in enumerate(["Phenotype_Assigned", "BRCA_status", "ER_status", "HRD"]):
    plt.colorbar(
        plt.cm.ScalarMappable(
            cmap=sns.color_palette(list(ann_tcga_colours[col].values()))
        ),
        ax=plt.gca(),
        orientation="horizontal",
        aspect=10,
        pad=0.05 + i * 0.05,
    )

plt.savefig("Figures/Figure1/TCGA_HRDclassificationHeatmapExtended.pdf")
plt.close()

print(pd.crosstab(ann_tcga["HRD"], ann_tcga["BRCA_status"] != "none"))
print(pd.crosstab(ann_tcga["HRD"], ann_tcga["BRCA_status"].isin(["BRCA1", "BRCA2"])))
print(pd.crosstab(ann_tcga["HRD"], ann_tcga["BRCA_status"]))

# Plot BRCA-defect/HRD status
ann_tcga["BRCA_status_broad"] = ann_tcga["BRCA_status"].apply(
    lambda x: "BRCA+" if x == "none" else "BRCA-defective"
)

ann_tcga_plot = (
    ann_tcga[ann_tcga["BRCA_status_broad"].notna()]
    .groupby(["HRD", "BRCA_status_broad"])
    .size()
    .reset_index(name="n")
)
g_hrd_brca = sns.barplot(
    x="BRCA_status_broad", y="n", hue="HRD", data=ann_tcga_plot, palette="Paired"
)
g_hrd_brca.set(xlabel="", ylabel="% Samples")
plt.legend(title="", loc="upper right")
plt.savefig("Figures/Figure1/TCGA_BRCASensitivity.pdf")
plt.close()

# BRCA type-specific HRD classification
hrd_brca1type = ["HRD_APOBEC", "HRD_ID6mid", "HRD_ID8", "HRD_SBS8"]
hrd_brca2type = ["HRD_ID6high"]
hrd_undefined = ["HRD_ID4", "HRD_ID9"]

ann_tcga["BRCAtype_HRD"] = "HR-proficient"
ann_tcga.loc[ann_tcga["Phenotype_Assigned"].isin(hrd_brca1type), "BRCAtype_HRD"] = (
    "BRCA1-type HRD"
)
ann_tcga.loc[ann_tcga["Phenotype_Assigned"].isin(hrd_brca2type), "BRCAtype_HRD"] = (
    "BRCA2-type HRD"
)
ann_tcga.loc[ann_tcga["Phenotype_Assigned"].isin(hrd_undefined), "BRCAtype_HRD"] = (
    "HRD unassigned"
)

ann_tcga["BRCA_defect_label"] = ann_tcga["BRCA_status"].fillna("BRCA+")

ann_brca_summary = (
    ann_tcga.groupby(["BRCA_defect_label", "BRCAtype_HRD"]).size().reset_index(name="n")
)
ann_brca_summary["BRCA_defect_label"] = pd.Categorical(
    ann_brca_summary["BRCA_defect_label"], categories=["BRCA1", "BRCA2", "BRCA+"]
)
ann_brca_summary["BRCAtype_HRD"] = pd.Categorical(
    ann_brca_summary["BRCAtype_HRD"],
    categories=["HR-proficient", "HRD unassigned", "BRCA2-type HRD", "BRCA1-type HRD"],
)

g_ann_brca_summary = sns.barplot(
    x="BRCA_defect_label",
    y="n",
    hue="BRCAtype_HRD",
    data=ann_brca_summary,
    palette=["grey90", "grey50", "red", "blue"],
)
g_ann_brca_summary.set(xlabel="", ylabel="")
plt.legend(title="", loc="upper right", ncol=2)
plt.savefig("Figures/Supp_ICGCBRCAclassify.pdf")
plt.close()

# CHORD:
print(pd.crosstab(ann_tcga["BRCA_status"], ann_tcga["CHORD_type"], dropna=False))


# notes

# found this alternative for maftools: https://github.com/dentearl/mafTools
# however I don't know how fleshed out it is.
# Implement or find Python equivalents for read_maf and sig_tally functions.
# Verify that the statistical analyses (like the Fisher's exact test) are producing equivalent results to the R version.

