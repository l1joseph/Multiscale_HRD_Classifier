import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load mutation tallies
mt_tally_brca_wgs = pd.read_pickle(
    "~/Projects/HRD_MutationalSignature/Data/BRCA_UKEU_mt_tally.pkl"
)
mut_complete = pd.concat(
    [mt_tally_brca_wgs["SBS_96"], mt_tally_brca_wgs["ID_83"]], axis=1
)

# Load annotation (ann)
ann = pd.read_pickle("Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.pkl")
mut_complete = mut_complete.loc[ann.index]

# Separate mutation counts into clusters
mut_byClust = {
    pheno: mut_complete[ann["Phenotype"] == pheno]
    for pheno in ann["Phenotype"].unique()
}


def collate_function(input_data, variant_type="ALL", collation="total"):
    sbs_index = input_data.columns.str.contains(">")
    id_index = input_data.columns.str.contains(":")

    if variant_type == "SBS":
        input_final = input_data.loc[:, sbs_index]
    elif variant_type == "ID":
        input_final = input_data.loc[:, id_index]
    else:
        input_final = input_data

    # Add pseudocount
    input_final = input_final + 1  # original
    if collation == "mean":
        dist_temp = input_final.div(input_final.sum(axis=1), axis=0)
        dist = dist_temp.mean()
    else:
        if collation != "total":
            print("Set collation to mean or total. Setting to total...")
        dist_temp = input_final.sum()
        dist = dist_temp / dist_temp.sum()

    return dist


mut_dists_total = pd.DataFrame()
mut_dists_mean = pd.DataFrame()

for pheno in mut_byClust.keys():
    mut_dists_total = mut_dists_total.append(
        collate_function(mut_byClust[pheno]), ignore_index=True
    )
    mut_dists_mean = mut_dists_mean.append(
        collate_function(mut_byClust[pheno], collation="mean"), ignore_index=True
    )

mut_dists_total.index = mut_byClust.keys()
mut_dists_mean.index = mut_byClust.keys()

# Save prior clusters
mut_dists_mean.to_pickle("../Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.pkl")
mut_dists_total.to_pickle(
    "../Data/ClusterLikelihoods/ICGC_clust20_mclust_totalCont.pkl"
)

# Plot prior clusters
mut_dists_mean_plot = pd.concat(
    [mut_dists_mean.iloc[:, :96], mut_dists_mean.iloc[:, 96:]], axis=1
)
mut_dists_mean_plot["Phenotype"] = mut_dists_mean_plot.index
mut_dists_mean_plot = mut_dists_mean_plot.melt(
    id_vars=["Phenotype"], var_name="Context", value_name="Contribution"
)
mut_dists_mean_plot["Phenotype"] = pd.Categorical(
    mut_dists_mean_plot["Phenotype"], categories=sorted(mut_dists_mean.index)
)

# Sort indel context order
signatures_id83 = pd.read_csv("~/Data/COSMIC_v3.3_ID_GRCh37.txt", sep="\t")
mut_dists_mean_plot["Context"] = pd.Categorical(
    mut_dists_mean_plot["Context"],
    categories=list(signatures_cosmic.columns) + list(signatures_id83["Type"]),
)
mut_dists_mean_plot["Type"] = np.where(
    mut_dists_mean_plot["Context"].str.contains(">"), "SBS", "indel"
)

g_meanPlot = sns.FacetGrid(
    mut_dists_mean_plot, col="Phenotype", col_wrap=4, height=3, aspect=1.5
)
g_meanPlot.map(sns.barplot, "Context", "Contribution", "Type", palette="Set1")
g_meanPlot.set_xticklabels(rotation=90)
g_meanPlot.add_legend()
plt.tight_layout()
plt.savefig("../Figures/SupplementaryFigures/Supp_LikelihoodDistributionsMeans.pdf")
plt.close()


# notes:

# The script assumes that signatures_cosmic is defined.
# Need to make sure this variable is properly loaded or defined before running this script.
