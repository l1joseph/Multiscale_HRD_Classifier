import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc
import datetime

# Load data
mut_complete = pd.read_pickle("Data/BRCA_UKEU_mt_tally.pkl")
mut_complete = pd.DataFrame(
    np.column_stack([mut_complete["SBS_96"], mut_complete["ID_83"]])
)

mut_dists_mean = pd.read_pickle(
    "Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.pkl"
)

ann = pd.read_pickle("Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.pkl")
mut_complete = mut_complete.loc[ann.index]
pheno_assigned = ann["Phenotype"]


def likelihood_calc(input_data, cluster_distributions, cluster_assign):
    log_likelihoods = pd.DataFrame(
        index=input_data.index, columns=cluster_distributions.index
    )
    for i, row in input_data.iterrows():
        log_likelihoods.loc[i] = cluster_distributions.apply(
            lambda x: np.sum(np.log10(x) * row)
        )

    marginal_probs = cluster_assign.value_counts() / len(cluster_assign)
    marginal_probs = marginal_probs[log_likelihoods.columns]

    log_posteriors = log_likelihoods.add(np.log10(marginal_probs), axis=1)

    final_probs = log_posteriors.apply(lambda x: 10 ** (x - x.max()), axis=1)
    final_probs = final_probs.div(final_probs.sum(axis=1), axis=0)

    return final_probs


def simulate_likelihood_calc(
    input_data,
    cluster_distributions,
    cluster_assign,
    mutation_types="ALL",
    sample_size="",
):

    # input_data: data frame of samples (rows) and 96 trinucloetide contexts (cols)
    # sample_size: number of mutations to be sampled from patient with replacement  (THIS NEEDS TO BE ASSIGNED BELOW, AND ENSURE THAT IT IS NOT THE PSEUDO COUNT ABOVE)
    # for some reason, if sample_size is set blank, it throws an error. So, I have set it to a blank string for now.

    if mutation_types == "SBS":
        input_data = input_data.filter(regex=">")
    elif mutation_types == "ID":
        input_data = input_data.filter(regex=":")

    sims = pd.DataFrame(0, index=input_data.index, columns=input_data.columns)

    for i in range(len(input_data)):
        s = pd.Series(
            np.random.choice(
                input_data.columns,
                size=sample_size,
                p=input_data.iloc[i] / input_data.iloc[i].sum(),
            )
        ).value_counts()
        sims.loc[input_data.index[i], s.index] = s

    posteriors = likelihood_calc(sims, cluster_distributions, cluster_assign)
    posteriors["PhenoTrue"] = cluster_assign

    auc_df = pd.DataFrame(
        {"Phenotype": posteriors["PhenoTrue"].unique(), "AUC": np.nan}
    )

    for pheno in posteriors["PhenoTrue"].unique():
        roc_full = roc_curve(posteriors["PhenoTrue"] == pheno, posteriors[pheno])
        auc_val = auc(roc_full[0], roc_full[1])
        auc_df.loc[auc_df["Phenotype"] == pheno, "AUC"] = auc_val

    return {"auc_df": auc_df, "posteriors": posteriors}


def run_likelihood_sims(
    input_data,
    sample_size,
    cluster_distributions,
    cluster_assign,
    mutation_types="ALL",
    n_simulations=100,
):
    results_mat = pd.DataFrame(columns=cluster_distributions.index)
    posteriors_list = []

    for i in range(n_simulations):
        print(f"Running simulation {i+1} of {n_simulations}...", end="\r")
        post_i = simulate_likelihood_calc(
            input_data,
            cluster_distributions,
            cluster_assign,
            mutation_types,
            sample_size,
        )
        posteriors_list.append(post_i["posteriors"].iloc[:, :-1])

        auc_i = post_i["auc_df"]["AUC"]
        results_mat = results_mat.append(
            pd.Series(auc_i.values, index=cluster_distributions.index),
            ignore_index=True,
        )

    results = pd.melt(
        results_mat.reset_index(), id_vars=["index"], var_name="Pheno", value_name="AUC"
    )
    results = results.drop("index", axis=1)

    return {"results": results, "posteriors_list": posteriors_list}


# Run likelihood simulations and save/plot output
np.random.seed(123)
sig_type = "ALL"  # one of ('ALL','SBS','ID')
n_simulations = 100

for sampleSize in [25, 50, 100]:
    print(f"Running simulations: Sample Size = {sampleSize}...", flush=True)

    res_simSampleSize = run_likelihood_sims(
        input_data=mut_complete,
        sample_size=sampleSize,
        cluster_distributions=mut_dists_mean,
        cluster_assign=pheno_assigned,
        mutation_types=sig_type,
        n_simulations=n_simulations,
    )
    res_simSampleSize["results"]["sample_size"] = sampleSize

    res_results = res_simSampleSize["results"]
    res_results.to_csv(
        f"Results/ICGC_simulations_AUCs_sims{sampleSize}.txt", sep="\t", index=False
    )

    res_results["Group"] = res_results["Pheno"].apply(
        lambda x: "ID_enriched" if "ID" in x else x.split("_")[0]
    )
    res_results["Pheno"] = pd.Categorical(
        res_results["Pheno"],
        categories=sorted(res_results["Pheno"].unique(), reverse=True),
    )

    plt.figure(figsize=(10, 6))
    g_auc = sns.boxplot(data=res_results, x="Pheno", y="AUC", hue="Group", orient="h")
    plt.title(f"Sample Size: {sampleSize}")
    plt.tight_layout()
    plt.savefig(f"Figures/Supp_ICGCsimulations_PhenoReassign_sim{sampleSize}_AUCs.pdf")
    plt.close()

    res_totalPosterior = pd.concat(res_simSampleSize["posteriors_list"])
    res_totalPosterior["Pheno_Assigned"] = res_totalPosterior.idxmax(axis=1)
    res_totalPosterior["Pheno_True"] = np.tile(ann["Phenotype"], n_simulations)

    res_totalPosterior_summary = pd.crosstab(
        res_totalPosterior["Pheno_True"],
        res_totalPosterior["Pheno_Assigned"],
        normalize="index",
    )

    plt.figure(figsize=(10, 8))
    sns.heatmap(res_totalPosterior_summary, annot=True, fmt=".2f", cmap="Blues")
    plt.title(f"Sample Size: {sampleSize}")
    plt.tight_layout()
    plt.savefig(
        f"Figures/Supp_ICGCsimulations_PhenoReassign_sim{sampleSize}_posteriorHeatmap.pdf"
    )
    plt.close()


# notes:

# The progress printing has been simplified. We might want to implement a more sophisticated progress tracking system for long-running simulations.
