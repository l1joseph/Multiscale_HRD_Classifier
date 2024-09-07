import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Load CCLE expression data
expr = pd.read_csv("~/Data/CCLE/Expression_Public_23Q2_subsetted.csv", index_col=0)

# Load signature centroids
signature_centroid_list = pd.read_pickle(
    "~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl"
)
signature_centroid_list = {
    "ElasticNet_alpha0.25": signature_centroid_list["ElasticNet_alpha0.25"]
}

signature_alternative_centroid_list = pd.read_pickle(
    "~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.pkl"
)
signature_alternative_centroid_list = {
    f"Alternative_{k}": v for k, v in signature_alternative_centroid_list.items()
}

signature_centroid_list.update(signature_alternative_centroid_list)

# Load drug sensitivity data
drugs = pd.read_csv("~/Data/CCLE/CCLE_breast_PRISM_drugSensitivity.csv", index_col=0)
parp_inhibitors = ["OLAPARIB", "TALAZOPARIB", "NIRAPARIB", "RUCAPARIB"]
drugs = drugs[
    [col for col in drugs.columns if any(parp in col for parp in parp_inhibitors)]
]
drugs.columns = [col.split("..")[0] for col in drugs.columns]
drugs["CellLine"] = drugs.index

# Initialize results matrices
res_dict = {
    drug: pd.DataFrame(columns=["Signature", "Drug", "cor", "pVal"])
    for drug in parp_inhibitors
}

# Correlate each signature with sensitivity to each PARP inhibitor
for i, (name, sig) in enumerate(signature_centroid_list.items()):
    print(name)

    genes_intersect = np.intersect1d(sig.index, expr.columns)
    sig = sig.loc[genes_intersect]
    expr_i = expr[genes_intersect]

    # Calculate HRD scores
    expr_i_hrdScores = pd.DataFrame(
        {
            "CellLine": expr_i.index,
            "HRD": expr_i.apply(lambda x: np.corrcoef(x, sig["HRD"])[0, 1], axis=1),
            "HR_proficient": expr_i.apply(
                lambda x: np.corrcoef(x, sig["HR_proficient"])[0, 1], axis=1
            ),
        }
    )
    expr_i_hrdScores["HRD_score"] = (
        expr_i_hrdScores["HRD"] - expr_i_hrdScores["HR_proficient"]
    )

    # Merge with drug data
    df_i = pd.merge(expr_i_hrdScores, drugs, on="CellLine")

    # Calculate correlations and save in relevant results tables
    for drug in parp_inhibitors:
        cor, p = stats.pearsonr(df_i["HRD_score"], df_i[drug])
        res_dict[drug] = res_dict[drug].append(
            {"Signature": name, "Drug": drug, "cor": cor, "pVal": p}, ignore_index=True
        )

# Combine results and calculate adjusted p-values
res = pd.concat(res_dict.values())
res["logP"] = -np.log10(res["pVal"])
res["group"] = res["Signature"].apply(lambda x: x.split("_")[0])

# Plot relative significance values
plt.figure(figsize=(12, 6))
g = sns.barplot(x="Signature", y="logP", hue="group", data=res)
g.set_xticklabels(g.get_xticklabels(), rotation=90)
plt.axhline(y=-np.log10(0.05), color="red", linestyle="--")
plt.legend(title="", loc="upper right")
plt.xlabel("Signatures")
plt.ylabel("-log10(p-value)")
plt.tight_layout()
plt.savefig("~/Projects/HRD_TranscriptionalSignature/Figures/Supp_CCLEsignificance.pdf")
plt.close()

# Plot correlate of ElasticNet_alpha0.25 score against PARPi response
sig_interest = signature_centroid_list["ElasticNet_alpha0.25"]

genes_intersect = np.intersect1d(sig_interest.index, expr.columns)
sig_interest = sig_interest.loc[genes_intersect]
expr_hrd = expr[genes_intersect]

# Calculate HRD scores and match drug sensitivity data
df_hrd = pd.DataFrame(
    {
        "CellLine": expr_hrd.index,
        "HRD_score": expr_hrd.apply(
            lambda x: np.corrcoef(x, sig_interest["HRD"])[0, 1]
            - np.corrcoef(x, sig_interest["HR_proficient"])[0, 1],
            axis=1,
        ),
    }
)
df_hrd = pd.merge(df_hrd, drugs, on="CellLine")

# Reshape df_hrd for plotting
df_hrd_melt = df_hrd.melt(
    id_vars=["CellLine", "HRD_score"],
    value_vars=parp_inhibitors,
    var_name="Drug",
    value_name="PRISM",
)

# Plot correlations
plt.figure(figsize=(12, 4))
g = sns.FacetGrid(df_hrd_melt, col="Drug", col_wrap=2, height=4, aspect=1.5)
g.map(sns.regplot, "HRD_score", "PRISM")
g.map(sns.scatterplot, "HRD_score", "PRISM")
g.set_axis_labels("HRD score", "PRISM")

# Add correlation coefficient and p-value to each subplot
for ax, drug in zip(g.axes.flat, parp_inhibitors):
    r, p = stats.pearsonr(df_hrd["HRD_score"], df_hrd[drug])
    ax.text(
        0.05,
        0.95,
        f"r = {r:.2f}\np = {p:.2e}",
        transform=ax.transAxes,
        verticalalignment="top",
    )

plt.tight_layout()
plt.savefig(
    "~/Projects/HRD_TranscriptionalSignature/Figures/Figure5/CCLE_PARPiResponse.pdf"
)
plt.close()


# notes:

# The correlation calculations are done using numpy's corrcoef function and scipy's pearsonr function, which are equivalent to R's cor function.
