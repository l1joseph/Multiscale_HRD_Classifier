import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gseapy import enrichr

# Load signature centroids
signature_centroid_list = pd.read_pickle(
    "~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl"
)

# Get the gene list from ElasticNet_alpha0.25 signature
gene_list = signature_centroid_list["ElasticNet_alpha0.25"].index.tolist()

# Define the databases to use
dbs = [
    "KEGG_2021_Human",
    "WikiPathways_2019_Human",
    "Reactome_2022",
    "GO_Biological_Process_2021",
]

# Run enrichr
enrichr_results = enrichr(
    gene_list=gene_list,
    gene_sets=dbs,
    description="HRD_signature_enrichment",
    outdir="~/Projects/HRD_TranscriptionalSignature/Results/enrichr_output",
    no_plot=True,
)

# Plot results for KEGG pathway
kegg_results = enrichr_results.results["KEGG_2021_Human"]
kegg_results = kegg_results.sort_values("Adjusted P-value").head(15)

plt.figure(figsize=(10, 8))
sns.barplot(
    x="-log10(Adjusted P-value)",
    y="Term",
    data=kegg_results,
    palette="YlOrRd",
    orient="h",
)
plt.title("KEGG Pathway Enrichment")
plt.tight_layout()
plt.savefig("~/Projects/HRD_TranscriptionalSignature/Figures/Supp_GSEAenrichR.pdf")
plt.close()

# Optionally, save the full results
for db, result in enrichr_results.results.items():
    result.to_csv(f"~/Projects/HRD_TranscriptionalSignature/Results/enrichr_{db}.csv")


# notes:

# Instead of the R package enrichR, I'm using the Python package gseapy, which provides similar functionality.
# The enrichr function from gseapy is equivalent to the enrichr function in R.
# The handling of multiple databases is done within the enrichr function, which returns results for all specified databases.
