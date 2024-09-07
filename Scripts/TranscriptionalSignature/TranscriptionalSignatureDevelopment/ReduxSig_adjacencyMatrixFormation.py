import pandas as pd
import numpy as np
from scipy import stats
import networkx as nx

# Load data
expr = pd.read_csv("../Downloads/TcgaBrca_exprHRDsig.csv", index_col=0).T

# Load important genes
imp_df = pd.read_csv("Data/imp_score_avg.csv")

genes_imp_hrd = imp_df.loc[imp_df["HRD"] > 0.7, "gene_names"].tolist()
genes_imp_hrprof = imp_df.loc[imp_df["HR.proficient"] > 0.7, "gene_names"].tolist()

genes_imp = genes_imp_hrd + genes_imp_hrprof
genes_imp_group = ["HRD"] * len(genes_imp_hrd) + ["HR-prof"] * len(genes_imp_hrprof)

expr_redux = expr[genes_imp]

# Calculate adjacency matrix
adj = pd.DataFrame(np.abs(expr_redux.corr()), index=genes_imp, columns=genes_imp)
np.fill_diagonal(adj.values, np.nan)
adj["source"] = adj.index
adj["group"] = genes_imp_group

# Melt the adjacency matrix
adj_melted = adj.melt(id_vars=["source", "group"], var_name="target", value_name="int")
adj_melted = adj_melted.dropna(subset=["int"])

# Filter connections
adj_melted = adj_melted[adj_melted["int"] > 0.0002]

# Create undirected edges
adj_melted["node1"] = adj_melted[["source", "target"]].min(axis=1)
adj_melted["node2"] = adj_melted[["source", "target"]].max(axis=1)
adj_melted["link"] = adj_melted["node1"] + "_" + adj_melted["node2"]
adj_melted = adj_melted.drop_duplicates(subset="link")

# Save adjacency list
adj_melted[["source", "target", "int"]].to_csv(
    "~/Projects/HRD_TranscriptionalSignature/Results/adjacency_reduxSig.csv",
    index=False,
)

# Save node information
adj_df = pd.DataFrame({"Gene": genes_imp, "Group": genes_imp_group})
adj_df.to_csv(
    "~/Projects/HRD_TranscriptionalSignature/Results/adjInfo_reduxSig.csv", index=False
)


# notes:
# I'm using pandas corr() function to calculate the correlation matrix instead of using the WGCNA package.

# Used melt() to reshape adj matrix into long format, which is the equivalent to pivot_longer() in R.
