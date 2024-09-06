import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
expr_data_qian2020 = pd.read_pickle("~/Data/scRNASeq/Qian2020/exprData_Qian2020.pkl")

meta_qian = pd.read_csv("~/Data/scRNASeq/Qian2020/2103-Breastcancer_metadata.csv")
meta_qian = meta_qian.set_index("Cell").loc[expr_data_qian2020.columns]
expr_cancer = expr_data_qian2020.loc[:, meta_qian["CellType"] == "Cancer"]

# Collate genes lists
genes_g0 = pd.read_csv("~/Data/QuiescenceBiomarkers.csv")
genes_g0 = list(set(genes_g0["Genes"]) & set(expr_cancer.index))

signature_centroid_list = pd.read_pickle(
    "~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl"
)
genes_hrd = signature_centroid_list["ElasticNet_alpha0.25"].index.tolist()

genes_ddr_table = pd.read_excel("../../../../Downloads/Pearl_DDRgenes.xlsx", skiprows=1)
genes_ddr = genes_ddr_table["Gene ID"].dropna().drop_duplicates().tolist()
genes_ddr = list(set(genes_ddr) & set(expr_cancer.index))


# Summarize
def summarize_gene_expression(expr_data, gene_list, title):
    expr_fraction = (expr_data.loc[gene_list] > 0).mean(axis=1)
    print(f"\n{title}")
    print(expr_fraction.describe())

    plt.figure(figsize=(10, 6))
    sns.histplot(expr_fraction, bins=50, kde=True)
    plt.xlim(0, 1)
    plt.title(f"{title} expression")
    plt.xlabel("Fraction of cells expressing gene")
    plt.ylabel("Number of genes")
    plt.savefig(f'{title.lower().replace(" ", "_")}_expression.png')
    plt.close()


summarize_gene_expression(expr_cancer, genes_g0, "G0 genes")
summarize_gene_expression(expr_cancer, genes_ddr, "DDR genes")
summarize_gene_expression(expr_cancer, genes_hrd, "HRD genes")

# Summarize number of HRD genes expressed per cell
hrd_genes_per_cell = (expr_cancer.loc[genes_hrd] > 0).sum(axis=0)
print("\nHRD genes expressed per cell:")
print(hrd_genes_per_cell.describe())

plt.figure(figsize=(10, 6))
sns.histplot(hrd_genes_per_cell, bins=50, kde=True)
plt.title("Number of HRD genes expressed per cell")
plt.xlabel("Number of HRD genes")
plt.ylabel("Number of cells")
plt.savefig("hrd_genes_per_cell.png")
plt.close()

# Optional: all genes expression summary
# all_genes_expr = (expr_cancer > 0).mean(axis=1)
# print("\nAll genes expression:")
# print(all_genes_expr.describe())
#
# plt.figure(figsize=(10, 6))
# sns.histplot(all_genes_expr, bins=50, kde=True)
# plt.xlim(0, 1)
# plt.title('All genes expression')
# plt.xlabel('Fraction of cells expressing gene')
# plt.ylabel('Number of genes')
# plt.savefig('all_genes_expression.png')
# plt.close()
