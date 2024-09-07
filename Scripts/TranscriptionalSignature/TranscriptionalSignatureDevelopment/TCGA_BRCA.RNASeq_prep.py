import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import anndata as ad
from scipy import sparse
import scanpy as sc
import bayesprism

# Set working directory
os.chdir("~/Daniel/TCGAExpressionDeconvolution/")

# 1. Load Qian et al. 2020 data and create matrix
bc_data = sc.read_10x_mtx("Data/Qian2020/BC_counts/")
bc_sc = bc_data.X.toarray().T

# 2. Match cell.type.labels from Qian et al. metadata
bc_met = pd.read_csv("Data/Qian2020/2103-Breastcancer_metadata.csv.gz")
bc_cellMeta = bc_met["CellType"]

# 3a. Load TCGA bulk data obtained from TCGAbiolinks
tcga_brca_rnaseq = ad.read_h5ad("Data/TCGA_BRCA.counts.SE.h5ad")
tcga_brca_rnaseq = tcga_brca_rnaseq[:, ~tcga_brca_rnaseq.var_names.duplicated()]
tcga_brca_rnaseq = tcga_brca_rnaseq[:, tcga_brca_rnaseq.var["gene_name"].notna()]
tcga_brca_rnaseq = tcga_brca_rnaseq[:, ~tcga_brca_rnaseq.obs["patient"].duplicated()]

# 3b. Match to HRD/BRCA-status output and remove PALB2/RAD51C defects
ann_tcga = pd.read_pickle("Data/TCGA_HRDclassification_BRCAannotation.pkl")
ann_tcga = ann_tcga[~ann_tcga["BRCA_status"].isna()]
ann_tcga["HRD_BRCAstatus"] = ann_tcga["BRCA_status"]
ann_tcga.loc[
    (ann_tcga["HRD"] == "HRD") & (ann_tcga["BRCA_status"] == "none"), "HRD_BRCAstatus"
] = "HRD_BRCA+"
ann_tcga.loc[ann_tcga["HRD_BRCAstatus"] == "none", "HRD_BRCAstatus"] = "HR-proficient"

ann_tcga = ann_tcga[~ann_tcga["HRD_BRCAstatus"].isin(["PALB2", "RAD51C"])]

# 3c. Separate into training and testing (preserving class proportions), and save both
patients_intersect = np.intersect1d(
    ann_tcga["Patient"], tcga_brca_rnaseq.obs["patient"]
)
ann_tcga = ann_tcga[ann_tcga["Patient"].isin(patients_intersect)]
tcga_brca_rnaseq = tcga_brca_rnaseq[
    tcga_brca_rnaseq.obs["patient"].isin(patients_intersect)
]

X_train, X_test, y_train, y_test = train_test_split(
    tcga_brca_rnaseq,
    ann_tcga["HRD_BRCAstatus"],
    test_size=1 / 3,
    stratify=ann_tcga["HRD_BRCAstatus"],
    random_state=1234,
)

X_train.obs["HRD_BRCAstatus"] = y_train.values
X_test.obs["HRD_BRCAstatus"] = y_test.values

X_train.write_h5ad("Data/TCGA_BRCA.counts.SE_050823_training.h5ad")
X_test.write_h5ad("Data/TCGA_BRCA.counts.SE_050823_testing.h5ad")

bc_bulk = X_train.X.T
bc_bulk = pd.DataFrame(bc_bulk, columns=X_train.var["gene_name"])

# 4. Preprocessing of Qian et al.

# Plot single cell outliers
sc.pp.calculate_qc_metrics(bc_data, inplace=True)
sc.pl.highest_expr_genes(bc_data, n_top=20, save="_Qian2020_outlierPlot.pdf")

# Filter genes expressed in <2% cancer cells
bc_sc_cancer = bc_sc[bc_cellMeta == "Cancer", :]
index_gene_propCancer = np.mean(bc_sc_cancer > 0, axis=0) > 0.02
bc_sc = bc_sc[:, index_gene_propCancer]

# Filter outlier genes from scRNA-seq data
bc_sc_filtered = bayesprism.cleanup_genes(
    bc_sc,
    species="hs",
    gene_group=["Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"],
)
print(bc_sc_filtered.shape)

# Check concordance of varying gene types between scRNA-seq and bulk
bayesprism.plot_bulk_vs_sc(bc_sc_filtered, bc_bulk)

# Subset protein coding genes
bc_sc_filtered_pc = bayesprism.select_gene_type(
    bc_sc_filtered, gene_type="protein_coding"
)

# Construct a prism object
myPrism = bayesprism.Prism(
    reference=bc_sc_filtered_pc,
    mixture=bc_bulk,
    cell_type_labels=bc_cellMeta,
    cell_state_labels=bc_cellMeta,
    key="Cancer",
    outlier_cut=0.01,
    outlier_fraction=0.1,
)

# Run BayesPrism
tcga_brca_bayesPrism = myPrism.run_prism(n_cores=20)

pd.to_pickle(tcga_brca_bayesPrism, "TCGA_BRCA.BayesPrism_050823_p0.79_training.pkl")

# Extract cancer-specific expression
Z_tumor_training = tcga_brca_bayesPrism.get_exp(
    state_or_type="type", cell_name="Cancer"
)
pd.to_pickle(Z_tumor_training, "TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.pkl")


# notes:

# Used scanpy instead of Seurat.
# Might need to implement BayesPrism or use an alternative package/method.
