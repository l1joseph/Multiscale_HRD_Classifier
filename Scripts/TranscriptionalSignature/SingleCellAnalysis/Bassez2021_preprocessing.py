import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata
import os

# Set working directory
os.chdir('/path/to/Data/scRNASeq/Bassez2021/')

# Load Bassez et al. 2021 data
bc_data = sc.read_h5ad('1863-counts_cells_cohort1.h5ad')  # Assuming the data is in h5ad format

# Extract pre-treatment cells
anno = pd.read_csv('1872-BIOKEY_metaData_cohort1_web.csv', header=0)
anno = anno[anno['timepoint'] == 'Pre']

bc_data = bc_data[:, bc_data.obs_names.isin(anno['Cell'])]

# Initialize the AnnData object with the raw (non-normalized) data
adata = sc.AnnData(bc_data.X, obs=bc_data.obs, var=bc_data.var)

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, inplace=True)

# Plot QC metrics
fig, axs = plt.subplots(2, 2, figsize=(15, 15))

sns.histplot(adata.obs['total_counts'], kde=True, ax=axs[0, 0])
axs[0, 0].set_title('Total counts per cell')
axs[0, 0].set_xlabel('Total counts')

sns.histplot(adata.obs['n_genes_by_counts'], kde=True, ax=axs[0, 1])
axs[0, 1].set_title('Genes per cell')
axs[0, 1].set_xlabel('Number of genes')

sns.histplot(adata.var['total_counts'], kde=True, ax=axs[1, 0])
axs[1, 0].set_title('Total counts per gene')
axs[1, 0].set_xlabel('Total counts')

sns.histplot(adata.var['n_cells_by_counts'], kde=True, ax=axs[1, 1])
axs[1, 1].set_title('Cells per gene')
axs[1, 1].set_xlabel('Number of cells')

plt.tight_layout()
plt.savefig('QC_metrics.pdf')
plt.close()

# Plot cells ranked by their number of detected genes
plt.figure(figsize=(10, 5))
plt.plot(range(adata.n_obs), sorted(adata.obs['n_genes_by_counts']))
plt.xlabel('Cell')
plt.ylabel('Number of genes')
plt.title('Genes per cell (ordered)')
plt.savefig('Genes_per_cell_ordered.pdf')
plt.close()

# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate percentage of mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Plot QC metrics after filtering
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

sns.violinplot(data=adata.obs, y='n_genes_by_counts', ax=axs[0])
axs[0].set_title('Number of genes')

sns.violinplot(data=adata.obs, y='total_counts', ax=axs[1])
axs[1].set_title('Total counts')

sns.violinplot(data=adata.obs, y='pct_counts_mt', ax=axs[2])
axs[2].set_title('Percentage of mitochondrial genes')

plt.tight_layout()
plt.savefig('QC_metrics_after_filtering.pdf')
plt.close()

# Plot correlations between QC metrics
sc.pl.scatter(adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_mt')
plt.savefig('QC_correlations.pdf')
plt.close()

# Filter out cells with high mitochondrial content and high number of genes
adata = adata[adata.obs['n_genes_by_counts'] < 6000, :]
adata = adata[adata.obs['pct_counts_mt'] < 15, :]

# Normalize data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save preprocessed data
adata.write('exprData_Bassez2021.h5ad')

# Create a pandas DataFrame for compatibility with other scripts
expr_data_bassez2021 = pd.DataFrame(adata.X.T, index=adata.var_names, columns=adata.obs_names)
expr_data_bassez2021.to_csv('exprData_Bassez2021.csv')

print(f"Final dataset shape: {adata.shape}")


# notes:
# The CellFromTumor, PatientNumber, and CellType metadata are not added in this script. 
# If these turn out to be important for downstream analysis, we can merge this information from the annotation file.