import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

# Set working directory
os.chdir('/path/to/Data/scRNASeq/Qian2020/')

# Load Qian et al. 2020 data
bc_data = sc.read_10x_mtx(
    '/path/to/Data/scRNASeq/Qian2020/export/BC_counts/',
    var_names='gene_symbols',
    cache=True
)

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

# Load metadata
anno = pd.read_csv('2103-Breastcancer_metadata.csv')
adata.obs['CellFromTumor'] = anno.loc[anno['Cell'].isin(adata.obs_names), 'CellFromTumor'].values
adata.obs['PatientNumber'] = anno.loc[anno['Cell'].isin(adata.obs_names), 'PatientNumber'].values
adata.obs['CellType'] = anno.loc[anno['Cell'].isin(adata.obs_names), 'CellType'].values

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
adata.write('exprData_Qian2020.h5ad')

# Create a pandas DataFrame for compatibility with other scripts
expr_data_qian2020 = pd.DataFrame(adata.X.T, index=adata.var_names, columns=adata.obs_names)
expr_data_qian2020.to_csv('exprData_Qian2020.csv')

print(f"Final dataset shape: {adata.shape}")