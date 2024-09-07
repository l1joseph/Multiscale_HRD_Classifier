import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import stats

# Set working directory
os.chdir('/path/to/Data/scRNASeq/Qian2020/')

# Load libraries and data
def load_rdata(file_path):
    # This is a placeholder function.
    # a proper R data loading mechanism, such as using the pyreadr
    pass

expr_data_qian2020 = load_rdata('exprData_Qian2020.Rdata')

# Load metadata
meta_qian = pd.read_csv('2103-Breastcancer_metadata.csv', header=0)
meta_qian = meta_qian[meta_qian['Cell'].isin(expr_data_qian2020.columns)]

# Load signature and subset expression data
centroids = load_rdata('/path/to/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid_hrd = centroids['signature.centroid.list']['ElasticNet_alpha0.25']
expr_hrd = expr_data_qian2020.loc[centroid_hrd.index, meta_qian['Cell']]

# Gene inclusion plotting
expr_nonzero = pd.DataFrame({
    'Cell': expr_hrd.columns,
    'Sample': [x.split('_')[0] for x in expr_hrd.columns],
    'prop_GenesExpressed': (expr_hrd > 0).mean() * 100
})

plt.figure(figsize=(10, 6))
sns.kdeplot(data=expr_nonzero, x='prop_GenesExpressed', hue='Sample', alpha=0.4)
plt.axvline(expr_nonzero['prop_GenesExpressed'].mean(), color='red', linestyle='dashed')
plt.xlabel('% Genes Expressed / Cell')
plt.title('Gene Expression Distribution')
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Supp_QianNonZeroGenes.pdf')
plt.close()

expr_nonzero_summary = expr_nonzero.groupby('Sample')['prop_GenesExpressed'].mean().reset_index()
plt.figure(figsize=(10, 6))
sns.histplot(data=expr_nonzero_summary, x='prop_GenesExpressed', bins=20, kde=True)
plt.xlabel('Mean % Genes Expressed Across Cells / Sample')
plt.title('Distribution of Mean Gene Expression')
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Supp_QianNonZeroGeneSummary.pdf')
plt.close()

# Calculate HRD scores across cells
def calculate_correlation(x, y):
    return stats.pearsonr(x, y)[0]

hrd_scores_qian2020 = pd.DataFrame({
    'Cell': expr_hrd.columns,
    'CellType': meta_qian['CellType'],
    'Sample': [x.split('_')[0] for x in expr_hrd.columns],
    'HRD': expr_hrd.apply(lambda x: calculate_correlation(x, centroid_hrd['HRD'])),
    'HR_proficient': expr_hrd.apply(lambda x: calculate_correlation(x, centroid_hrd['HR_proficient']))
})

hrd_scores_qian2020 = hrd_scores_qian2020.dropna(subset=['HRD', 'HR_proficient'])
hrd_scores_qian2020['HRD_score'] = hrd_scores_qian2020['HRD'] - hrd_scores_qian2020['HR_proficient']

# Factor results in decreasing average HRD score
qian2020_hrd_summary = hrd_scores_qian2020[hrd_scores_qian2020['CellType'] == 'Cancer'].groupby('Sample')['HRD_score'].mean().sort_values(ascending=False)
hrd_scores_qian2020['Sample'] = pd.Categorical(hrd_scores_qian2020['Sample'], categories=qian2020_hrd_summary.index, ordered=True)
hrd_scores_qian2020['Cancer'] = hrd_scores_qian2020['CellType'] == 'Cancer'

mu = hrd_scores_qian2020.groupby(['Sample', 'Cancer'])['HRD_score'].median().reset_index()

plt.figure(figsize=(8, 8))
sns.boxplot(data=mu, x='Cancer', y='HRD_score')
sns.pointplot(data=mu, x='Cancer', y='HRD_score', hue='Sample', 
              palette='deep', dodge=0.3, join=True, scale=0.5)
plt.ylabel('median HRD score')
plt.title('Qian et al.')
plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDinTME_Qian.pdf')
plt.close()

# Map clinical data and proportion of cells with HRD > 0
qian2020_hrd_props = hrd_scores_qian2020[hrd_scores_qian2020['Cancer']].groupby('Sample').apply(lambda x: (x['HRD_score'] > 0).mean() * 100).reset_index()
qian2020_hrd_props.columns = ['Sample', 'prop_HRD']

bc_subtypes = ['HER2', 'TN', 'B1_TN', 'TN', 'TN', 'HER2',
               'TN', 'HER2', 'Lum_HER2', 'TN', 'TN', 'TN', 'LumB', 'LumA']
qian2020_hrd_props['BC_subtype'] = pd.Categorical(bc_subtypes, 
                                                  categories=['B1_TN', 'Lum_HER2', 'HER2', 'TN', 'LumA', 'LumB'])
qian2020_hrd_props = qian2020_hrd_props.sort_values('prop_HRD')
qian2020_hrd_props['Sample'] = pd.Categorical(qian2020_hrd_props['Sample'], 
                                              categories=qian2020_hrd_props['Sample'], ordered=True)

plt.figure(figsize=(8, 8))
sns.barplot(data=qian2020_hrd_props, x='Sample', y='prop_HRD', hue='BC_subtype')
plt.xlabel('')
plt.ylabel('% HRD cells')
plt.xticks([])
plt.legend(title='BC subtype', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Qian2020_propHRDscores.pdf')
plt.close()

# Plot density plots
hrd_scores_qian2020 = hrd_scores_qian2020.merge(qian2020_hrd_props[['Sample', 'BC_subtype']], on='Sample')
plt.figure(figsize=(16, 10))
g = sns.FacetGrid(hrd_scores_qian2020[hrd_scores_qian2020['Cancer']], col='Sample', col_wrap=4, hue='BC_subtype')
g.map(sns.kdeplot, 'HRD_score', shade=True)
g.add_legend()
plt.tight_layout()
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDdensities_Qian.pdf')
plt.close()

# Qian2020 UMAP plotting
adata = sc.AnnData(expr_data_qian2020.loc[:, meta_qian['Cell'][meta_qian['CellType'] == 'Cancer']].T)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)

umap_data = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'], index=adata.obs_names)
umap_data = umap_data.merge(hrd_scores_qian2020[['Sample', 'HRD_score']], left_index=True, right_on='Cell')
umap_data = umap_data.merge(qian2020_hrd_props[['Sample', 'BC_subtype']], on='Sample')

plt.figure(figsize=(8, 8))
sns.scatterplot(data=umap_data, x='UMAP_1', y='UMAP_2', hue='BC_subtype', s=1)
plt.title('UMAP: BC Subtype')
plt.legend(title='BC subtype', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Qian2020_UMAP_BC_subtype.pdf')
plt.close()

plt.figure(figsize=(8, 8))
plt.scatter(umap_data['UMAP_1'], umap_data['UMAP_2'], c=umap_data['HRD_score'], s=1, cmap='RdBu_r')
plt.colorbar(label='HRD score')
plt.title('UMAP: HRD Score')
plt.tight_layout()
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Qian2020_UMAP_HRDscores.pdf')
plt.close()