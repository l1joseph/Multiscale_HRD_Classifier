import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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
meta_qian = meta_qian.loc[meta_qian['Cell'].isin(expr_data_qian2020.columns)]
expr_cancer = expr_data_qian2020.loc[:, meta_qian['CellType'] == 'Cancer']

# Load signature centroids
centroids = load_rdata('/path/to/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid_hrd = centroids['signature.centroid.list']['ElasticNet_alpha0.25']

# Intersect genes
genes_intersect = list(set(centroid_hrd.index) & set(expr_cancer.index))
centroid_hrd = centroid_hrd.loc[genes_intersect]
expr_cancer = expr_cancer.loc[genes_intersect]

# Calculate HRD scores across cells
def calculate_correlation(x, y):
    return stats.pearsonr(x, y)[0]

hrd_scores_qian2020 = pd.DataFrame({
    'Cell': expr_cancer.columns,
    'Sample': [x.split('_')[0] for x in expr_cancer.columns],
    'HRD': expr_cancer.apply(lambda x: calculate_correlation(x, centroid_hrd['HRD'])),
    'HR_proficient': expr_cancer.apply(lambda x: calculate_correlation(x, centroid_hrd['HR_proficient']))
})

hrd_scores_qian2020 = hrd_scores_qian2020.dropna(subset=['HRD', 'HR_proficient'])
hrd_scores_qian2020['HRD_score'] = hrd_scores_qian2020['HRD'] - hrd_scores_qian2020['HR_proficient']

# Factor results in decreasing average HRD score
qian2020_hrd_summary = hrd_scores_qian2020.groupby('Sample')['HRD_score'].mean().sort_values(ascending=False)
hrd_scores_qian2020['Sample'] = pd.Categorical(hrd_scores_qian2020['Sample'], categories=qian2020_hrd_summary.index, ordered=True)

# Plot density plots
plt.figure(figsize=(12, 8))
for sample in hrd_scores_qian2020['Sample'].cat.categories:
    sns.kdeplot(data=hrd_scores_qian2020[hrd_scores_qian2020['Sample'] == sample], x='HRD_score', label=sample)
plt.title('HRD Score Distribution by Sample')
plt.xlabel('HRD Score')
plt.ylabel('Density')
plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('HRD_score_distribution.png')
plt.close()

# Compare BRCA1-/- vs Luminal A samples
hrd_scores_plot = hrd_scores_qian2020.copy()
hrd_scores_plot['label'] = np.where(hrd_scores_plot['Sample'] == 'sc5rJUQ033', 'BRCA1-/- TNBC',
                                    np.where(hrd_scores_plot['Sample'] == 'sc5rJUQ064', 'Lum A-like', np.nan))
hrd_scores_plot = hrd_scores_plot.dropna(subset=['label'])

plt.figure(figsize=(8, 6))
sns.kdeplot(data=hrd_scores_plot, x='HRD_score', hue='label', shade=True, alpha=0.4)
plt.title('HRD Score Distribution: BRCA1-/- TNBC vs Lum A-like')
plt.xlabel('HRD Score')
plt.ylabel('Density')
plt.legend(title='Sample Type')
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Qian2020_BRCA1_LumA.pdf')
plt.close()

# Prep for CellphoneDB analysis
cpdb_meta = pd.DataFrame({
    'Cell': meta_qian['Cell'],
    'CellType': meta_qian['CellType']
})

hrd_scores_qian2020['HRD_group'] = np.where(hrd_scores_qian2020['HRD_score'] > 0, 'HRD', 'HR-proficient')
cpdb_meta = cpdb_meta.merge(hrd_scores_qian2020[['Cell', 'HRD_group']], on='Cell', how='left')
cpdb_meta['HRD_group'] = cpdb_meta['HRD_group'].fillna('')
cpdb_meta['cell_type'] = cpdb_meta['CellType'] + cpdb_meta['HRD_group']
cpdb_meta = cpdb_meta[cpdb_meta['cell_type'] != 'Cancer']

cpdb_meta[['Cell', 'cell_type']].to_csv('/path/to/Projects/HRD_TranscriptionalSignature/CellphoneDB/Qian2020/Qian2020_meta.txt', 
                                        sep='\t', index=False)

cpdb_expr = expr_data_qian2020.loc[:, cpdb_meta['Cell']]
cpdb_expr = cpdb_expr.reset_index().rename(columns={'index': 'Gene'})
cpdb_expr = cpdb_expr[['Gene'] + list(cpdb_expr.columns[1:])]

cpdb_expr.to_csv('/path/to/Projects/HRD_TranscriptionalSignature/CellphoneDB/Qian2020/Qian2020_counts.txt', 
                 sep='\t', index=False)