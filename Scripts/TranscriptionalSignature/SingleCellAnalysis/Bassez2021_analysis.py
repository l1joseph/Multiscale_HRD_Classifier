import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Set working directory
os.chdir('/path/to/Data/scRNASeq/Bassez2021/')

# Load libraries and data
def load_rdata(file_path):
    # This is a placeholder function.
    # a proper R data loading mechanism, such as using the pyreadr
    pass

expr_data_bassez2021 = load_rdata('exprData_Bassez2021.Rdata')

# Load expression and metadata and subset for cancer cells
meta_bassez = pd.read_csv('1872-BIOKEY_metaData_cohort1_web.csv', header=0)
meta_bassez = meta_bassez[meta_bassez['Cell'].isin(expr_data_bassez2021.columns)]
expr_cancer = expr_data_bassez2021.loc[:, meta_bassez['cellType'] == 'Cancer_cell']

# Load signature and subset expression data
centroids = load_rdata('/path/to/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid_hrd = centroids['signature.centroid.list']['ElasticNet_alpha0.25']

genes_intersect = list(set(centroid_hrd.index) & set(expr_cancer.index))
centroid_hrd = centroid_hrd.loc[genes_intersect]
expr_cancer = expr_cancer.loc[genes_intersect]

nonzero_genes = (expr_cancer > 0).sum()
plt.figure(figsize=(10, 6))
plt.hist(nonzero_genes, bins=50)
plt.title('Distribution of Non-Zero Genes per Cell')
plt.xlabel('Number of Non-Zero Genes')
plt.ylabel('Frequency')
plt.savefig('nonzero_genes_distribution.pdf')
plt.close()

# Calculate HRD scores across cells
def calculate_correlation(x, y):
    return stats.pearsonr(x, y)[0]

hrd_scores_bassez2021 = pd.DataFrame({
    'Cell': expr_cancer.columns,
    'Sample': [x.split('_')[0] for x in expr_cancer.columns],
    'HRD': expr_cancer.apply(lambda x: calculate_correlation(x, centroid_hrd['HRD'])),
    'HR_proficient': expr_cancer.apply(lambda x: calculate_correlation(x, centroid_hrd['HR_proficient']))
})

hrd_scores_bassez2021 = hrd_scores_bassez2021.dropna(subset=['HRD', 'HR_proficient'])
hrd_scores_bassez2021['HRD_score'] = hrd_scores_bassez2021['HRD'] - hrd_scores_bassez2021['HR_proficient']

# Prep for CellphoneDB analysis
cpdb_meta = pd.DataFrame({
    'Cell': meta_bassez['Cell'],
    'CellType': meta_bassez['cellType']
})

hrd_scores_bassez2021['HRD_group'] = np.where(hrd_scores_bassez2021['HRD_score'] > 0, 'HRD', 'HR-proficient')
cpdb_meta = cpdb_meta.merge(hrd_scores_bassez2021[['Cell', 'HRD_group']], on='Cell', how='left')
cpdb_meta['HRD_group'] = cpdb_meta['HRD_group'].fillna('')
cpdb_meta['cell_type'] = cpdb_meta['CellType'] + cpdb_meta['HRD_group']
cpdb_meta = cpdb_meta[cpdb_meta['cell_type'] != 'Cancer_cell']

cpdb_meta[['Cell', 'cell_type']].to_csv('/path/to/Projects/HRD_TranscriptionalSignature/CellphoneDB/Bassez2021/Bassez2021_meta.txt', 
                                        sep='\t', index=False)

cpdb_expr = expr_data_bassez2021.loc[:, cpdb_meta['Cell']]
cpdb_expr = cpdb_expr.reset_index().rename(columns={'index': 'Gene'})
cpdb_expr = cpdb_expr[['Gene'] + list(cpdb_expr.columns[1:])]

cpdb_expr.to_csv('/path/to/Projects/HRD_TranscriptionalSignature/CellphoneDB/Bassez2021/Bassez2021_counts.txt', 
                 sep='\t', index=False)