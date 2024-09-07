import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# Load libraries and data
def load_rdata(file_path):
    # This is a placeholder function.
    # a proper R data loading mechanism, such as using the pyreadr
    pass

# Load signature
centroids = load_rdata('/path/to/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid_hrd = centroids['signature.centroid.list']['ElasticNet_alpha0.25']

# Load expression data, match with signature, and log2-normalise
expr_chung = pd.read_csv('/path/to/Data/scRNASeq/Chung2017/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt', sep='\t')
expr_chung = expr_chung.set_index('gene_name')

genes_intersect = list(set(centroid_hrd.index) & set(expr_chung.index))

centroid_hrd = centroid_hrd.loc[genes_intersect]
expr_chung = expr_chung.loc[genes_intersect]

expr_chung = np.log2(expr_chung + 1)

# Separate into bulk and scRNAseq
expr_chung_bulk = expr_chung.loc[:, expr_chung.columns.str.contains('Pooled')]
expr_chung_sc = expr_chung.iloc[:, 14:]

# Extract only tumour cells from single cell data
chung_info = pd.read_csv('/path/to/Data/scRNASeq/Chung2017/GSE75688_final_sample_information.txt', sep='\t')
cells_tumor = chung_info.loc[(chung_info['type'] == 'SC') & (chung_info['index'] == 'Tumor'), 'sample']

expr_chung_sc = expr_chung_sc.loc[:, expr_chung_sc.columns.isin(cells_tumor)]

# Calculate HRD scores
def calculate_correlation(x, y):
    return stats.pearsonr(x, y)[0]

hrd_scores_bk = pd.DataFrame({
    'sample': [x.split('_')[0] for x in expr_chung_bulk.columns],
    'HRD': expr_chung_bulk.apply(lambda x: calculate_correlation(x, centroid_hrd['HRD'])),
    'HR_proficient': expr_chung_bulk.apply(lambda x: calculate_correlation(x, centroid_hrd['HR_proficient']))
})
hrd_scores_bk['HRD_score_bulk'] = hrd_scores_bk['HRD'] - hrd_scores_bk['HR_proficient']

hrd_scores_sc = pd.DataFrame({
    'cell': expr_chung_sc.columns,
    'HRD': expr_chung_sc.apply(lambda x: calculate_correlation(x, centroid_hrd['HRD'])),
    'HR_proficient': expr_chung_sc.apply(lambda x: calculate_correlation(x, centroid_hrd['HR_proficient']))
})
hrd_scores_sc['HRD_score_sc'] = hrd_scores_sc['HRD'] - hrd_scores_sc['HR_proficient']
hrd_scores_sc['sample'] = [x.split('_')[0] for x in hrd_scores_sc['cell']]

# Match bulk and single-cell HRD scores
hrd_scores_sc_summary = hrd_scores_sc.groupby('sample')['HRD_score_sc'].mean().reset_index()

hrd_scores_df = pd.merge(hrd_scores_bk[['sample', 'HRD_score_bulk']], 
                         hrd_scores_sc_summary, 
                         on='sample')

# Plot results
plt.figure(figsize=(8, 6))
sns.regplot(x='HRD_score_bulk', y='HRD_score_sc', data=hrd_scores_df)
plt.xlabel('Bulk HRD')
plt.ylabel('Mean single-cell HRD')
plt.title('Bulk vs Single-cell HRD Scores')

# Add correlation coefficient and p-value to the plot
corr, p_value = stats.pearsonr(hrd_scores_df['HRD_score_bulk'], hrd_scores_df['HRD_score_sc'])
plt.text(0.05, 0.95, f'r = {corr:.2f}\np = {p_value:.2e}', 
         transform=plt.gca().transAxes, verticalalignment='top')

plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Chung_BulkSingleCell.pdf')
plt.close()