import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score

# Load SMC HRD results
results_smc_df = pd.read_pickle('~/Projects/HRD_MutationalSignature/Results/SMC_HRD_resultsSummary.pkl')

# Load BRCA status data
data_brca = pd.read_excel('~/Data/SMC_BRCA/SMC_BRCA.BrcaStatus.xlsx', skiprows=2)

# Process sample IDs
results_smc_df['sample_id'] = results_smc_df['Patient'].apply(lambda x: x[14:])
results_smc_df = pd.merge(results_smc_df, data_brca[['sample_id', 'gene_symbol']], on='sample_id', how='left')

# Assign groups
results_smc_df['group'] = results_smc_df['gene_symbol']
results_smc_df.loc[(results_smc_df['HRD_prob'] > 0.79) & (results_smc_df['gene_symbol'].isna()), 'group'] = 'HRD_BRCA+'
results_smc_df['group'] = results_smc_df['group'].fillna('HR-proficient')
results_smc_df['group'] = pd.Categorical(results_smc_df['group'], 
                                         categories=['HR-proficient', 'HRD_BRCA+', 'BRCA1', 'BRCA2'])

# Load RNA-seq data
smc_rnaseq = pd.read_csv('~/Data/SMC_BRCA/data_mrna_seq_tpm.txt', sep='\t')
smc_rnaseq = smc_rnaseq.drop_duplicates(subset='Hugo_Symbol')
smc_rnaseq = smc_rnaseq.set_index('Hugo_Symbol')
smc_rnaseq = smc_rnaseq.iloc[:, 2:]

# Match samples
samples_intersect = list(set(smc_rnaseq.columns) & set(results_smc_df['Patient']))
smc_rnaseq = smc_rnaseq[samples_intersect]
results_smc_df = results_smc_df[results_smc_df['Patient'].isin(samples_intersect)]

# Load signature centroids
signature_centroid_list = pd.read_pickle('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl')
sig_interest = signature_centroid_list['ElasticNet_alpha0.25']

# Calculate HRD and BRCA-defect scores
genes_intersect = list(set(sig_interest.index) & set(smc_rnaseq.index))
sig_interest = sig_interest.loc[genes_intersect]
smc_rnaseq_hrd = smc_rnaseq.loc[genes_intersect]

df_hrd = pd.DataFrame({
    'Patient': smc_rnaseq_hrd.columns,
    'HRD_score': smc_rnaseq_hrd.apply(lambda x: stats.pearsonr(x, sig_interest['HRD'])[0] - 
                                      stats.pearsonr(x, sig_interest['HR_proficient'])[0]),
    'BRCA1': smc_rnaseq_hrd.apply(lambda x: stats.pearsonr(x, sig_interest['BRCA1'])[0]),
    'BRCA2': smc_rnaseq_hrd.apply(lambda x: stats.pearsonr(x, sig_interest['BRCA2'])[0]),
    'HRD_BRCApos': smc_rnaseq_hrd.apply(lambda x: stats.pearsonr(x, sig_interest['HRD_BRCApos'])[0]),
    'HR_BRCA_proficient': smc_rnaseq_hrd.apply(lambda x: stats.pearsonr(x, sig_interest['HR_BRCA_proficient'])[0])
})

# Match with HRD classification
df_hrd = pd.merge(df_hrd, results_smc_df[['Patient', 'HRD', 'group']], on='Patient')

# Plot HRD by HRD
plt.figure(figsize=(6, 6))
sns.boxplot(x='HRD', y='HRD_score', data=df_hrd, palette='Set2')
plt.title('HRD Score by HRD Status')
plt.savefig('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_HRDvsHRD.pdf')
plt.close()

# Plot HRD vs BRCA
plt.figure(figsize=(8, 6))
sns.boxplot(x='group', y='HRD_score', data=df_hrd, palette='Set2')
plt.title('HRD Score by BRCA Status')
plt.savefig('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_HRDvsBRCA.pdf')
plt.close()

# Plot BRCA vs BRCA
df_hrd_plot = df_hrd.melt(id_vars=['Patient', 'group'], 
                          value_vars=['BRCA1', 'BRCA2', 'HRD_BRCApos', 'HR_BRCA_proficient'], 
                          var_name='Signature', value_name='Correlation')

plt.figure(figsize=(12, 8))
sns.boxplot(x='group', y='Correlation', hue='Signature', data=df_hrd_plot, palette='Set2')
plt.title('BRCA Scores by BRCA Status')
plt.savefig('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_BRCAvsBRCA.pdf')
plt.close()

# Calculate AUC for HRD prediction
auc_hrd = roc_auc_score(df_hrd['HRD'] == 'HRD', df_hrd['HRD_score'])
print(f"AUC for HRD prediction: {auc_hrd:.3f}")