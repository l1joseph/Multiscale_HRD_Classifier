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
smc_rnaseq = smc_rnaseq.T  # Transpose to match R script

# Match samples
samples_intersect = list(set(smc_rnaseq.index) & set(results_smc_df['Patient']))
smc_rnaseq = smc_rnaseq.loc[samples_intersect]
results_smc_df = results_smc_df[results_smc_df['Patient'].isin(samples_intersect)]

# Load signature centroids
signature_centroid_list = pd.read_pickle('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl')
signature_alternative_centroid_list = pd.read_pickle('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.pkl')

# Combine signature lists
signature_alternative_centroid_list = {f'Alternative_{k}': v for k, v in signature_alternative_centroid_list.items()}
signature_centroid_list = {**signature_centroid_list, **signature_alternative_centroid_list}

# Initialize AUC dataframe
auc_df = pd.DataFrame(columns=['Model', 'BRCA1', 'BRCA2', 'HRD_BRCApos', 'HR_BRCA_proficient', 'HRD'])

# Calculate AUCs for each model
for model_name, sig in signature_centroid_list.items():
    print(f"Processing {model_name}")
    
    genes_include = list(set(sig.index) & set(smc_rnaseq.columns))
    sig = sig.loc[genes_include]
    
    # Calculate BRCA scores
    results_brca_testing = pd.DataFrame({
        'group': results_smc_df['group'],
        'BRCA1': smc_rnaseq[genes_include].apply(lambda x: stats.pearsonr(x, sig['BRCA1'])[0], axis=1),
        'BRCA2': smc_rnaseq[genes_include].apply(lambda x: stats.pearsonr(x, sig['BRCA2'])[0], axis=1),
        'HRD_BRCApos': smc_rnaseq[genes_include].apply(lambda x: stats.pearsonr(x, sig['HRD_BRCApos'])[0], axis=1),
        'HR_proficient': smc_rnaseq[genes_include].apply(lambda x: stats.pearsonr(x, sig['HR_BRCA_proficient'])[0], axis=1)
    })
    
    # Calculate HRD scores
    results_hrd_testing = pd.DataFrame({
        'HRD_status': results_smc_df['HRD'],
        'HRD': smc_rnaseq[genes_include].apply(lambda x: stats.pearsonr(x, sig['HRD'])[0], axis=1),
        'HR_proficient': smc_rnaseq[genes_include].apply(lambda x: stats.pearsonr(x, sig['HR_proficient'])[0], axis=1)
    })
    results_hrd_testing['HRD_score'] = results_hrd_testing['HRD'] - results_hrd_testing['HR_proficient']
    
    # Calculate AUCs
    auc_values = {
        'Model': model_name,
        'BRCA1': roc_auc_score(results_brca_testing['group'] == 'BRCA1', results_brca_testing['BRCA1']),
        'BRCA2': roc_auc_score(results_brca_testing['group'] == 'BRCA2', results_brca_testing['BRCA2']),
        'HRD_BRCApos': roc_auc_score(results_brca_testing['group'] == 'HRD_BRCA+', results_brca_testing['HRD_BRCApos']),
        'HR_BRCA_proficient': roc_auc_score(results_brca_testing['group'] == 'HR-proficient', results_brca_testing['HR_proficient']),
        'HRD': roc_auc_score(results_hrd_testing['HRD_status'] == 'HRD', results_hrd_testing['HRD_score'])
    }
    
    auc_df = auc_df.append(auc_values, ignore_index=True)

# Calculate AUCs for individual gene markers
gene_markers = ['BRCA1', 'BRCA2', 'POLQ', 'PARP1']
for gene in gene_markers:
    auc_values = {
        'Model': f'Gene_{gene}',
        'BRCA1': roc_auc_score(results_smc_df['group'] == 'BRCA1', smc_rnaseq[gene]),
        'BRCA2': roc_auc_score(results_smc_df['group'] == 'BRCA2', smc_rnaseq[gene]),
        'HRD_BRCApos': roc_auc_score(results_smc_df['group'] == 'HRD_BRCA+', smc_rnaseq[gene]),
        'HR_BRCA_proficient': roc_auc_score(results_smc_df['group'] == 'HR-proficient', smc_rnaseq[gene]),
        'HRD': roc_auc_score(results_smc_df['HRD'] == 'HRD', smc_rnaseq[gene])
    }
    auc_df = auc_df.append(auc_values, ignore_index=True)

# Sort AUC dataframe by HRD AUC
auc_df = auc_df.sort_values('HRD', ascending=False)
auc_df['Model'] = pd.Categorical(auc_df['Model'], categories=auc_df['Model'], ordered=True)

# Plot AUC results
plt.figure(figsize=(12, 6))
sns.barplot(x='Model', y='HRD', data=auc_df, palette='Spectral')
plt.ylim(0, 1)
plt.xticks(rotation=90)
plt.title('AUC for HRD prediction by different models')
plt.tight_layout()
plt.savefig('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_AUCresults_HRD.pdf')
plt.close()

# Plot AUC results for BRCA-specific scores
auc_df_long = pd.melt(auc_df, id_vars=['Model'], value_vars=['BRCA1', 'BRCA2', 'HRD_BRCApos', 'HR_BRCA_proficient'],
                      var_name='Signature', value_name='AUC')
plt.figure(figsize=(16, 6))
sns.barplot(x='Model', y='AUC', hue='Signature', data=auc_df_long, palette='Spectral')
plt.ylim(0, 1)
plt.xticks(rotation=90)
plt.title('AUC for BRCA-specific scores by different models')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_AUCresults_BRCA.pdf')
plt.close()