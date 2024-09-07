import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score

# Load data
ann_tcga = pd.read_pickle('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl')
ann_tcga = ann_tcga.set_index('Patient')[['HRD', 'BRCA_status']]
ann_tcga = ann_tcga.dropna(subset=['BRCA_status'])
ann_tcga = ann_tcga[~ann_tcga['BRCA_status'].isin(['PALB2', 'RAD51C'])]

ann_tcga['group'] = ann_tcga['BRCA_status']
ann_tcga.loc[(ann_tcga['HRD'] == 'HRD') & (ann_tcga['BRCA_status'] == 'none'), 'group'] = 'HRD_BRCA+'
ann_tcga.loc[ann_tcga['group'] == 'none', 'group'] = 'HR-proficient'

ann_tcga['HRD'] = pd.Categorical(ann_tcga['HRD'], categories=['HR-proficient', 'HRD'])
ann_tcga['group'] = pd.Categorical(ann_tcga['group'], categories=['HR-proficient', 'HRD_BRCA+', 'BRCA1', 'BRCA2'])

# Load reference testing data
Z_tumor_testing = pd.read_pickle('~/Data/TCGA/TCGA_BRCA.ExprDeconvolution_050823_p0.79_testing.pkl')

# Define test annotations
samples_intersect = list(set(Z_tumor_testing.index.str[:12]) & set(ann_tcga.index))
ann_tcga_test = ann_tcga.loc[samples_intersect]

# Load TCGA expression data
import gzip
import pickle

with gzip.open('~/Data/TCGA/TCGA_BRCA_expr_test.pkl.gz', 'rb') as f:
    expr_tumor_testing = pickle.load(f)

expr_tumor_testing = expr_tumor_testing.loc[ann_tcga_test.index]

# Load complete set of signatures
signature_centroid_list = pd.read_pickle('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl')
signature_alternative_centroid_list = pd.read_pickle('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.pkl')

signature_alternative_centroid_list = {f'Alternative_{k}': v for k, v in signature_alternative_centroid_list.items()}
signature_centroid_list.update(signature_alternative_centroid_list)

# Initialize AUC dataframe
auc_df = pd.DataFrame(columns=['Model', 'BRCA1', 'BRCA2', 'HRD_BRCApos', 'HR_BRCA_proficient', 'HRD'])

# Calculate AUCs for each model
for model_name, sig in signature_centroid_list.items():
    print(f"Processing {model_name}")
    genes_include = list(set(sig.index) & set(expr_tumor_testing.columns))
    
    results_brca_testing = pd.DataFrame({
        'group': ann_tcga_test['group'],
        'BRCA1': expr_tumor_testing[genes_include].apply(lambda x: stats.pearsonr(x, sig['BRCA1'])[0], axis=1),
        'BRCA2': expr_tumor_testing[genes_include].apply(lambda x: stats.pearsonr(x, sig['BRCA2'])[0], axis=1),
        'HRD_BRCApos': expr_tumor_testing[genes_include].apply(lambda x: stats.pearsonr(x, sig['HRD_BRCApos'])[0], axis=1),
        'HR_proficient': expr_tumor_testing[genes_include].apply(lambda x: stats.pearsonr(x, sig['HR_BRCA_proficient'])[0], axis=1)
    })
    
    results_hrd_testing = pd.DataFrame({
        'HRD_status': ann_tcga_test['HRD'],
        'HRD': expr_tumor_testing[genes_include].apply(lambda x: stats.pearsonr(x, sig['HRD'])[0], axis=1),
        'HR_proficient': expr_tumor_testing[genes_include].apply(lambda x: stats.pearsonr(x, sig['HR_proficient'])[0], axis=1)
    })
    results_hrd_testing['HRD_score'] = results_hrd_testing['HRD'] - results_hrd_testing['HR_proficient']
    
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
        'BRCA1': roc_auc_score(ann_tcga_test['group'] == 'BRCA1', expr_tumor_testing[gene]),
        'BRCA2': roc_auc_score(ann_tcga_test['group'] == 'BRCA2', expr_tumor_testing[gene]),
        'HRD_BRCApos': roc_auc_score(ann_tcga_test['group'] == 'HRD_BRCA+', expr_tumor_testing[gene]),
        'HR_BRCA_proficient': roc_auc_score(ann_tcga_test['group'] == 'HR-proficient', expr_tumor_testing[gene]),
        'HRD': roc_auc_score(ann_tcga_test['HRD'] == 'HRD', expr_tumor_testing[gene])
    }
    auc_df = auc_df.append(auc_values, ignore_index=True)

# Prepare data for plotting
auc_df_plot = auc_df.loc[[0] + list(range(4, len(auc_df))), :]
auc_df_plot['Signature'] = auc_df_plot['Model'].apply(lambda x: x.split('_')[1] if '_' in x else 'ElasticNet')
auc_df_plot = auc_df_plot.sort_values('HRD', ascending=False)
auc_df_plot['Signature'] = pd.Categorical(auc_df_plot['Signature'], categories=auc_df_plot['Signature'].unique())

# Plot HRD AUCs
plt.figure(figsize=(10, 6))
sns.barplot(x='Signature', y='HRD', data=auc_df_plot, palette='Spectral')
plt.xticks(rotation=90)
plt.ylim(0.5, 1)
plt.title('AUC for HRD Prediction')
plt.tight_layout()
plt.savefig('~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_AUCbyHRD.pdf')
plt.close()

# Plot BRCA-specific AUCs
auc_df_plot2 = pd.melt(auc_df_plot, id_vars=['Signature'], value_vars=['BRCA1', 'BRCA2', 'HRD_BRCApos', 'HR_BRCA_proficient'], var_name='group', value_name='AUC')
auc_df_plot2['group'] = pd.Categorical(auc_df_plot2['group'], categories=['BRCA1', 'BRCA2', 'HRD_BRCApos', 'HR_BRCA_proficient'])

plt.figure(figsize=(12, 6))
sns.barplot(x='Signature', y='AUC', hue='group', data=auc_df_plot2, palette='Spectral')
plt.xticks(rotation=90)
plt.ylim(0.5, 1)
plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('AUC for BRCA-specific Predictions')
plt.tight_layout()
plt.savefig('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_TCGA_AUCbyBRCA.pdf')
plt.close()