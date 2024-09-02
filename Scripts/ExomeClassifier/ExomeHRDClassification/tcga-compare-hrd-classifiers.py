import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc

# Set working directory
os.chdir('~/Projects/HRD_MutationalSignature/Results')

# Load data
results_tcga_df = pd.read_pickle('TCGA_HRD_resultsSummary.pkl')
results_tcga_df = pd.DataFrame({
    'Tumor': results_tcga_df['Patient'],
    'p_HRDi': results_tcga_df['HRD_prob'],
    'HRDi': results_tcga_df['HRD']
})

sigma = pd.read_csv('SigMA_TCGA_output.csv')
sigma['Tumor'] = sigma['tumor'].str[:12]
sigma = sigma[['Tumor', 'Signature_3_mva', 'pass_mva']]
sigma.columns = ['Tumor', 'p_HRD_SigMA', 'SigMA']

sbs3 = pd.read_csv('TCGA_BRCA_deconstructSigs.txt', sep='\t')
sbs3['Tumor'] = sbs3.index.str[:12]
sbs3['SBS3_dominant'] = sbs3.idxmax(axis=1) == 'SBS3'
sbs3 = sbs3[['Tumor', 'SBS3', 'SBS3_dominant']]
sbs3.columns = ['Tumor', 'p_HRD_SBS3', 'SBS3_dominant']

mac_full = pd.read_csv('~/Data/TCGA/Macintyre_TCGAExposures.txt', sep='\t')
mac_full = mac_full[mac_full['Cancer'] == 'BRCA']
mac_full = mac_full.pivot(index='Sample', columns='Signature', values='Exposure')
mac_full['Tumor'] = mac_full.index
mac_full['p_HRD_CX3'] = mac_full['CX3']
mac_full['CX3_dominant'] = mac_full.idxmax(axis=1) == 'CX3'
mac_full = mac_full[['Tumor', 'p_HRD_CX3', 'CX3_dominant']]

marq = pd.read_csv('~/Data/TCGA/Marquard_HRDScores.txt', sep='\t')
marq['HRD_index'] = marq[['NtAI', 'LST', 'HRD.LOH']].sum(axis=1)
marq['HRDindex_42'] = marq['HRD_index'] >= 42
marq['HRDindex_63'] = marq['HRD_index'] > 63
marq = marq[['Tumor', 'HRD_index', 'HRDindex_42', 'HRDindex_63']]

# Merge all dataframes
df_full = results_tcga_df.merge(sigma, on='Tumor', how='left')
df_full = df_full.merge(sbs3, on='Tumor', how='left')
df_full = df_full.merge(mac_full, on='Tumor', how='left')
df_full = df_full.merge(marq, on='Tumor', how='left')

# Add Valieris BRCA defects
valieris = pd.read_excel('~/Data/TCGA/cancers-12-03687-s001/BRCA.xlsx', sheet_name='class-original')
valieris = valieris[valieris['BRCA1_somatic_null'] != 'NA']
valieris = valieris[['sample', 'event.BRCA1', 'event.BRCA2', 'event.RAD51C', 'event.PALB2']]
valieris['BRCA1'] = valieris['event.BRCA1'] != 0
valieris['BRCA2'] = ~valieris['event.BRCA2'].isin([0, 'Mono-allelic-inactivation'])
valieris['RAD51C'] = valieris['event.RAD51C'] != 0
valieris['PALB2'] = valieris['event.PALB2'] != 0

valieris['BRCA_status'] = 'none'
valieris.loc[valieris['PALB2'], 'BRCA_status'] = 'PALB2'
valieris.loc[valieris['RAD51C'], 'BRCA_status'] = 'RAD51C'
valieris.loc[valieris['BRCA2'], 'BRCA_status'] = 'BRCA2'
valieris.loc[valieris['BRCA1'], 'BRCA_status'] = 'BRCA1'

valieris = valieris[['sample', 'BRCA_status']]
valieris.columns = ['Tumor', 'BRCA_defective']
valieris['BRCA_defective'] = valieris['BRCA_defective'] != 'none'

df_full = df_full.merge(valieris, on='Tumor', how='left')

# Calculate performance metrics
def calculate_metrics(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)
    f_score = 2 * (precision * recall) / (precision + recall)
    return pd.Series({'TP': tp, 'TN': tn, 'FP': fp, 'FN': fn, 'Recall': recall, 'Precision': precision, 'F_score': f_score})

methods = ['HRDi', 'SigMA', 'SBS3_dominant', 'CX3', 'HRDindex_42', 'HRDindex_63']
results = pd.DataFrame(index=methods)

for method in methods:
    results.loc[method] = calculate_metrics(df_full['BRCA_defective'], df_full[method])

# Plot results
plt.figure(figsize=(9, 4.5))
sns.barplot(data=results.reset_index().melt(id_vars='index', value_vars=['Recall', 'Precision', 'F_score']),
            x='index', y='value', hue='variable')
plt.title('Performance Metrics for HRD Classifiers')
plt.xlabel('Method')
plt.ylabel('Score')
plt.xticks(rotation=45)
plt.legend(title='Metric')
plt.tight_layout()
plt.savefig('../Figures/Supp_TCGAcomparisons_SensSpec.pdf')
plt.close()

# Plot ROC curves
plt.figure(figsize=(6, 6))
for method in ['p_HRDi', 'p_HRD_SigMA', 'p_HRD_SBS3', 'p_HRD_CX3', 'HRD_index']:
    fpr, tpr, _ = roc_curve(df_full['BRCA_defective'], df_full[method])
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, label=f'{method} (AUC = {roc_auc:.2f})')

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.savefig('../Figures/Supp_TCGAcomparisons_ROC.pdf')
plt.close()
