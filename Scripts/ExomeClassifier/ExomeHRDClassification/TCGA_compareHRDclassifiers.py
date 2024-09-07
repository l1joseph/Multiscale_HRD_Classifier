import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc



# Load HRD classification results
results_tcga_df = pd.read_pickle('Results/TCGA_HRD_resultsSummary.pkl')

# Load SigMA results (run using webtool)
sigma = pd.read_csv('Results/SigMA_TCGA_output.csv')
sigma['Tumor'] = sigma['tumor'].str[:12]
sigma = sigma[['Tumor', 'Signature_3_mva', 'pass_mva']]
sigma = sigma.rename(columns={'Signature_3_mva': 'p_HRD_SigMA', 'pass_mva': 'SigMA'})

# Load deconstructSigs results (run manually)
sbs3 = pd.read_csv('Results/TCGA_BRCA_deconstructSigs.txt', sep='\t', index_col=0)
sbs3['Tumor'] = sbs3.index.str[:12]
sbs3['p_HRD_SBS3'] = sbs3['SBS3']
sbs3['SBS3_dominant'] = (sbs3.iloc[:, :-2].idxmax(axis=1) == 'SBS3').astype(int)
sbs3 = sbs3[['Tumor', 'p_HRD_SBS3', 'SBS3_dominant']]

# Load CX3 Copy Number Signatures (taken from Drews et al. Nature 2022)
mac_full = pd.read_csv('~/Data/TCGA/Macintyre_TCGAExposures.txt', sep='\t')
mac_full = mac_full[mac_full['Cancer'] == 'BRCA']
mac_full = mac_full.pivot(index='Sample', columns='Signature', values='Exposure')
mac_full['p_HRD_CX3'] = mac_full['CX3']
mac_full['CX3_dominant'] = (mac_full.iloc[:, :-1].idxmax(axis=1) == 'CX3').astype(int)
mac_full = mac_full.reset_index()[['Sample', 'p_HRD_CX3', 'CX3_dominant']]
mac_full = mac_full.rename(columns={'Sample': 'Tumor'})

# Load HRD index scores (taken from Marquard et al. 2015)
marq = pd.read_csv('~/Data/TCGA/Marquard_HRDScores.txt', sep='\t')
marq['HRD_index'] = marq[['NtAI', 'LST', 'HRD-LOH']].sum(axis=1)
marq['HRDindex_42'] = (marq['HRD_index'] >= 42).astype(int)
marq['HRDindex_63'] = (marq['HRD_index'] > 63).astype(int)
marq = marq[['Tumor', 'HRD_index', 'HRDindex_42', 'HRDindex_63']]

# Join all dataframes
df_full = results_tcga_df.merge(sigma, on='Tumor', how='left')
df_full = df_full.merge(sbs3, on='Tumor', how='left')
df_full = df_full.merge(mac_full, on='Tumor', how='left')
df_full = df_full.merge(marq, on='Tumor', how='left')

# Add Valieris BRCA defects
valieris = pd.read_excel('~/Data/TCGA/cancers-12-03687-s001/BRCA.xlsx', sheet_name='class-original')
valieris = valieris[valieris['BRCA1_somatic_null'] != 'NA']
valieris['BRCA1'] = valieris['event.BRCA1'] != 0
valieris['BRCA2'] = ~(valieris['event.BRCA2'].isin([0, 'Mono-allelic-inactivation']))
valieris['RAD51C'] = valieris['event.RAD51C'] != 0
valieris['PALB2'] = valieris['event.PALB2'] != 0

valieris['BRCA_status'] = 'none'
valieris.loc[valieris['PALB2'], 'BRCA_status'] = 'PALB2'
valieris.loc[valieris['RAD51C'], 'BRCA_status'] = 'RAD51C'
valieris.loc[valieris['BRCA2'], 'BRCA_status'] = 'BRCA2'
valieris.loc[valieris['BRCA1'], 'BRCA_status'] = 'BRCA1'

valieris = valieris[['sample', 'BRCA_status']]
valieris['BRCA_defective'] = (valieris['BRCA_status'] != 'none').astype(int)
valieris = valieris.rename(columns={'sample': 'Tumor'})

df_full = df_full.merge(valieris, on='Tumor', how='left')

# Calculate sensitivities, specificities, and F-scores
def break_table_func(column):
    t = pd.crosstab(df_full[column], df_full['BRCA_defective'])
    true_negative = t.loc[0, 0]
    false_negative = t.loc[0, 1]
    false_positive = t.loc[1, 0]
    true_positive = t.loc[1, 1]
    return pd.Series([true_positive, true_negative, false_positive, false_negative])

results = pd.DataFrame({
    'Method': ['HRDi', 'SigMA', 'SBS3_dominant', 'CX3', 'HRDindex_42', 'HRDindex_63'],
    'TP': np.nan, 'TN': np.nan, 'FP': np.nan, 'FN': np.nan
})

for i, method in enumerate(results['Method']):
    results.iloc[i, 1:] = break_table_func(method)

# Calculate sensitivity, specificity, and balanced F-score of each method
weight = 1

results['recall'] = results['TP'] / (results['TP'] + results['FN'])
results['precision'] = results['TP'] / (results['TP'] + results['FP'])
results['F_score'] = ((1 + weight**2) * results['precision'] * results['recall']) / ((weight**2 * results['precision']) + results['recall'])

results['Method'] = pd.Categorical(results['Method'], 
                                   categories=results.sort_values('F_score', ascending=False)['Method'])

# Plot results
plt.figure(figsize=(6, 6))
sns.barplot(x='Method', y='F_score', data=results, palette='Paired')
plt.title('F-scores of HRD Classifiers')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('../Figures/Figure1/TCGA_HRDcomparisonsFscores.pdf')
plt.close()

# Plot sensitivity and precision
results_plot = results.melt(id_vars='Method', value_vars=['recall', 'precision'], 
                            var_name='Measure', value_name='value')
plt.figure(figsize=(8, 4.8))
sns.barplot(x='Method', y='value', hue='Measure', data=results_plot, palette='Paired')
plt.title('Sensitivity and Precision of HRD Classifiers')
plt.xticks(rotation=45)
plt.legend(title='')
plt.tight_layout()
plt.savefig('../Figures/Supp_TCGAcomparisonsSensSpec.pdf')
plt.close()

# Plot comparative AUC curves
plt.figure(figsize=(8, 6))

for column in ['p_HRDi', 'p_HRD_SigMA', 'p_HRD_SBS3', 'p_HRD_CX3', 'HRD_index']:
    fpr, tpr, _ = roc_curve(df_full['BRCA_defective'], df_full[column])
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, label=f'{column} (AUC = {roc_auc:.2f})')

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig('../Figures/Supp_TCGAcomparisonsROC.pdf')
plt.close()