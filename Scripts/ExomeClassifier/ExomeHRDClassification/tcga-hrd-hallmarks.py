import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Set working directory
os.chdir('~/Projects/HRD_MutationalSignature/')

# Load TCGA HRD classifications and hallmarks, and combine
ann_tcga = pd.read_pickle('Results/TCGA_HRDclassification_BRCAannotation.pkl')
tcga_HRDHallmarks = pd.read_csv('~/Data/TCGA/HRD_hallmarks.txt', sep='\t')

resultsHRD = pd.merge(ann_tcga[['Patient', 'Phenotype_Assigned', 'HRD']], 
                      tcga_HRDHallmarks, how='left')
resultsHRD['HRD'] = pd.Categorical(resultsHRD['HRD'], categories=['HR-proficient', 'HRD'], ordered=True)

# MYC amplification
resultsHRD['MYC_status'] = pd.Categorical(resultsHRD['MYC_status'], 
                                          categories=['Deletion', 'Normal', 'Amplification'], 
                                          ordered=True)
res_MYC = resultsHRD.groupby(['HRD', 'MYC_status']).size().reset_index(name='n')
res_MYC = res_MYC.dropna()

plt.figure(figsize=(6, 3))
sns.barplot(data=res_MYC, x='HRD', y='n', hue='MYC_status')
plt.ylabel('% Samples')
plt.legend(title='')
plt.savefig('Figures/Figure2/HRD_vs_MYCamplification.pdf')
plt.close()

# POLQ expression
plt.figure(figsize=(4, 4))
sns.boxplot(data=resultsHRD, x='HRD', y='log_POLQ_FPKM')
plt.xlabel('')
statannot.add_stat_annotation(data=resultsHRD, x='HRD', y='log_POLQ_FPKM', 
                              box_pairs=[('HR-proficient', 'HRD')],
                              test='t-test_ind', text_format='p={:.3e}',
                              loc='outside', verbose=2)
plt.savefig('Figures/Figure2/HRD_vs_POLQexpression.pdf')
plt.close()

# HRD index scores
plt.figure(figsize=(4, 4))
sns.boxplot(data=resultsHRD, x='HRD', y='HRD_index')
plt.xlabel('')
statannot.add_stat_annotation(data=resultsHRD, x='HRD', y='HRD_index', 
                              box_pairs=[('HR-proficient', 'HRD')],
                              test='t-test_ind', text_format='p={:.3e}',
                              loc='outside', verbose=2)
plt.savefig('Figures/Figure2/HRD_vs_HRDscore.pdf')
plt.close()

# Individual HRD index scores
df_hrdIndex = resultsHRD[['HRD', 'NtAI', 'LST', 'HRD.LOH']].melt(id_vars=['HRD'], 
                                                                 var_name='HRD_index', 
                                                                 value_name='score')
plt.figure(figsize=(8, 4))
sns.boxplot(data=df_hrdIndex, x='HRD', y='score', hue='HRD')
plt.xlabel('')
plt.legend(title='')
plt.facet_wrap('HRD_index')
statannot.add_stat_annotation(data=df_hrdIndex, x='HRD', y='score', 
                              hue='HRD_index', box_pairs=[('HR-proficient', 'HRD')],
                              test='t-test_ind', text_format='p={:.3e}',
                              loc='outside', verbose=2)
plt.savefig('Figures/Supp_HRDvsHRDscore_ind.pdf')
plt.close()

# CX3 Copy Number Signature Exposure
plt.figure(figsize=(4, 4))
sns.boxplot(data=resultsHRD, x='HRD', y='CX3')
plt.xlabel('')
statannot.add_stat_annotation(data=resultsHRD, x='HRD', y='CX3', 
                              box_pairs=[('HR-proficient', 'HRD')],
                              test='t-test_ind', text_format='p={:.3e}',
                              loc='outside', verbose=2)
plt.savefig('Figures/Figure2/HRD_vs_CX3Exposure.pdf')
plt.close()

# Quiescence Score Comparison
plt.figure(figsize=(4, 4))
sns.boxplot(data=resultsHRD, x='HRD', y='ProliferativeCapacity')
plt.xlabel('')
statannot.add_stat_annotation(data=resultsHRD, x='HRD', y='ProliferativeCapacity', 
                              box_pairs=[('HR-proficient', 'HRD')],
                              test='t-test_ind', text_format='p={:.3e}',
                              loc='outside', verbose=2)
plt.savefig('Figures/Figure2/HRD_vs_Quiescence.pdf')
plt.close()

# Breast Cancer Subtypes
res_BRCA = resultsHRD[['Patient', 'HRD', 'er_status_by_ihc', 'her2_status_by_ihc', 'pr_status_by_ihc']]

res_BRCA['BRCA_subtype'] = np.where(res_BRCA['er_status_by_ihc'] == 'Positive', 'ER+',
                           np.where(res_BRCA['her2_status_by_ihc'] == 'Positive', 'HER2+',
                           np.where((res_BRCA['er_status_by_ihc'] == 'Positive') & 
                                    (res_BRCA['her2_status_by_ihc'] == 'Positive'), 'ER+HER2+',
                           np.where((res_BRCA['er_status_by_ihc'] == 'Negative') & 
                                    (res_BRCA['her2_status_by_ihc'] == 'Negative') & 
                                    (res_BRCA['pr_status_by_ihc'] == 'Negative'), 'TNBC', np.nan))))

res_BRCA_plot = res_BRCA.dropna(subset=['BRCA_subtype']).groupby(['HRD', 'BRCA_subtype']).size().reset_index(name='n')

plt.figure(figsize=(4, 4))
sns.barplot(data=res_BRCA_plot, x='BRCA_subtype', y='n', hue='HRD')
plt.xlabel('')
plt.ylabel('% Samples')
plt.legend(title='')
plt.savefig('Figures/Figure2/HRDvsBRCAsubtype.pdf')
plt.close()

plt.figure(figsize=(6, 3))
sns.barplot(