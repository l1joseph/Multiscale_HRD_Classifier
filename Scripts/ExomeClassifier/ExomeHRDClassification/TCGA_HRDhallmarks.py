import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Load required data
ann_tcga = pd.read_pickle('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl')
tcga_HRDHallmarks = pd.read_csv('~/Data/TCGA/HRD_hallmarks.txt', sep='\t')

resultsHRD = pd.merge(ann_tcga[['Patient', 'Phenotype_Assigned', 'HRD']], 
                      tcga_HRDHallmarks, left_on='Patient', right_on='Patient', how='left')
resultsHRD['HRD'] = pd.Categorical(resultsHRD['HRD'], categories=['HR-proficient', 'HRD'])

# MYC amplification
resultsHRD['MYC_status'] = pd.Categorical(resultsHRD['MYC_status'], 
                                          categories=['Deletion', 'Normal', 'Amplification'])
res_MYC = resultsHRD.groupby(['HRD', 'MYC_status']).size().reset_index(name='n')
res_MYC = res_MYC.dropna()

plt.figure(figsize=(6, 3))
g_myc = sns.barplot(x='HRD', y='n', hue='MYC_status', data=res_MYC, 
                    palette=[sns.color_palette("Zissou1")[3], 'gray90', sns.color_palette("Zissou1")[1]])
plt.xlabel('')
plt.ylabel('% Samples')
plt.legend(title='', loc='upper right')
plt.gca().invert_yaxis()
plt.savefig('Figures/Figure2/HRD_vs_MYCamplification.pdf', bbox_inches='tight')
plt.close()

# POLQ expression
plt.figure(figsize=(4, 4))
sns.boxplot(x='HRD', y='log_POLQ_FPKM', data=resultsHRD, palette=sns.color_palette("GrandBudapest1"))
plt.xlabel('')
statannot.add_stat_annotation(plt.gca(), data=resultsHRD, x='HRD', y='log_POLQ_FPKM', 
                              test='t-test_ind', text_format='star', loc='inside', verbose=2)
plt.savefig('Figures/Figure2/HRD_vs_POLQexpression.pdf', bbox_inches='tight')
plt.close()

# HRD index scores
plt.figure(figsize=(4, 4))
sns.boxplot(x='HRD', y='HRD_index', data=resultsHRD, palette=sns.color_palette("GrandBudapest1"))
plt.xlabel('')
statannot.add_stat_annotation(plt.gca(), data=resultsHRD, x='HRD', y='HRD_index', 
                              test='t-test_ind', text_format='star', loc='inside', verbose=2)
plt.savefig('Figures/Figure2/HRD_vs_HRDscore.pdf', bbox_inches='tight')
plt.close()

# Individual HRD index scores
df_hrdIndex = resultsHRD.melt(id_vars=['HRD'], value_vars=['NtAI', 'LST', 'HRD.LOH'], 
                              var_name='HRD_index', value_name='score')
plt.figure(figsize=(8, 4))
g_hrdIndex_individual = sns.boxplot(x='HRD', y='score', hue='HRD_index', data=df_hrdIndex, 
                                    palette=sns.color_palette("GrandBudapest1"))
plt.xlabel('')
plt.legend(title='')
for i in range(2):
    statannot.add_stat_annotation(plt.gca(), data=df_hrdIndex[df_hrdIndex['HRD_index'] == df_hrdIndex['HRD_index'].unique()[i]], 
                                  x='HRD', y='score', test='t-test_ind', text_format='star', loc='inside', verbose=2)
plt.savefig('Figures/Supp_HRDvsHRDscore_ind.pdf', bbox_inches='tight')
plt.close()

# CX3 Copy Number Signature Exposure
plt.figure(figsize=(4, 4))
sns.boxplot(x='HRD', y='CX3', data=resultsHRD, palette=sns.color_palette("GrandBudapest1"))
plt.xlabel('')
statannot.add_stat_annotation(plt.gca(), data=resultsHRD, x='HRD', y='CX3', 
                              test='t-test_ind', text_format='star', loc='inside', verbose=2)
plt.savefig('Figures/Figure2/HRD_vs_CX3Exposure.pdf', bbox_inches='tight')
plt.close()

# Quiescence Score Comparison
plt.figure(figsize=(4, 4))
sns.boxplot(x='HRD', y='ProliferativeCapacity', data=resultsHRD, palette=sns.color_palette("GrandBudapest1"))
plt.xlabel('')
statannot.add_stat_annotation(plt.gca(), data=resultsHRD, x='HRD', y='ProliferativeCapacity', 
                              test='t-test_ind', text_format='star', loc='inside', verbose=2)
plt.savefig('Figures/Figure2/HRD_vs_Quiescence.pdf', bbox_inches='tight')
plt.close()

# Breast Cancer Subtypes
res_BRCA = resultsHRD[['Patient', 'HRD', 'er_status_by_ihc', 'her2_status_by_ihc', 'pr_status_by_ihc']]

res_BRCA['BRCA_subtype'] = np.nan
res_BRCA.loc[res_BRCA['er_status_by_ihc'] == 'Positive', 'BRCA_subtype'] = 'ER+'
res_BRCA.loc[res_BRCA['her2_status_by_ihc'] == 'Positive', 'BRCA_subtype'] = 'HER2+'
res_BRCA.loc[(res_BRCA['er_status_by_ihc'] == 'Positive') & 
             (res_BRCA['her2_status_by_ihc'] == 'Positive'), 'BRCA_subtype'] = 'ER+HER2+'
res_BRCA.loc[(res_BRCA['er_status_by_ihc'] == 'Negative') & 
             (res_BRCA['her2_status_by_ihc'] == 'Negative') & 
             (res_BRCA['pr_status_by_ihc'] == 'Negative'), 'BRCA_subtype'] = 'TNBC'

res_BRCA_plot = res_BRCA.dropna(subset=['BRCA_subtype']).groupby(['HRD', 'BRCA_subtype']).size().reset_index(name='n')

plt.figure(figsize=(4, 4))
g_resBRCA = sns.barplot(x='BRCA_subtype', y='n', hue='HRD', data=res_BRCA_plot, 
                        palette=sns.color_palette("GrandBudapest1"))
plt.xlabel('')
plt.ylabel('% Samples')
plt.legend(title='', loc='upper right')
plt.savefig('Figures/Figure2/HRDvsBRCAsubtype.pdf', bbox_inches='tight')
plt.close()

plt.figure(figsize=(6, 3))
g_resBRCA2 = sns.barplot(x='HRD', y='n', hue='BRCA_subtype', data=res_BRCA_plot, 
                         palette=sns.color_palette("Royal2")[3::-1])
plt.xlabel('')
plt.ylabel('% Samples')
plt.legend(title='', loc='upper right')
plt.gca().invert_yaxis()
plt.savefig('Figures/Figure2/BRCAsubtypevsHRD.pdf', bbox_inches='tight')
plt.close()


# notes

# using the statannot library, which provides similar functionality to ggpubr's stat_compare_means