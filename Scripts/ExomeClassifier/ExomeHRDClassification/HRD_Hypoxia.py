import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Set working directory

# Load required data
ann_tcga = pd.read_pickle('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl')
hypoxia = pd.read_csv('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_clinical_supp_hypoxia.txt', sep='\t')

ann_tcga = ann_tcga[ann_tcga['BRCA_status'].notna() & 
                    ~ann_tcga['BRCA_status'].isin(['PALB2', 'RAD51C'])]
ann_tcga = ann_tcga[ann_tcga['ER_status'].isin(['Negative', 'Positive'])]

df = pd.merge(ann_tcga[['Patient', 'HRD', 'BRCA_status', 'ER_status']], 
              hypoxia[['PATIENT_ID', 'BUFFA_HYPOXIA_SCORE']], 
              left_on='Patient', right_on='PATIENT_ID')
df = df.rename(columns={'BUFFA_HYPOXIA_SCORE': 'Hypoxia_Buffa'})
df['ER_status'] = df['ER_status'].map({'Negative': 'ER-', 'Positive': 'ER+'})
df['ER_status'] = pd.Categorical(df['ER_status'], categories=['ER+', 'ER-'])

df['Group'] = df['BRCA_status']
df.loc[(df['BRCA_status'] == 'none') & (df['HRD'] == 'HRD'), 'Group'] = 'HRD_BRCA+'
df.loc[df['Group'] == 'none', 'Group'] = 'HR-proficient'
df['Group'] = pd.Categorical(df['Group'], categories=['HR-proficient', 'HRD_BRCA+', 'BRCA1', 'BRCA2'])

# Plotting
plt.figure(figsize=(8, 4))
sns.violinplot(data=df, x='Group', y='Hypoxia_Buffa', hue='ER_status', split=True)
plt.title('Hypoxia Score (Buffa)')
plt.xlabel('')
plt.legend(title='')

# Add statistical comparisons
comps = [('HR-proficient', 'HRD_BRCA+'), ('HRD_BRCA+', 'BRCA2'), ('HRD_BRCA+', 'BRCA1')]
y_max = df['Hypoxia_Buffa'].max()
y_increment = (y_max - df['Hypoxia_Buffa'].min()) * 0.05

for i, comp in enumerate(comps):
    group1 = df[df['Group'] == comp[0]]['Hypoxia_Buffa']
    group2 = df[df['Group'] == comp[1]]['Hypoxia_Buffa']
    t_stat, p_val = stats.ttest_ind(group1, group2)
    y_pos = y_max + (i+1) * y_increment
    plt.plot([df['Group'].cat.categories.get_loc(comp[0]), df['Group'].cat.categories.get_loc(comp[1])], 
             [y_pos, y_pos], 'k-')
    plt.text((df['Group'].cat.categories.get_loc(comp[0]) + df['Group'].cat.categories.get_loc(comp[1]))/2, 
             y_pos, f'p = {p_val:.3f}', ha='center', va='bottom')

plt.savefig('~/Projects/HRD_MutationalSignature/Figures/Figure2/HRD_Hypoxia.pdf', bbox_inches='tight')
plt.close()

# ER+ only plot
plt.figure(figsize=(6, 4))
sns.boxplot(data=df[df['ER_status'] == 'ER+'], x='Group', y='Hypoxia_Buffa')
plt.title('Hypoxia Score (Buffa) - ER+ only')
plt.xlabel('')

# Add statistical comparisons for ER+ only
for i, comp in enumerate(comps):
    group1 = df[(df['Group'] == comp[0]) & (df['ER_status'] == 'ER+')]['Hypoxia_Buffa']
    group2 = df[(df['Group'] == comp[1]) & (df['ER_status'] == 'ER+')]['Hypoxia_Buffa']
    t_stat, p_val = stats.ttest_ind(group1, group2)
    y_pos = y_max + (i+1) * y_increment
    plt.plot([df['Group'].cat.categories.get_loc(comp[0]), df['Group'].cat.categories.get_loc(comp[1])], 
             [y_pos, y_pos], 'k-')
    plt.text((df['Group'].cat.categories.get_loc(comp[0]) + df['Group'].cat.categories.get_loc(comp[1]))/2, 
             y_pos, f'p = {p_val:.3f}', ha='center', va='bottom')

plt.savefig('~/Projects/HRD_MutationalSignature/Figures/Figure2/HRD_Hypoxia_ERpositive.pdf', bbox_inches='tight')
plt.close()

# Linear model
from statsmodels.formula.api import ols

model = ols('Hypoxia_Buffa ~ C(ER_status) + C(BRCA_status) * C(HRD)', data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)