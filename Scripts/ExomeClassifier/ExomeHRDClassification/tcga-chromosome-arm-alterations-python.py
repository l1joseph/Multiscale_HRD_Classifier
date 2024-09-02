import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# Set working directory
os.chdir('~/Projects/HRD_MutationalSignature/')

# Load arm-level CNA data and organize
cna = pd.read_csv('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_armlevel_cna.txt', sep='\t')
cna = cna.set_index('NAME').iloc[:, 3:].T
cna['Patient'] = cna.index.str.split('-').str[:3].str.join('-')

# Load HRD/HR-proficiency labels and merge with CNA data
ann_tcga = pd.read_pickle('Results/TCGA_HRDclassification_BRCAannotation.pkl')
ann_tcga = ann_tcga[ann_tcga.BRCA_status.notna()]
ann_tcga = ann_tcga[ann_tcga.HRD == 'HRD']
ann_tcga['group'] = np.where(ann_tcga.BRCA_status == 'none', 'BRCA+', 'BRCA-defective')

df_ann = ann_tcga[['Patient', 'group']]

df = pd.merge(df_ann, cna, on='Patient')

# Initialize data frames to track gain/loss enrichments
arms = df.columns[2:]
fishers_gain = pd.DataFrame({'Arm': arms, 'Estimate': np.nan, 'pVal': np.nan})
fishers_loss = fishers_gain.copy()

# Conduct Fisher's exact tests
for arm in arms:
    df_i = pd.DataFrame({
        'group': df['group'],
        'Gain': df[arm] == 'Gain',
        'Loss': df[arm] == 'Loss'
    }).dropna()
    
    df_i['group'] = pd.Categorical(df_i['group'], categories=['BRCA-defective', 'BRCA+'])
    df_i['Gain'] = pd.Categorical(df_i['Gain'], categories=[False, True])
    df_i['Loss'] = pd.Categorical(df_i['Loss'], categories=[False, True])
    
    # Gains
    table_gain = pd.crosstab(df_i['group'], df_i['Gain'])
    odds_ratio, p_value = fisher_exact(table_gain)
    fishers_gain.loc[fishers_gain['Arm'] == arm, ['Estimate', 'pVal']] = [odds_ratio, p_value]
    
    # Losses
    table_loss = pd.crosstab(df_i['group'], df_i['Loss'])
    odds_ratio, p_value = fisher_exact(table_loss)
    fishers_loss.loc[fishers_loss['Arm'] == arm, ['Estimate', 'pVal']] = [odds_ratio, p_value]

# Process results
for df in [fishers_gain, fishers_loss]:
    df['padj'] = multipletests(df['pVal'], method='fdr_bh')[1]
    df['Estimate'] = np.where(df['Estimate'] == np.inf, df['Estimate'].replace(np.inf, np.nan).max(), df['Estimate'])
    df['Estimate'] = np.where(df['Estimate'] == 0, df['Estimate'].replace(0, np.nan).min(), df['Estimate'])
    df['l2fc'] = np.log2(df['Estimate'])
    df['logp'] = -np.log10(df['padj'])
    df['label'] = np.where(df['padj'] <= 0.05, df['Arm'], '')

# Merge log-fold changes
df_fishers = pd.merge(fishers_gain[['Arm', 'l2fc', 'label']], 
                      fishers_loss[['Arm', 'l2fc', 'label']], 
                      on='Arm', suffixes=('_GAIN', '_LOSS'))

df_fishers['Label'] = 'none'
df_fishers.loc[df_fishers['label_GAIN'] != '', 'Label'] = 'GAIN'
df_fishers.loc[df_fishers['label_LOSS'] != '', 'Label'] = 'LOSS'
df_fishers.loc[(df_fishers['label_GAIN'] != '') & (df_fishers['label_LOSS'] != ''), 'Label'] = 'GAIN+LOSS'
df_fishers['Arm_Label'] = np.where(df_fishers['Label'] != 'none', df_fishers['Arm'], '')

# Plot
plt.figure(figsize=(8, 8))
sns.scatterplot(data=df_fishers, x='l2fc_GAIN', y='l2fc_LOSS', hue='Label', style='Label')
for _, row in df_fishers[df_fishers['Arm_Label'] != ''].iterrows():
    plt.annotate(row['Arm'], (row['l2fc_GAIN'], row['l2fc_LOSS']))
plt.axhline(y=0, color='k', linestyle='--')
plt.axvline(x=0, color='k', linestyle='--')
plt.xlabel('log2(Fold Change) - Gain')
plt.ylabel('log2(Fold Change) - Loss')
plt.title('Chromosome Arm Alterations in HRD vs HR-proficient Samples')
plt.legend(title='', loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig('Figures/Supp_ChrArmEnrich_GainVsLoss_HRDonly.pdf')
plt.close()
