import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Set working directory
os.chdir('~/Projects/HRD_MutationalSignature/')

# Load arm-level CNA data
cna = pd.read_csv('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_armlevel_cna.txt', sep='\t')
cna.set_index('NAME', inplace=True)
cna = cna.iloc[:, 3:].T
cna['Patient'] = cna.index.str.split('.').str[:3].str.join('-')

# Load HRD/HR-proficiency labels
results_tcga_df = pd.read_pickle('Results/TCGA_HRD_resultsSummary.pkl')
df_ann = pd.DataFrame({
    'Patient': results_tcga_df['Patient'],
    'HRD': results_tcga_df['HRD']
})

df = pd.merge(df_ann, cna, on='Patient')

# Initialize dataframes to track gain/loss enrichments
arms = df.columns[2:]
fishers_gain = pd.DataFrame({'Arm': arms, 'Estimate': np.nan, 'pVal': np.nan})
fishers_loss = fishers_gain.copy()

# Perform Fisher's exact tests
for arm in arms:
    df_i = pd.DataFrame({
        'HRD': df['HRD'],
        'Gain': df[arm] == 'Gain',
        'Loss': df[arm] == 'Loss'
    }).dropna()
    
    df_i['HRD'] = pd.Categorical(df_i['HRD'], categories=['HR-proficient', 'HRD'])
    
    # Gains
    table_gain = pd.crosstab(df_i['HRD'], df_i['Gain'])
    _, p_value_gain = fisher_exact(table_gain)
    fishers_gain.loc[fishers_gain['Arm'] == arm, 'Estimate'] = table_gain.iloc[1, 1] * table_gain.iloc[0, 0] / (table_gain.iloc[0, 1] * table_gain.iloc[1, 0])
    fishers_gain.loc[fishers_gain['Arm'] == arm, 'pVal'] = p_value_gain
    
    # Losses
    table_loss = pd.crosstab(df_i['HRD'], df_i['Loss'])
    _, p_value_loss = fisher_exact(table_loss)
    fishers_loss.loc[fishers_loss['Arm'] == arm, 'Estimate'] = table_loss.iloc[1, 1] * table_loss.iloc[0, 0] / (table_loss.iloc[0, 1] * table_loss.iloc[1, 0])
    fishers_loss.loc[fishers_loss['Arm'] == arm, 'pVal'] = p_value_loss

# Process results
for df in [fishers_gain, fishers_loss]:
    df['padj'] = multipletests(df['pVal'], method='fdr_bh')[1]
    df['Estimate'] = df['Estimate'].replace([np.inf, 0], [df['Estimate'].max(), df['Estimate'].min()])
    df['l2fc'] = np.log2(df['Estimate'])
    df['logp'] = -np.log10(df['padj'])
    df['label'] = df['Arm']
    df.loc[df['padj'] > 0.01, 'label'] = np.nan

# Merge results
df_fishers = pd.merge(
    fishers_gain[['Arm', 'l2fc', 'label']].rename(columns={'l2fc': 'l2fc_GAIN', 'label': 'label_GAIN'}),
    fishers_loss[['Arm', 'l2fc', 'label']].rename(columns={'l2fc': 'l2fc_LOSS', 'label': 'label_LOSS'}),
    on='Arm'
)

df_fishers['Label'] = 'none'
df_fishers.loc[df_fishers['label_GAIN'].notna(), 'Label'] = 'GAIN'
df_fishers.loc[df_fishers['label_LOSS'].notna(), 'Label'] = 'LOSS'
df_fishers.loc[(df_fishers['label_GAIN'].notna()) & (df_fishers['label_LOSS'].notna()), 'Label'] = 'GAIN+LOSS'
df_fishers['Arm_Label'] = np.where(df_fishers['Label'] != 'none', df_fishers['Arm'], np.nan)

# Plot results
plt.figure(figsize=(8, 8))
sns.scatterplot(data=df_fishers, x='l2fc_GAIN', y='l2fc_LOSS', hue='Label', palette='colorblind')
for _, row in df_fishers.iterrows():
    if row['Arm_Label'] is not None:
        plt.annotate(row['Arm_Label'], (row['l2fc_GAIN'], row['l2fc_LOSS']))
plt.axhline(y=0, color='gray', linestyle='--')
plt.axvline(x=0, color='gray', linestyle='--')
plt.xlabel('log2(Fold Change) - Gain')
plt.ylabel('log2(Fold Change) - Loss')
plt.legend(title='', loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4)
plt.savefig('Figures/Figure2/ChrArmEnrich_GainVsLoss.pdf', bbox_inches='tight')
plt.close()
