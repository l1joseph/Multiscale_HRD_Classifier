import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

# Set working directory

# Load arm-level CNA data, and organize to remove metadata
cna = pd.read_csv('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_armlevel_cna.txt', sep='\t')
cna = cna.set_index('NAME')
cna = cna.drop(['Chromosome', 'Arm'], axis=1)
cna = cna.T
cna['Patient'] = cna.index.map(lambda x: '-'.join(x.split('-')[:3]))

# Load HRD/HR-proficiency labels and merge with CNA data
ann_tcga = pd.read_pickle('Results/TCGA_HRDclassification_BRCAannotation.pkl')
ann_tcga = ann_tcga[ann_tcga['BRCA_status'].notna()]
ann_tcga = ann_tcga[ann_tcga['HRD'] == 'HRD']
ann_tcga['group'] = ann_tcga['BRCA_status'].apply(lambda x: 'BRCA+' if x == 'none' else 'BRCA-defective')

df_ann = ann_tcga[['Patient', 'group']]

df = pd.merge(df_ann, cna, on='Patient')

# Initialize data frames to track gain/loss enrichments
fishers_gain = pd.DataFrame({'Arm': df.columns[2:], 'Estimate': np.nan, 'pVal': np.nan})
fishers_loss = pd.DataFrame({'Arm': df.columns[2:], 'Estimate': np.nan, 'pVal': np.nan})

# For each chromosome arm:
#   Conduct a Fisher's exact test for HRD and HR-proficient samples against:
#     1) Gain vs not-Gain
#     2) Loss vs not-Loss
#   Save odds ratios and p-values in initialized datasets
for i, arm in enumerate(df.columns[2:]):
    arm_i = arm
    
    df_i = pd.DataFrame({
        'group': df['group'],
        'Gain': df[arm_i] == 'Gain',
        'Loss': df[arm_i] == 'Loss'
    })
    df_i = df_i.dropna(subset=['Gain'])
    
    df_i['group'] = pd.Categorical(df_i['group'], categories=['BRCA-defective', 'BRCA+'])
    df_i['Gain'] = pd.Categorical(df_i['Gain'], categories=[False, True])
    df_i['Loss'] = pd.Categorical(df_i['Loss'], categories=[False, True])
    
    # Gains
    table_gain_i = pd.crosstab(df_i['group'], df_i['Gain'])
    fishers_gain_i = stats.fisher_exact(table_gain_i)
    
    fishers_gain.loc[i, 'Estimate'] = fishers_gain_i[0]
    fishers_gain.loc[i, 'pVal'] = fishers_gain_i[1]
    
    # Losses
    table_loss_i = pd.crosstab(df_i['group'], df_i['Loss'])
    fishers_loss_i = stats.fisher_exact(table_loss_i)
    
    fishers_loss.loc[i, 'Estimate'] = fishers_loss_i[0]
    fishers_loss.loc[i, 'pVal'] = fishers_loss_i[1]

# Add information to gain/loss results:
#   - Adjust p-values
#   - If Estimate == 'Inf', set to maximum estimate
#   - If Estimate == 0, set to minimum estimate
#   - log2-normalize estimates and -log10(p-adjust)
#   - Add labels for arms with enrichment significance < .01
#   - 'Volcano plot' displaying results

# Gains
fishers_gain['padj'] = multipletests(fishers_gain['pVal'], method='fdr_bh')[1]
fishers_gain.loc[fishers_gain['Estimate'] == np.inf, 'Estimate'] = fishers_gain['Estimate'].replace(np.inf, np.nan).max()
fishers_gain.loc[fishers_gain['Estimate'] == 0, 'Estimate'] = fishers_gain['Estimate'].replace(0, np.nan).min()
fishers_gain['l2fc'] = np.log2(fishers_gain['Estimate'])
fishers_gain['logp'] = -np.log10(fishers_gain['padj'])
fishers_gain['label'] = fishers_gain['Arm']
fishers_gain.loc[fishers_gain['padj'] > 0.05, 'label'] = np.nan

g_gains = sns.scatterplot(data=fishers_gain, x='l2fc', y='logp')
for _, row in fishers_gain.dropna(subset=['label']).iterrows():
    g_gains.text(row['l2fc'], row['logp'], row['label'], fontsize=8)
g_gains.set_title('Gain vs !Gain')
plt.savefig('Figures/Supplementary/ChrArmEnrich_GainVsNotGain.pdf')
plt.close()

# Losses
fishers_loss['padj'] = multipletests(fishers_loss['pVal'], method='fdr_bh')[1]
fishers_loss.loc[fishers_loss['Estimate'] == np.inf, 'Estimate'] = fishers_loss['Estimate'].replace(np.inf, np.nan).max()
fishers_loss.loc[fishers_loss['Estimate'] == 0, 'Estimate'] = fishers_loss['Estimate'].replace(0, np.nan).min()
fishers_loss['l2fc'] = np.log2(fishers_loss['Estimate'])
fishers_loss['logp'] = -np.log10(fishers_loss['padj'])
fishers_loss['label'] = fishers_loss['Arm']
fishers_loss.loc[fishers_loss['padj'] > 0.05, 'label'] = np.nan

g_losses = sns.scatterplot(data=fishers_loss, x='l2fc', y='logp')
for _, row in fishers_loss.dropna(subset=['label']).iterrows():
    g_losses.text(row['l2fc'], row['logp'], row['label'], fontsize=8)
g_losses.set_title('Loss vs !Loss')
plt.savefig('Figures/Supplementary/ChrArmEnrich_LossVsNotLoss.pdf')
plt.close()

# Merge log-fold changes (+ labels, showing significant arms for each test)
df_fishers = pd.merge(fishers_gain[['Arm', 'l2fc', 'label']], 
                      fishers_loss[['Arm', 'l2fc', 'label']], 
                      on='Arm', suffixes=('_GAIN', '_LOSS'))

# Add GAIN/LOSS/GAIN+LOSS labels for which tests are significant
df_fishers['Label'] = 'none'
df_fishers.loc[df_fishers['label_GAIN'].notna(), 'Label'] = 'GAIN'
df_fishers.loc[df_fishers['label_LOSS'].notna(), 'Label'] = 'LOSS'
df_fishers.loc[(df_fishers['label_GAIN'].notna()) & (df_fishers['label_LOSS'].notna()), 'Label'] = 'GAIN+LOSS'
df_fishers['Arm_Label'] = np.nan
df_fishers.loc[df_fishers['Label'] != 'none', 'Arm_Label'] = df_fishers.loc[df_fishers['Label'] != 'none', 'Arm']

# Plot
plt.figure(figsize=(4, 4))
sns.scatterplot(data=df_fishers, x='l2fc_GAIN', y='l2fc_LOSS', hue='Label', palette=['gray'])
for _, row in df_fishers.dropna(subset=['Arm_Label']).iterrows():
    plt.text(row['l2fc_GAIN'], row['l2fc_LOSS'], row['Arm_Label'], fontsize=8)
plt.axhline(y=0, color='black', linestyle='--')
plt.axvline(x=0, color='black', linestyle='--')
plt.xlabel('log2(Fold Change) - Gain')
plt.ylabel('log2(Fold Change) - Loss')
plt.legend(title='', loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4)
plt.savefig('Figures/Supp_ChrArmEnrich_GainVsLoss_HRDonly.pdf', bbox_inches='tight')
plt.close()