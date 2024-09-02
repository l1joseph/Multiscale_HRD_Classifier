import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def perform_cna_enrichment_analysis(hrd_only=False):
    # Set working directory
    os.chdir('~/Projects/HRD_MutationalSignature/')

    # Load Gene-level CNA data
    cna = pd.read_csv('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_cna.txt', sep='\t')
    cna = cna.drop_duplicates(subset='Hugo_Symbol')

    # Extract cancer genes
    cgc_data = pd.read_csv('~/Data/Census_allMon Jul  3 15_56_40 2023.csv')
    cgc_data_genes = cgc_data['Gene.Symbol'].tolist()
    cna = cna[cna['Hugo_Symbol'].isin(cgc_data_genes)]

    # Sort dataframe
    cna = cna.set_index('Hugo_Symbol').iloc[:, 2:].T
    cna['Patient'] = cna.index.str.split('.').str[:3].str.join('-')

    # Load HRD/HR-proficiency labels
    ann_tcga = pd.read_pickle('Results/TCGA_HRDclassification_BRCAannotation.pkl')
    ann_tcga = ann_tcga[~ann_tcga['BRCA_status'].isna()]

    if hrd_only:
        ann_tcga = ann_tcga[ann_tcga['HRD'] == 'HRD']
        ann_tcga['group'] = ann_tcga['BRCA_status'].apply(lambda x: 'BRCA-defective' if x != 'none' else 'BRCA+')
    else:
        ann_tcga['group'] = ann_tcga['HRD']

    df_ann = pd.DataFrame({
        'Patient': ann_tcga['Patient'],
        'group': ann_tcga['group']
    })

    df = df_ann.merge(cna, on='Patient')

    # Initialize dataframes to track gain/loss enrichments
    genes = df.columns[2:]
    fishers_gain = pd.DataFrame({'Gene': genes, 'Estimate': np.nan, 'pVal': np.nan})
    fishers_loss = fishers_gain.copy()

    # Perform Fisher's exact tests
    for gene in genes:
        df_i = pd.DataFrame({
            'group': df['group'],
            'Gain': df[gene] > 0,
            'Loss': df[gene] < 0
        }).dropna()
        
        df_i['group'] = pd.Categorical(df_i['group'], categories=['BRCA-defective', 'BRCA+'] if hrd_only else ['HR-proficient', 'HRD'])
        
        # Gains
        table_gain = pd.crosstab(df_i['group'], df_i['Gain'])
        _, p_value_gain = fisher_exact(table_gain)
        fishers_gain.loc[fishers_gain['Gene'] == gene, 'Estimate'] = table_gain.iloc[1, 1] * table_gain.iloc[0, 0] / (table_gain.iloc[0, 1] * table_gain.iloc[1, 0])
        fishers_gain.loc[fishers_gain['Gene'] == gene, 'pVal'] = p_value_gain
        
        # Losses
        table_loss = pd.crosstab(df_i['group'], df_i['Loss'])
        _, p_value_loss = fisher_exact(table_loss)
        fishers_loss.loc[fishers_loss['Gene'] == gene, 'Estimate'] = table_loss.iloc[1, 1] * table_loss.iloc[0, 0] / (table_loss.iloc[0, 1] * table_loss.iloc[1, 0])
        fishers_loss.loc[fishers_loss['Gene'] == gene, 'pVal'] = p_value_loss

    # Process results
    for df in [fishers_gain, fishers_loss]:
        df['padj'] = multipletests(df['pVal'], method='fdr_bh')[1]
        df['Estimate'] = df['Estimate'].replace([np.inf, 0], [df['Estimate'].max(), df['Estimate'].min()])
        df['l2fc'] = np.log2(df['Estimate'])
        df['logp'] = -np.log10(df['padj'])
        df['label'] = df['Gene']
        df.loc[df['padj'] > 0.05, 'label'] = np.nan

    # Merge results
    df_fishers = pd.merge(
        fishers_gain[['Gene', 'l2fc', 'label']].rename(columns={'l2fc': 'l2fc_GAIN', 'label': 'label_GAIN'}),
        fishers_loss[['Gene', 'l2fc', 'label']].rename(columns={'l2fc': 'l2fc_LOSS', 'label': 'label_LOSS'}),
        on='Gene'
    )

    df_fishers['Label'] = 'none'
    df_fishers.loc[df_fishers['label_GAIN'].notna(), 'Label'] = 'GAIN'
    df_fishers.loc[df_fishers['label_LOSS'].notna(), 'Label'] = 'LOSS'
    df_fishers.loc[(df_fishers['label_GAIN'].notna()) & (df_fishers['label_LOSS'].notna()), 'Label'] = 'GAIN+LOSS'
    df_fishers['Gene_Label'] = np.where(df_fishers['Label'] != 'none', df_fishers['Gene'], np.nan)

    # Plot results
    plt.figure(figsize=(8, 8))
    sns.scatterplot(data=df_fishers, x='l2fc_GAIN', y='l2fc_LOSS', hue='Label', palette='colorblind')
    for _, row in df_fishers.iterrows():
        if pd.notna(row['Gene_Label']):
            plt.annotate(row['Gene_Label'], (row['l2fc_GAIN'], row['l2fc_LOSS']))
    plt.axhline(y=0, color='gray', linestyle='--')
    plt.axvline(x=0, color='gray', linestyle='--')
    plt.xlabel('log2(Fold Change) - Gain')
    plt.ylabel('log2(Fold Change) - Loss')
    plt.legend(title='', loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4)
    plt.title('Gene CNA Enrichment: ' + ('HRD only' if hrd_only else 'All samples'))
    plt.savefig(f'Figures/Supp_GeneCnaEnrich_GainVsLoss{"_HRDonly" if hrd_only else ""}.pdf', bbox_inches='tight')
    plt.close()

# Run analysis for both cases
perform_cna_enrichment_analysis(hrd_only=False)
perform_cna_enrichment_analysis(hrd_only=True)
