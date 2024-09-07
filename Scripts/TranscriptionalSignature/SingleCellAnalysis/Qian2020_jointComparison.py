import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set working directory
os.chdir('/path/to/Data/scRNASeq/Qian2020/')

# Load libraries and data
def load_rdata(file_path):
    # This is a placeholder function. 
    # a proper R data loading mechanism, such as using the pyreadr 
    pass

expr_data_qian2020 = load_rdata('exprData_Qian2020.Rdata')

# Load expression and metadata and subset for cancer cells
meta_qian = pd.read_csv('2103-Breastcancer_metadata.csv', header=0)
meta_qian = meta_qian[meta_qian['Cell'].isin(expr_data_qian2020.columns)]
expr_cancer = expr_data_qian2020.loc[:, meta_qian['CellType'] == 'Cancer']

# Load signature centroids
centroids = load_rdata('/path/to/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
signature_centroid_list = centroids['signature.centroid.list']

# Uncomment these lines if you need to load alternative signatures
# alternative_centroids = load_rdata('/path/to/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.Rdata')
# signature_alternative_centroid_list = alternative_centroids['signature_alternative.centroid.list']
# signature_alternative_centroid_list = {f'Alternative_Sig_{k}': v for k, v in signature_alternative_centroid_list.items()}
# signature_centroid_list.update(signature_alternative_centroid_list)

stats_df = pd.DataFrame(columns=['Model', 'prop_nonZeroCells', 'median_GenesExpressedPerCell', 'mean_GenesExpressedPerCell'])

for model, sig in signature_centroid_list.items():
    print(model)
    
    genes_intersect = list(set(sig.index) & set(expr_cancer.index))
    
    sig = sig.loc[genes_intersect]
    expr_cancer_i = expr_cancer.loc[genes_intersect]
    
    # Collate relevant results
    stats_df = stats_df.append({
        'Model': model,
        'prop_nonZeroCells': (expr_cancer_i.sum() > 0).mean(),
        'median_GenesExpressedPerCell': (expr_cancer_i > 0).sum().median(),
        'mean_GenesExpressedPerCell': (expr_cancer_i > 0).sum().mean()
    }, ignore_index=True)

# Plotting
stats_df_plot = pd.melt(stats_df, id_vars=['Model'], var_name='Measure', value_name='Value')

g_dropout = sns.catplot(
    data=stats_df_plot, 
    x='Model', y='Value', 
    col='Measure', 
    kind='bar', 
    height=4, aspect=1,
    col_wrap=3
)

g_dropout.set_xticklabels(rotation=90)
g_dropout.fig.suptitle('Signature Comparison', y=1.02)
plt.tight_layout()
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Supp_QianSignatureComparison.pdf', 
            bbox_inches='tight')
plt.close()