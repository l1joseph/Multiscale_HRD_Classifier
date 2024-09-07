import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# Load and process data
Z_tumor_training = pd.read_pickle('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.pkl')

samples = Z_tumor_training.index.str[:12]
Z_tumor_training = np.log2(Z_tumor_training + 1)
Z_tumor_training = stats.zscore(Z_tumor_training, axis=0)
Z_tumor_training = pd.DataFrame(Z_tumor_training, index=samples)

# Load annotation data
ann_tcga = pd.read_pickle('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl')
ann_tcga = ann_tcga.set_index('Patient')[['BRCA_status', 'HRD']]
ann_tcga_train = ann_tcga.loc[Z_tumor_training.index]

# Load signature and subset expression data
signature_centroid_list = pd.read_pickle('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl')
sig = signature_centroid_list['ElasticNet_alpha0.25']

Z_tumor_training = Z_tumor_training[sig.index]

# Define color palette
ann_cols = {
    'BRCA_status': {'BRCA1': 'blue', 'BRCA2': 'red', 'none': 'white'},
    'HRD': {'HRD': sns.color_palette("GrandBudapest1")[1], 'HR-proficient': sns.color_palette("GrandBudapest1")[0]}
}

# Set color range
palette_length = 100
my_color = sns.color_palette("coolwarm", palette_length)
my_breaks = np.concatenate([
    np.linspace(Z_tumor_training.min().min(), -3, num=palette_length//4 + 1),
    np.linspace(-3, 0, num=palette_length//4)[1:],
    np.linspace(0, 3, num=palette_length//4)[1:],
    np.linspace(3, Z_tumor_training.max().max(), num=palette_length//4)[1:]
])

# Create the heatmap
plt.figure(figsize=(20, 15))
g = sns.clustermap(Z_tumor_training.T, 
                   cmap=my_color,
                   vmin=my_breaks.min(), vmax=my_breaks.max(),
                   col_colors=ann_tcga_train,
                   col_cluster=False,
                   row_cluster=True,
                   xticklabels=False, yticklabels=False,
                   cbar_kws={"ticks": my_breaks[::len(my_breaks)//5]},
                   method='average')

# Add color bar for annotations
for label in ann_cols:
    for status, color in ann_cols[label].items():
        g.ax_col_dendrogram.bar(0, 0, color=color, label=status, linewidth=0)
g.ax_col_dendrogram.legend(title='', loc='center', ncol=2)

plt.title('TCGA-BRCA Training Cohort Heatmap')
plt.tight_layout()
plt.savefig('~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_Heatmap.pdf')
plt.close()