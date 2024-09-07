import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pyvis.network import Network

# Set working directory
os.chdir('/path/to/Projects/HRD_TranscriptionalSignature/CellphoneDB/Results/')

# Load data
bassez = pd.read_csv('Bassez2021/significant_means.txt', sep='\t')
bassez.columns = [col.replace('HR.proficient', 'HR_proficient') for col in bassez.columns]

# Select relevant columns
bassez = bassez.iloc[:, [1] + list(range(12, bassez.shape[1]))]

# Reshape data
bassez_long = bassez.melt(id_vars='interacting_pair', var_name='cell_int', value_name='interaction')
bassez_long['SOURCE'] = bassez_long['cell_int'].apply(lambda x: x.split('.')[0])
bassez_long['TARGET'] = bassez_long['cell_int'].apply(lambda x: x.split('.')[1])

# Summarize interactions
bassez_summary = bassez_long.groupby(['SOURCE', 'TARGET']).agg({'interaction': lambda x: sum(pd.notnull(x))}).reset_index()

# Create heatmap data
bassez_heatmap = bassez_summary.pivot(index='SOURCE', columns='TARGET', values='interaction')

# Define cell type order
cellType_order = ['Fibroblast', 'Endothelial_cell', 'Myeloid_cell', 'pDC', 'Mast_cell',
                  'T_cell', 'B_cell', 'Cancer_cellHR_proficient', 'Cancer_cellHRD']

# Reorder heatmap data
bassez_heatmap = bassez_heatmap.reindex(index=cellType_order, columns=cellType_order)

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(bassez_heatmap, cmap='RdBu_r', annot=True, fmt='d')
plt.title('Interaction Counts')
plt.tight_layout()
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Bassez_Heatmap.pdf')
plt.close()

# Prepare data for Sankey diagram
bassez_sankey = bassez_summary[(~bassez_summary['SOURCE'].str.contains('Cancer')) & 
                               (bassez_summary['TARGET'].str.contains('Cancer'))]

# Create nodes dataframe
nodes = pd.DataFrame({'name': pd.concat([bassez_sankey['SOURCE'], bassez_sankey['TARGET']]).unique()})
nodes['id'] = range(len(nodes))

# Add source and target ids to sankey dataframe
bassez_sankey = bassez_sankey.merge(nodes, left_on='SOURCE', right_on='name').rename(columns={'id': 'source'})
bassez_sankey = bassez_sankey.merge(nodes, left_on='TARGET', right_on='name').rename(columns={'id': 'target'})

# Create Sankey diagram using pyvis
net = Network(height='750px', width='100%', directed=True)
net.from_pandas(bassez_sankey, source='source', target='target', value='interaction')
net.show('sankey_diagram.html')

# Plot bar charts
bassez_sankey['TARGET'] = pd.Categorical(bassez_sankey['TARGET'], 
                                         categories=['Cancer_cellHR_proficient', 'Cancer_cellHRD'])

plt.figure(figsize=(12, 6))
g = sns.catplot(data=bassez_sankey, x='TARGET', y='interaction', hue='TARGET',
                col='SOURCE', kind='bar', height=4, aspect=.7, col_wrap=4)
g.set_axis_labels('', 'Interaction Count')
g.set_titles('{col_name}')
plt.tight_layout()
plt.savefig('/path/to/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Bassez_countBarplot.pdf')
plt.close()


# notes:

# We can change the sankey plot to be static if deemed necessary.