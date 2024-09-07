import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from TCGAutils import TCGAutils  # Assuming a Python equivalent exists; if not, we'll need to implement custom functions
import anndata
from decoupler import decouple, get_progeny  # Assuming a Python package for decoupleR exists

# Set working directory
os.chdir('~/Data/TCGA')

# Load TCGA data
query = TCGAutils.GDCquery(
    project='TCGA-BRCA',
    data_category='Transcriptome Profiling',
    data_type='Gene Expression Quantification',
    workflow_type='STAR - Counts'
)
# TCGAutils.GDCdownload(query)
data = TCGAutils.GDCprepare(query=query)
data = data[data.sample_type == 'Primary Tumor']
data = data[~data.patient.duplicated()]

# Match with HRD classification
ann_tcga = pd.read_pickle('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl')
ann_tcga['HRD'] = ann_tcga['HRD_prob'].apply(lambda x: 'HRD' if x >= 0.79 else 'HR-proficient')
ann_tcga = ann_tcga[ann_tcga['ER_status'].notna()]

patients_intersect = list(set(ann_tcga['Patient']) & set(data.patient))
ann_tcga = ann_tcga[ann_tcga['Patient'].isin(patients_intersect)]
data = data[data.patient.isin(patients_intersect)]

data.obs['HRD_status'] = ann_tcga['HRD'].values

# Construct a DESeqDataSet data object
# Note: This is a simplified version of constructing a DESeqDataSet object
# For a more robust analysis, we should consider using a package like PyDESeq2.
# I also have my own DeSeq2 equivilent package that I've implemented,
# but it needs to be fleshed out further before its ready to use in this project.
counts = data.to_df()

# Gene count filtering
counts = counts.loc[:, counts.sum() > 10]

# Normalization (simplified)
size_factors = counts.sum(axis=0) / counts.sum(axis=0).median()
counts_normalized = counts / size_factors

# DIFFERENTIAL GENE EXPRESSION ANALYSIS
# Note: This is a simplified version of differential expression analysis
# For a more robust analysis, consider using a package like PyDESeq2

# Calculate mean expression for each group
hrd_mean = counts_normalized.loc[:, data.obs['HRD_status'] == 'HRD'].mean(axis=1)
hr_prof_mean = counts_normalized.loc[:, data.obs['HRD_status'] == 'HR-proficient'].mean(axis=1)

# Calculate log2 fold change
log2_fold_change = np.log2(hrd_mean / hr_prof_mean)

# Perform t-test
t_stat, p_val = stats.ttest_ind(
    counts_normalized.loc[:, data.obs['HRD_status'] == 'HRD'],
    counts_normalized.loc[:, data.obs['HRD_status'] == 'HR-proficient'],
    axis=1
)

dds_statVals = pd.DataFrame({
    'EnsemblID': counts.index,
    'stat': t_stat,
    'ID': data.var['gene_name']
})
dds_statVals = dds_statVals.drop_duplicates(subset='ID')
dds_statVals = dds_statVals.set_index('ID')

deg = dds_statVals[['stat']].T

counts_logNorm = np.log2(counts + 1)
counts_logNorm.columns = [col[:12] for col in counts_logNorm.columns]

design = ann_tcga[['Patient', 'HRD_status']]
design.columns = ['sample', 'condition']

# Run progeny
net = get_progeny(organism='human', top=100)

# Run mlm
contrast_acts = decouple.run_mlm(
    mat=deg,
    net=net,
    source='source',
    target='target',
    weight='weight',
    minsize=5
)

# Plot
g_progeny = sns.barplot(data=contrast_acts, x='source', y='score', 
                        order=contrast_acts.sort_values('score', ascending=False)['source'])
g_progeny.set_xticklabels(g_progeny.get_xticklabels(), rotation=45, ha='right')
plt.title('Enrichment in HRD Samples')
plt.xlabel('Pathways')
plt.ylabel('Enrichment Score')
plt.tight_layout()
plt.savefig('~/Projects/HRD_MutationalSignature/Figures/Figure2/HRD_decoupleR.pdf', dpi=300)
plt.close()



#notes:

# Need to find a python equivalent for decoupleR, TCGAutils, DeSeq2(PyDeq2).