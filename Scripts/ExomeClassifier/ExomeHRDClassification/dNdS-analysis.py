import os
import pandas as pd
import numpy as np
from TCGAbiolinks import GDCquery, GDCdownload, GDCprepare
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
import seaborn as sns

# Set working directory
os.chdir('~/Projects/HRD_MutationalSignature/')

# Load required R packages
pandas2ri.activate()
maftools = importr('maftools')
dndscv = importr('dndscv')

# Load dN/dS references
r('load("~/Data/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda")')
RefCDS = r('RefCDS')

# Load TCGA data
os.chdir('~/Data/TCGA')
query = GDCquery(
    project='TCGA-BRCA',
    data_category='Simple Nucleotide Variation',
    data_type='Masked Somatic Mutation'
)
# GDCdownload(query)
mutations = GDCprepare(query=query)

data = maftools.read_maf(maf=mutations, isTCGA=True)
data = pd.concat([data.data, data.maf_silent])
dat_small = data[['Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']]
dat_small.columns = ['sampleID', 'chr', 'pos', 'ref', 'mut']

# Load HRD/BRCA groups
os.chdir('~/Projects/HRD_MutationalSignature/')
r('load("Results/TCGA_HRDclassification_BRCAannotation.Rdata")')
ann_tcga = r('ann_tcga')
ann_tcga = ann_tcga[~ann_tcga.BRCA_status.isna()]
ann_tcga = ann_tcga[~ann_tcga.BRCA_status.isin(['RAD51C', 'PALB2'])]
dat_small = pd.merge(dat_small, ann_tcga[['Patient', 'HRD', 'BRCA_status']], 
                     left_on='sampleID', right_on='Patient')

# Separate data into different groups
data_HRD = dat_small[dat_small.HRD == 'HRD'].iloc[:, :5]
data_HRproficient = dat_small[dat_small.HRD == 'HR-proficient'].iloc[:, :5]
data_BRCA1 = dat_small[dat_small.BRCA_status == 'BRCA1'].iloc[:, :5]
data_BRCA2 = dat_small[dat_small.BRCA_status == 'BRCA2'].iloc[:, :5]
data_RAD51C = dat_small[dat_small.BRCA_status == 'RAD51C'].iloc[:, :5]
data_HRDBRCApos = dat_small[(dat_small.HRD == 'HRD') & (dat_small.BRCA_status == 'none')].iloc[:, :5]

# Run dN/dS analysis on each group
def run_dndscv(data):
    result = dndscv.dndscv(data, cv=None, refdb=RefCDS)
    return result.rx2('sel_cv')

sel_cvBRCA1 = run_dndscv(data_BRCA1)
sel_cvBRCA2 = run_dndscv(data_BRCA2)
sel_cvRAD51C = run_dndscv(data_RAD51C)
sel_cvHRDBRCApos = run_dndscv(data_HRDBRCApos)
sel_cvHRD = run_dndscv(data_HRD)
sel_cvHRprof = run_dndscv(data_HRproficient)

# Plot results
def plot_dndscv(sel_cv, title):
    df = pd.DataFrame({
        'gene_name': sel_cv.rx2('gene_name'),
        'wind_cv': sel_cv.rx2('wind_cv'),
        'qind_cv': sel_cv.rx2('qind_cv')
    })
    df = df[df.qind_cv < 0.1]
    plt.figure(figsize=(10, 6))
    plt.scatter(np.log2(df.wind_cv), -np.log10(df.qind_cv))
    for _, row in df.iterrows():
        plt.annotate(row.gene_name, (np.log2(row.wind_cv), -np.log10(row.qind_cv)))
    plt.xlabel('log2(dN/dS)')
    plt.ylabel('-log10(q-value)')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(f'Figures/Figure2/dNdS_{title.replace(" ", "_")}.pdf')
    plt.close()

plot_dndscv(sel_cvHRD, 'HRD')
plot_dndscv(sel_cvHRprof, 'HR-proficient')
plot_dndscv(sel_cvBRCA1, 'BRCA1')
plot_dndscv(sel_cvBRCA2, 'BRCA2')
plot_dndscv(sel_cvRAD51C, 'RAD51C')
plot_dndscv(sel_cvHRDBRCApos, 'HRD BRCA+')
