import pandas as pd
import numpy as np
from TCGAbiolinks import TCGAbiolinks
import anndata as ad

# Load TCGA transcriptional data
# Templates will be formed from FPKM-normalized training data
# Therefore, we must also load the training data to obtain sample IDs

Z_tumor_training = pd.read_pickle('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.pkl')
barcodes_training = Z_tumor_training.index

# Load TCGA data
query = TCGAbiolinks.GDCquery(
    project="TCGA-BRCA",
    data_category="Transcriptome Profiling",
    data_type="Gene Expression Quantification",
    workflow_type="STAR - Counts",
    barcode=barcodes_training
)
expr_train = TCGAbiolinks.GDCprepare(query=query)

# Filter and process expression data
expr_train = expr_train[expr_train.sample_type == 'Primary Tumor']
expr_train = expr_train[~expr_train.gene_name.duplicated() & ~expr_train.gene_name.isna()]
genes_tcga = expr_train.gene_name

# Get signatures
signatures_alternative = {}

parpi7 = ['BRCA1', 'MRE11', 'NBN', 'TDG', 'XPA', 'CHEK2', 'MAPKAPK2']
cin70 = ['TPX2','PRC1','FOXM1','CDK1','TGIF2','MCM2','H2AZ1','TOP2A','PCNA','UBE2C',
         'MELK','TRIP13','CEP250','MCM7','RNASEH2A','RAD51AP1','KIF20A','CDC45','MAD2L1','ESPL1',
         'CCNB2','FEN1','TTK','CCT5','RFC4','ATAD2','CKAP5','NUP205','CDC20','CKS2',
         'RRM2','ELAVL1','CCNB1','RRM1','AURKB','MSH6','EZH2','CTPS1','DKC1','OIP5',
         'CDCA8','PTTG1','CEP55','H2AX','CMAS','NCAPH','MCM10','LSM4','NCAPG2','ASF1B',
         'ZWINT','PBK','ZWILCH','CDCA3','ECT2','CDC6','UNG','MTCH2','RAD21','ACTL6A',
         'PDCD2L','SRSF2','HDGF','NXT1','NEK2','DHCR7','AURKA','NDUFAB1','NEMP1','KIF4A']

signatures_alternative['PARPi7'] = parpi7
signatures_alternative['CIN70'] = cin70

# Severson
sev_init = pd.read_excel('../../../Downloads/BRCA1nessSignature_Severson2017.xlsx')
sev = sev_init.iloc[:, 0].tolist()
sev = [gene if gene in genes_tcga.tolist() else replacement for gene, replacement in 
       zip(sev, ['JAML','FAM241B','PIMREG','HIF1A','PLAAT1','IRAG2'])]
signatures_alternative['Severson'] = sev

# Peng 2014
peng = pd.read_excel('../../../Downloads/HRDSignature_Peng2014.xlsx', skiprows=1)
peng = peng['Gene Symbol'].tolist()
peng_replacements = ['MCMBP','FAM170B','DDIAS','SKA3','CEP128','TICRR',
                     'TEDC2','METTL22','ATAD5','KIZ','ISM1','SMIM14',
                     'SNHG32','DSCC1','DEFB1','DDX39A','HJURP','DLGAP5',
                     'DNA2','RETREG1','H1-2','H2BC5','H2AC18','H2BC21',
                     'HSP90AA2P','CREBRF','LOC554223','LOC649679','LOC729843','LOC91431',
                     'VWA5A','ETFRF1','CENPU','MTARC1','BEX3','LRR1',
                     'SRSF2','EPB41L4A-AS1','SLC35G1','TUBB4B','TUBB7P','WRAP53']
peng = [gene if gene in genes_tcga.tolist() else replacement for gene, replacement in zip(peng, peng_replacements)]
peng = [gene for gene in peng if gene in genes_tcga.tolist()]
signatures_alternative['Peng'] = peng

# Save alternative signatures
pd.to_pickle(signatures_alternative, '~/Projects/HRD_TranscriptionalSignature/AlternativeSignatures.pkl')

# Collate centroids for alternative signatures
signature_alternative_centroid_list = {}

# Prepare FPKM-normalized TCGA training data
expr_train_fpkm = expr_train.fpkm_unstrand.T
expr_train_fpkm.index = expr_train_fpkm.index.str[:12]
expr_train_fpkm = np.log2(expr_train_fpkm + 1)
expr_train_fpkm = (expr_train_fpkm - expr_train_fpkm.mean()) / expr_train_fpkm.std()

# Add mutational signature/BRCA defect classification
ann_tcga = pd.read_pickle('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl')
ann_tcga = ann_tcga.set_index('Patient')
ann_tcga = ann_tcga[~ann_tcga['BRCA_status'].isna()]
ann_tcga = ann_tcga[~ann_tcga['BRCA_status'].isin(['PALB2', 'RAD51C'])]

ann_tcga['group'] = ann_tcga['BRCA_status']
ann_tcga.loc[(ann_tcga['HRD'] == 'HRD') & (ann_tcga['BRCA_status'] == 'none'), 'group'] = 'HRD_BRCA+'
ann_tcga.loc[ann_tcga['group'] == 'none', 'group'] = 'HR-proficient'
ann_tcga['group'] = pd.Categorical(ann_tcga['group'], categories=['HR-proficient', 'HRD_BRCA+', 'BRCA1', 'BRCA2'])

samples_intersect = expr_train_fpkm.index.intersection(ann_tcga.index)

ann_tcga = ann_tcga.loc[samples_intersect]
expr_train_fpkm = expr_train_fpkm.loc[samples_intersect]

# For the remaining three signatures:
# Generate average centroids and add to signature list
for signature in ['Severson', 'PARPi7', 'CIN70', 'Peng']:
    sig_genes = signatures_alternative[signature]
    
    centroid_signature = pd.DataFrame({
        'HRD': expr_train_fpkm.loc[ann_tcga['HRD'] == 'HRD', sig_genes].mean(),
        'HR_proficient': expr_train_fpkm.loc[ann_tcga['HRD'] == 'HR-proficient', sig_genes].mean(),
        'BRCA1': expr_train_fpkm.loc[ann_tcga['group'] == 'BRCA1', sig_genes].mean(),
        'BRCA2': expr_train_fpkm.loc[ann_tcga['group'] == 'BRCA2', sig_genes].mean(),
        'HRD_BRCApos': expr_train_fpkm.loc[ann_tcga['group'] == 'HRD_BRCA+', sig_genes].mean(),
        'HR_BRCA_proficient': expr_train_fpkm.loc[ann_tcga['group'] == 'HR-proficient', sig_genes].mean()
    })
    
    if signature == 'Severson':
        centroid_signature['HRD'] = sev_init['BRCA1ness template Pearson correlations']
        centroid_signature['HR_proficient'] = sev_init['non-BRCAness template Pearson correlations']
    
    signature_alternative_centroid_list[signature] = centroid_signature

pd.to_pickle(signature_alternative_centroid_list, '~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.pkl')


# notes:

# Need to make sure that there is Python equivalent for TCGAbiolinks. Pretty sure there is an API.

