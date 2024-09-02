import os
import pandas as pd
import numpy as np
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
import seaborn as sns

# Set working directory
os.chdir('~/Data/SMC_BRCA/')

# Load required R packages
pandas2ri.activate()
maftools = importr('maftools')
sigminer = importr('sigminer')

# Load SMC data and tally mutation contributions
smc_mutations = pd.read_csv('data_mutations.txt', sep='\t')

mut_maf = maftools.read_maf(maf=smc_mutations, 
                            vc_nonSyn=list(smc_mutations.Variant_Classification.unique()))

mt_tally_smc = sigminer.sig_tally(
    mut_maf,
    ref_genome='BSgenome.Hsapiens.UCSC.hg19',
    mode='ALL',
    use_syn=True
)

smc_muts = pd.concat([mt_tally_smc.rx2('SBS_96'), mt_tally_smc.rx2('ID_83')], axis=1)

# Load relevant data for classifier
os.chdir('~/Projects/HRD_MutationalSignature/')

# Prior cluster mean distributions
r('load("Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.Rdata")')
mut_dists_mean = r('mut.dists_mean')

# Signature Phenotype Assignment
r('load("Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.Rdata")')
pheno_assigned = r('ann$Phenotype')

# Likelihood function
def likelihood_calc(input_data, cluster_distributions, cluster_assign):
    log_likelihoods = np.zeros((input_data.shape[0], cluster_distributions.shape[0]))
    for i in range(input_data.shape[0]):
        log_likelihoods[i, :] = np.sum(np.log10(cluster_distributions) * input_data.iloc[i, :], axis=1)
    
    marginal_probs = pd.Series(cluster_assign).value_counts() / len(cluster_assign)
    marginal_probs = marginal_probs[log_likelihoods.columns]
    
    log_posteriors = log_likelihoods + np.log10(marginal_probs.values)
    
    final_probs = 10**(log_posteriors - np.max(log_posteriors, axis=1, keepdims=True))
    final_probs = final_probs / np.sum(final_probs, axis=1, keepdims=True)
    
    return pd.DataFrame(final_probs, index=input_data.index, columns=cluster_distributions.index)

# Apply log-likelihood approach
results_smc_loglik = likelihood_calc(input_data=smc_muts, 
                                     cluster_distributions=mut_dists_mean,
                                     cluster_assign=pheno_assigned)

results_smc_df = pd.DataFrame({
    'Patient': results_smc_loglik.index,
    'Phenotype_Assigned': results_smc_loglik.idxmax(axis=1),
    'Phenotype_Assigned.prob': results_smc_loglik.max(axis=1),
    'HRD_prob': results_smc_loglik.filter(regex='HRD').sum(axis=1)
})
results_smc_df['HRD'] = np.where(results_smc_df.HRD_prob >= 0.79, 'HRD', 'HR-proficient')

results_smc_df.to_pickle('Results/SMC_HRD_resultsSummary.pkl')

# Match with BRCA status
data_brca = pd.read_excel('~/Data/SMC_BRCA/SMC_BRCA.BrcaStatus.xlsx', skiprows=2)

results_smc_df['sample_id'] = results_smc_df.Patient.str[14:]
results_smc_df = pd.merge(results_smc_df, data_brca[['sample_id', 'gene_symbol']], how='left')

# Match with additional clinical data
smc_clinic = pd.read_csv('~/Data/SMC_BRCA/data_clinical_sample.txt', sep='\t', skiprows=4)

results_smc_df = pd.merge(results_smc_df, smc_clinic[['PATIENT_ID', 'SUBTYPE_CONSENSUS']], 
                          left_on='Patient', right_on='PATIENT_ID')
results_smc_df = results_smc_df.rename(columns={'SUBTYPE_CONSENSUS': 'Subtype', 'gene_symbol': 'BRCA_defect'})

# Plot results
ann_smc = results_smc_df[['BRCA_defect', 'HRD_prob', 'HRD', 'Phenotype_Assigned', 'Subtype']]
ann_smc = ann_smc.set_index(results_smc_df.Patient)
ann_smc = ann_smc.sort_values('Phenotype_Assigned')

results_smc_plot = results_smc_loglik.loc[ann_smc.index, sorted(results_smc_loglik.columns)]

# Plotting
plt.figure(figsize=(12, 8))
sns.heatmap(results_smc_plot.T, cmap='Blues', cbar_kws={'label': 'Probability'})
plt.title('SMC HRD Classification Heatmap')
plt.tight_layout()
plt.savefig('Figures/Supp_SMCheatmap.pdf')
plt.close()

# Compare with BRCA subtype
print(pd.crosstab(ann_smc.HRD_prob >= 0.79, ann_smc.BRCA_defect, margins=True))

# Plotting barplots
ann_smc['BRCA_status'] = np.where(ann_smc.BRCA_defect.notna(), 'BRCA-defective', 'BRCA+')
ann_smc['BRCA_status'] = pd.Categorical(ann_smc.BRCA_status, categories=['BRCA-defective', 'BRCA+'])

ann_smc['HRDgroup'] = pd.cut(ann_smc.HRD_prob, 
                             bins=[-np.inf, 0.5, 0.79, np.inf], 
                             labels=['HR-proficient', 'HRD > 0.5', 'HRD > 0.79'])

# BRCA status plot
plt.figure(figsize=(6, 4))
sns.countplot(x='BRCA_status', hue='HRDgroup', data=ann_smc)
plt.title('HRD Classification by BRCA Status')
plt.tight_layout()
plt.savefig('Figures/Supp_SMCbrcaClassification.pdf')
plt.close()

# Subtype plot
plt.figure(figsize=(8, 4))
sns.countplot(x='Subtype', hue='HRDgroup', data=ann_smc)
plt.title('HRD Classification by Subtype')
plt.tight_layout()
plt.savefig('Figures/Supp_SMCsubtypeClassification.pdf')
plt.close()
