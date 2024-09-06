import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import roc_curve, auc


# Load ICGC deconstructSigs data
sigs_complete = pd.read_pickle('Results/ICGC_BRCA_deconstructSigs_Cutoff0.01_SBSandIDnormalised.pkl')
sigs_complete = sigs_complete.filter(regex='SBS')

# Run mixture modelling using KMeans
# Note: This is a simplification. We should consider using sklearn.mixture.GaussianMixture
kmeans = KMeans(n_clusters=3, random_state=42).fit(sigs_complete)
sigs_clust3 = pd.DataFrame({'Cluster_num': kmeans.labels_})
sigs_clust3['Phenotype'] = sigs_clust3['Cluster_num'].map({0: 'SBS5', 1: 'SBS3', 2: 'APOBEC'})

# Load in ICGC mutation tallies
mt_tally_brca_wgs = pd.read_pickle('Data/BRCA_UKEU_mt_tally.pkl')
mut_complete = pd.concat([mt_tally_brca_wgs['SBS_96'], mt_tally_brca_wgs['ID_83']], axis=1)

# Separate out mut_complete based on cluster assignment
mut_apobec = mut_complete[sigs_clust3['Phenotype'] == 'APOBEC']
mut_sbs3 = mut_complete[sigs_clust3['Phenotype'] == 'SBS3']
mut_sbs5 = mut_complete[sigs_clust3['Phenotype'] == 'SBS5']

def collate_function(input_data, variant_type='ALL'):
    if variant_type == 'SBS':
        input_final = input_data.filter(regex='>')
    elif variant_type == 'ID':
        input_final = input_data.filter(regex=':')
    else:
        input_final = input_data
    
    total = input_final.sum()
    dist = total / total.sum()
    return dist

mut_prob_apobec = collate_function(mut_apobec)
mut_prob_sbs3 = collate_function(mut_sbs3)
mut_prob_sbs5 = collate_function(mut_sbs5)

mut_prob = pd.DataFrame([mut_prob_apobec, mut_prob_sbs3, mut_prob_sbs5])
mut_prob.index = ['APOBEC', 'SBS3', 'SBS5']

def likelihood_calc(input_data, cluster_distributions, cluster_assign):
    likelihoods = pd.DataFrame(index=input_data.index, columns=cluster_distributions.index)
    for i, row in input_data.iterrows():
        likelihoods.loc[i] = cluster_distributions.apply(lambda x: np.prod(x ** row))
    
    marginal_probs = cluster_assign.value_counts() / len(cluster_assign)
    
    posteriors = likelihoods.mul(marginal_probs, axis=1)
    posteriors = posteriors.div(posteriors.sum(axis=1), axis=0)
    return posteriors

def simulate_likelihood_calc(input_data, cluster_distributions, cluster_assign, sample_size, indel_prop=0):
    data_sbs = input_data.filter(regex='>')
    data_id = input_data.filter(regex=':')
    
    sims = pd.DataFrame(0, index=input_data.index, columns=input_data.columns)
    
    for i in range(len(input_data)):
        s_sbs = pd.Series(np.random.choice(data_sbs.columns, size=int(sample_size * (1 - indel_prop)),
                                           p=data_sbs.iloc[i] / data_sbs.iloc[i].sum())).value_counts()
        sims.loc[input_data.index[i], s_sbs.index] = s_sbs
        
        s_id = pd.Series(np.random.choice(data_id.columns, size=int(sample_size * indel_prop),
                                          p=data_id.iloc[i] / data_id.iloc[i].sum())).value_counts()
        sims.loc[input_data.index[i], s_id.index] += s_id
    
    posteriors = likelihood_calc(sims, cluster_distributions, cluster_assign)
    posteriors['PhenoTrue'] = cluster_assign
    
    roc_apobec = roc_curve(posteriors['PhenoTrue'] == 'APOBEC', posteriors['APOBEC'])
    roc_sbs3 = roc_curve(posteriors['PhenoTrue'] == 'SBS3', posteriors['SBS3'])
    roc_sbs5 = roc_curve(posteriors['PhenoTrue'] == 'SBS5', posteriors['SBS5'])
    
    auc_values = [auc(roc_apobec[0], roc_apobec[1]),
                  auc(roc_sbs3[0], roc_sbs3[1]),
                  auc(roc_sbs5[0], roc_sbs5[1])]
    
    return {'auc_values': auc_values, 'posteriors': posteriors}

def run_likelihood_sims(input_data, sample_size, indel_prop, cluster_distributions, cluster_assign, n_simulations):
    results_mat = []
    posterior_df = pd.DataFrame()
    
    for i in range(n_simulations):
        np.random.seed(i)
        
        print(f'Sample size = {sample_size}, indel_prop = {indel_prop}, Running simulation {i+1} of {n_simulations}...')
        run_i = simulate_likelihood_calc(input_data, cluster_distributions, cluster_assign, sample_size, indel_prop)
        results_mat.append(run_i['auc_values'])
        
        run_i['posteriors']['Run'] = i
        posterior_df = pd.concat([posterior_df, run_i['posteriors']])
    
    posterior_df['sample_size'] = sample_size
    posterior_df['indel_prop'] = indel_prop
    
    results = pd.DataFrame({
        'Pheno': np.repeat(['APOBEC', 'SBS3', 'SBS5'], len(results_mat)),
        'AUC': np.concatenate(results_mat),
        'sample_size': sample_size,
        'indel_prop': indel_prop
    })
    
    return {'results': results, 'posteriors': posterior_df}

nSim = 100
res_full_25 = pd.DataFrame()
res_full_50 = pd.DataFrame()
res_full_100 = pd.DataFrame()
full_posteriors = pd.DataFrame()

for indel_props in np.arange(0, 0.55, 0.05):
    res_i_25 = run_likelihood_sims(mut_complete, 25, indel_props, mut_prob, sigs_clust3['Phenotype'], nSim)
    res_full_25 = pd.concat([res_full_25, res_i_25['results']])
    full_posteriors = pd.concat([full_posteriors, res_i_25['posteriors']])
    
    print(f'Indel proportion = {indel_props}, sample size = 50')
    
    res_i_50 = run_likelihood_sims(mut_complete, 50, indel_props, mut_prob, sigs_clust3['Phenotype'], nSim)
    res_full_50 = pd.concat([res_full_50, res_i_50['results']])
    full_posteriors = pd.concat([full_posteriors, res_i_50['posteriors']])
    
    print(f'Indel proportion = {indel_props}, sample size = 100')
    
    res_i_100 = run_likelihood_sims(mut_complete, 100, indel_props, mut_prob, sigs_clust3['Phenotype'], nSim)
    res_full_100 = pd.concat([res_full_100, res_i_100['results']])
    full_posteriors = pd.concat([full_posteriors, res_i_100['posteriors']])

myComp = [('0', '0.05'), ('0', '0.1'), ('0', '0.15'), ('0', '0.2'), ('0', '0.25')]

res_full_25['best_AUC'] = res_full_25['indel_prop'] == 0.2
plt.figure(figsize=(12, 4))
sns.boxplot(data=res_full_25[res_full_25['Pheno'] == 'SBS3'], x='indel_prop', y='AUC', hue='best_AUC')
plt.axvline(x=1+(1/0.05)*0.0675, linestyle='dashed', color='blue')
plt.legend([],[], frameon=False)
for i, comp in enumerate(myComp):
    x1, x2 = [list(res_full_25['indel_prop'].unique()).index(float(comp[0])),
               list(res_full_25['indel_prop'].unique()).index(float(comp[1]))]
    y, h, col = res_full_25['AUC'].max() + 0.02, 0.02, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, f"p={stats.ttest_ind(res_full_25[(res_full_25['Pheno'] == 'SBS3') & (res_full_25['indel_prop'] == float(comp[0]))]['AUC'],
                                                   res_full_25[(res_full_25['Pheno'] == 'SBS3') & (res_full_25['indel_prop'] == float(comp[1]))]['AUC']).pvalue:.2e}",
             ha='center', va='bottom', color=col)
plt.title('sample size = 25')
plt.xlabel('Indel Proportions')
plt.savefig('Figures/Supp_ICGCsimulations_indelProps_sim25.pdf')
plt.close()

# Similar plots for sample sizes 50 and 100


# notes:

# In the R file, mixture modeling was done using mclust. However for simplicity's sake, we used KMeans in here.
# To use mixture modeling, we can use sklearn.mixture.GaussianMixture.
# We also might want to create a t-test function if we plan on repeating it multiple times
# Need to create plots for 50 and 100 sample size scenarios.
# Also need to verify random seed provides reproducible results.