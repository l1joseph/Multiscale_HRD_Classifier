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

def simulate_likelihood_calc(input_data, cluster_distributions, cluster_assign, sample_size):
    sims = pd.DataFrame(0, index=input_data.index, columns=input_data.columns)
    
    for i in range(len(input_data)):
        s = pd.Series(np.random.choice(input_data.columns, size=sample_size,
                                       p=input_data.iloc[i] / input_data.iloc[i].sum())).value_counts()
        sims.loc[input_data.index[i], s.index] = s
    
    posteriors = likelihood_calc(sims, cluster_distributions, cluster_assign)
    posteriors['PhenoTrue'] = cluster_assign
    
    roc_apobec = roc_curve(posteriors['PhenoTrue'] == 'APOBEC', posteriors['APOBEC'])
    roc_sbs3 = roc_curve(posteriors['PhenoTrue'] == 'SBS3', posteriors['SBS3'])
    roc_sbs5 = roc_curve(posteriors['PhenoTrue'] == 'SBS5', posteriors['SBS5'])
    
    auc_values = [auc(roc_apobec[0], roc_apobec[1]),
                  auc(roc_sbs3[0], roc_sbs3[1]),
                  auc(roc_sbs5[0], roc_sbs5[1])]
    
    return {'auc_values': auc_values, 'posteriors': posteriors}

cluster_distribution_list = []
indel_alterations = [1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5]
for indel_alter in indel_alterations:
    mut_prob_sbs = mut_prob.filter(regex='>')
    mut_prob_id = mut_prob.filter(regex=':')
    
    mut_prob_id *= indel_alter
    for j in range(3):
        mut_prob_sbs.iloc[j] = mut_prob_sbs.iloc[j] * (1-mut_prob_id.iloc[j].sum()) / mut_prob_sbs.iloc[j].sum()
    
    mut_prob_altered = pd.concat([mut_prob_sbs, mut_prob_id], axis=1)
    cluster_distribution_list.append(mut_prob_altered)

def run_likelihood_sims(input_data, sample_size, cluster_distributions_index, cluster_assign, n_simulations):
    cluster_distributions = cluster_distribution_list[cluster_distributions_index]
    
    results_mat = []
    posterior_df = pd.DataFrame()
    
    for i in range(n_simulations):
        np.random.seed(123*i)
        
        print(f'Sample size = {sample_size}, mut.prob_index = {indel_alterations[cluster_distributions_index]}, Running simulation {i+1} of {n_simulations}...{pd.Timestamp.now()}')
        
        run_i = simulate_likelihood_calc(input_data, cluster_distributions, cluster_assign, sample_size)
        results_mat.append(run_i['auc_values'])
        
        run_i['posteriors']['Run'] = i
        posterior_df = pd.concat([posterior_df, run_i['posteriors']])
    
    posterior_df['sample_size'] = sample_size
    posterior_df['indel_alter'] = indel_alterations[cluster_distributions_index]
    
    results = pd.DataFrame({
        'Pheno': np.repeat(['APOBEC', 'SBS3', 'SBS5'], len(results_mat)),
        'AUC': np.concatenate(results_mat),
        'sample_size': sample_size,
        'indel_alteration': indel_alterations[cluster_distributions_index]
    })
    
    return {'results': results, 'posteriors': posterior_df}

nSim = 100
sampleSize = 25
res_full_25 = pd.DataFrame()

np.random.seed(123)
for i in range(len(cluster_distribution_list)):
    print(f'Indel alteration = {indel_alterations[i]}')
    
    res_i_25 = run_likelihood_sims(mut_complete, sampleSize, i, sigs_clust3['Phenotype'], nSim)
    res_full_25 = pd.concat([res_full_25, res_i_25['results']])

res_25_sbs3 = res_full_25[res_full_25['Pheno'] == 'SBS3']

my_comparisons = [('0.2', '1'), ('1', '5')]

plt.figure(figsize=(12, 4))
sns.boxplot(data=res_25_sbs3, x='indel_alteration', y='AUC')
plt.ylim(0.65, 1.05)
plt.title('sample size = 25')
plt.xlabel('Indel alteration')
plt.xticks(range(9), ['1/5', '1/4', '1/3', '1/2', '1', '2', '3', '4', '5'])

for i, comp in enumerate(my_comparisons):
    x1, x2 = [list(res_25_sbs3['indel_alteration'].unique()).index(float(comp[0])),
               list(res_25_sbs3['indel_alteration'].unique()).index(float(comp[1]))]
    y, h, col = res_25_sbs3['AUC'].max() + 0.02, 0.02, 'k'
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, f"p={stats.ttest_ind(res_25_sbs3[res_25_sbs3['indel_alteration'] == float(comp[0])]['AUC'],
                                                   res_25_sbs3[res_25_sbs3['indel_alteration'] == float(comp[1])]['AUC']).pvalue:.2e}",
             ha='center', va='bottom', color=col)

plt.savefig('Figures/Supp_SimulationsIndelLikelihoods_size25.pdf')
plt.close()

# Similar code for sample sizes 50 and 100

# notes:

# In the R file, mixture modeling was done using mclust. However for simplicity's sake, we used KMeans in here.
# To use mixture modeling, we can use sklearn.mixture.GaussianMixture.
# We also might want to create a t-test function if we plan on repeating it multiple times
# Need to create plots for 50 and 100 sample size scenarios.
# Also need to verify random seed provides reproducible results.