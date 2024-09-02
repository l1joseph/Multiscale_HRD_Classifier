plt.figure(figsize=(6, 3))
sns.barplot(data=res_BRCA_plot, x='HRD', y='n', hue='BRCA_subtype')
plt.xlabel('')
plt.ylabel('% Samples')
plt.legend(title='')
plt.gca().invert_yaxis()
plt.savefig('Figures/Figure2/BRCAsubtypevsHRD.pdf')
plt.close()

## Plot supplementary figures

# 1. HRD hallmarks: HR-proficient vs HRD_BRCA+ vs HRD_BRCA-
ann_tcga2 = ann_tcga.copy()
ann_tcga2 = ann_tcga2[ann_tcga2['BRCA_status'].notna()]
ann_tcga2['group'] = ann_tcga2['HRD']
ann_tcga2.loc[(ann_tcga2['HRD'] == 'HRD') & (ann_tcga2['BRCA_status'] == 'none'), 'group'] = 'HRD_BRCA+'
ann_tcga2.loc[(ann_tcga2['HRD'] == 'HRD') & (ann_tcga2['BRCA_status'] != 'none'), 'group'] = 'HRD_BRCA-'
ann_tcga2['group'] = pd.Categorical(ann_tcga2['group'], 
                                    categories=['HR-proficient', 'HRD_BRCA+', 'HRD_BRCA-'], 
                                    ordered=True)

df_supp1 = pd.merge(ann_tcga2[['Patient', 'group']], 
                    tcga_HRDHallmarks[['Patient', 'HRD_index', 'CX3', 'log_POLQ_FPKM', 'ProliferativeCapacity']])
df_supp1_melted = df_supp1.melt(id_vars=['Patient', 'group'], 
                                var_name='Hallmark', 
                                value_name='score')

plt.figure(figsize=(10, 5))
sns.boxplot(data=df_supp1_melted, x='group', y='score', hue='group')
plt.xlabel('')
plt.ylabel('')
plt.legend(title='')
plt.xticks([])
plt.facet_grid(cols='Hallmark', scales='free')
statannot.add_stat_annotation(data=df_supp1_melted, x='group', y='score', hue='Hallmark',
                              box_pairs=[('HR-proficient', 'HRD_BRCA+'),
                                         ('HRD_BRCA+', 'HRD_BRCA-'),
                                         ('HR-proficient', 'HRD_BRCA-')],
                              test='t-test_ind', text_format='p={:.3e}',
                              loc='outside', verbose=2)
plt.savefig('Figures/Supp_TCGAhallmarksBRCAPositive.pdf')
plt.close()

# 2. HRD hallmarks: HR-proficient vs p(HRD) > 0.5 vs p(HRD) > 0.79
ann_tcga3 = ann_tcga[['Patient', 'HRD', 'HRD_prob']].copy()
ann_tcga3['group'] = pd.cut(ann_tcga3['HRD_prob'], 
                            bins=[-np.inf, 0.5, 0.79, np.inf],
                            labels=['HR-proficient', 'p(HRD) > 0.5', 'p(HRD) > 0.79'])

df_supp2 = pd.merge(ann_tcga3[['Patient', 'group']], 
                    tcga_HRDHallmarks[['Patient', 'HRD_index', 'CX3', 'log_POLQ_FPKM', 'ProliferativeCapacity']])
df_supp2_melted = df_supp2.melt(id_vars=['Patient', 'group'], 
                                var_name='Hallmark', 
                                value_name='score')

plt.figure(figsize=(10, 5))
sns.boxplot(data=df_supp2_melted, x='group', y='score', hue='group')
plt.xlabel('')
plt.ylabel('')
plt.legend(title='')
plt.xticks([])
plt.facet_grid(cols='Hallmark', scales='free')
statannot.add_stat_annotation(data=df_supp2_melted, x='group', y='score', hue='Hallmark',
                              box_pairs=[('HR-proficient', 'p(HRD) > 0.5'),
                                         ('p(HRD) > 0.5', 'p(HRD) > 0.79'),
                                         ('HR-proficient', 'p(HRD) > 0.79')],
                              test='t-test_ind', text_format='p={:.3e}',
                              loc='outside', verbose=2)
plt.savefig('Figures/Supp_TCGAhallmarksHRDthresholds.pdf')
plt.close()

print("All plots have been generated and saved.")
