import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set working directory
os.chdir('/path/to/Projects/HRD_TranscriptionalSignature/CellphoneDB/Results/')

# Define cell type of interest
cell_type_of_interest = 'Myeloid_cell'

# Process Qian data
qian = pd.read_csv('Qian2020/significant_means.txt', sep='\t')
qian = qian[['interacting_pair', 
             f'{cell_type_of_interest}.CancerHR.proficient',
             f'{cell_type_of_interest}.CancerHRD']]

qian['status'] = np.select([
    qian.iloc[:, 1].notna() & qian.iloc[:, 2].isna(),
    qian.iloc[:, 1].isna() & qian.iloc[:, 2].notna(),
    qian.iloc[:, 1].notna() & qian.iloc[:, 2].notna()
], [
    'ONLY HR-proficient',
    'ONLY HRD',
    'Both'
], default=np.nan)

qian = qian.dropna(subset=['status'])
qian['Dataset'] = 'Qian2020'

# Process Bassez data
bassez = pd.read_csv('Bassez2021/significant_means.txt', sep='\t')
bassez = bassez[['interacting_pair', 
                 f'{cell_type_of_interest}.Cancer_cellHR.proficient',
                 f'{cell_type_of_interest}.Cancer_cellHRD']]

bassez['status'] = np.select([
    bassez.iloc[:, 1].notna() & bassez.iloc[:, 2].isna(),
    bassez.iloc[:, 1].isna() & bassez.iloc[:, 2].notna(),
    bassez.iloc[:, 1].notna() & bassez.iloc[:, 2].notna()
], [
    'ONLY HR-proficient',
    'ONLY HRD',
    'Both'
], default=np.nan)

bassez = bassez.dropna(subset=['status'])
bassez['Dataset'] = 'Bassez2021'

# Merge datasets
bassez.columns = qian.columns
df = pd.concat([qian, bassez])
df = df.rename(columns={df.columns[1]: 'int_HRproficient', df.columns[2]: 'int_HRD'})

# Reshape data for plotting
df_long = df.melt(id_vars=['interacting_pair', 'status', 'Dataset'],
                  value_vars=['int_HRproficient', 'int_HRD'],
                  var_name='CancerStatus', value_name='Interaction')

df_long['group'] = df_long['Dataset'] + '_' + df_long['CancerStatus']
df_long['group'] = pd.Categorical(df_long['group'],
                                  categories=['Qian2020_int_HRD', 'Bassez2021_int_HRD',
                                              'Qian2020_int_HRproficient', 'Bassez2021_int_HRproficient'])

# Plot
plt.figure(figsize=(20, 8))
g = sns.scatterplot(data=df_long, x='interacting_pair', y='group', size='Interaction', 
                    hue='status', size_norm=(20, 200), legend='brief')

# Customize the plot
plt.title(f'Cell type: {cell_type_of_interest}')
plt.axhline(y=1.5, color='red', linestyle='--')
plt.legend(title='', loc='upper left', bbox_to_anchor=(1, 1))
plt.xlabel('Interacting Pair')
plt.ylabel('')

# Rotate x-axis labels
plt.xticks(rotation=90)

# Adjust layout and save
plt.tight_layout()
plt.savefig(f'/path/to/Projects/HRD_TranscriptionalSignature/Figures/Supp_CPDB_Interactions_{cell_type_of_interest}.pdf',
            bbox_inches='tight')
plt.close()