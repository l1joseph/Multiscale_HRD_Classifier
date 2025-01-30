import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LinearRegression, LogisticRegression, ElasticNet, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import roc_curve, auc, confusion_matrix, classification_report, mean_squared_error, r2_score
from sklearn.model_selection import ParameterGrid, train_test_split
from scipy.stats import pearsonr
from sklearn.model_selection import ParameterSampler
import joblib

# metadata
ann_tcga = pd.read_csv('../data/toga.breast.brca.status.txt', sep='\t', index_col=0)
hrd_scores = pd.read_excel('../data/tcga.hrdscore.xlsx', index_col=0)

# rnaseq versions
fpkm = pd.read_csv('../data/tcga.brca.rnaseq.unstranded.fpkm.counts.matrix.txt', sep='\t', index_col=0)
deconvo = pd.read_csv('../data/Deconvo2.csv',  index_col=0)
tpm = pd.read_csv('../data/tpm.csv', index_col=0, low_memory=False)

# ann_tcga[ann_tcga['event.PALB2'].ne('0')]
# ann_tcga[ann_tcga['event.RAD51C'].ne('0')]
metadata = ann_tcga[~ann_tcga['event.RAD51C'].ne('0')]
metadata = metadata[~metadata['event.PALB2'].ne('0')]
metadata = metadata[metadata['event.BRCA1'] != '1']
metadata.index = metadata.index.str.replace('.', '-', regex=False)
hrd_scores.index = hrd_scores.index.map(lambda x: x[:12])
metadata = metadata.merge(hrd_scores[['HRD-sum']], left_index=True, right_index=True, how='inner')


fpkm = fpkm.set_index('Case ID')
fpkm = fpkm[fpkm['Sample Type'] == 'Primary Tumor']
fpkm = fpkm.drop(columns=["Sample ID","Sample Type"])


dick = {}
for i, v in enumerate(fpkm.columns):
    # print(i,v)
    dick[v.split('|')[0]] =v.split('|')[1]
tpm = tpm.rename(columns=dick)
tpm = tpm[~tpm.index.duplicated(keep='first')]

# Keeps only Primary Tumor
tpm = tpm.loc[tpm.index.str[13:15] == "01"]

# deconvo = deconvo.sort_index()
# deconvo = np.log2(deconvo + 1)
fpkm.columns = [col.split('|')[1] if '|' in col else col for col in fpkm.columns]
fpkm.index = fpkm.index.map(lambda x: x[:12])
fpkm = fpkm.loc[fpkm.index.intersection(metadata.index)]
fpkm.sort_index(inplace=True)
fpkm = fpkm.rename_axis("fpkm", axis="index")
fpkm = fpkm.apply(pd.to_numeric, errors='coerce')
fpkm.fillna(0, inplace=True)

deconvo.index = deconvo.index.map(lambda x: x[:12])
deconvo = deconvo.loc[deconvo.index.intersection(metadata.index)]
deconvo.sort_index(inplace=True)
deconvo = deconvo.rename_axis("deconvo", axis="index")
deconvo.fillna(0, inplace=True)

tpm.index = tpm.index.map(lambda x: x[:12])
tpm = tpm.loc[tpm.index.intersection(metadata.index)]
tpm.sort_index(inplace=True)
tpm = tpm.rename_axis("tpm", axis="index")
tpm = tpm.apply(pd.to_numeric, errors='coerce')
tpm.fillna(0, inplace=True)



print(f"fpkm shape{fpkm.shape}")
print(f"deconvo shape{deconvo.shape}")
print(f"tpm shape{tpm.shape}")


# plt.figure(figsize=(10, 6))
# sns.histplot(data=metadata, x='HRD-sum', hue='event.PAM50', multiple='layer', kde=True)
# plt.title('Histogram of HRD-sum with PAM50 Subtypes')
# plt.xlabel('HRD-sum')
# plt.ylabel('Frequency')
# plt.show()
pam50_counts = metadata['event.PAM50'].value_counts()
print(pam50_counts)


def downsampling_lumA(metadata, lumA_cutoff):
    lumA_HRD = metadata[(metadata['event.PAM50'] == 'LumA') & (metadata['HRD-sum'] >= lumA_cutoff)]
    lumA_HRP = metadata[(metadata['event.PAM50'] == 'LumA') & (metadata['HRD-sum'] < lumA_cutoff)]
    # print(lumA_HRP.shape, lumA_HRD.shape)
    if lumA_HRP.shape[0] < lumA_HRD.shape[0]:
        print(f"Not enough HRP samples ({lumA_HRP.shape[0]}) to match HRD count ({lumA_HRD.shape[0]}). Skipping downsampling.")
        return metadata, pd.DataFrame()
    lumA_HRP_downsampled = lumA_HRP.sample(n=lumA_HRD.shape[0], random_state=42)
    df_downsampled = pd.concat([lumA_HRD, lumA_HRP_downsampled])
    df_downsampled = pd.concat([df_downsampled, metadata[metadata['event.PAM50'] != 'LumA']])
    # print(df_downsampled['event.PAM50'].value_counts())
    print(f"number of samples left: {df_downsampled.shape}")
    unused_majority = lumA_HRP.loc[~lumA_HRP.index.isin(lumA_HRP_downsampled.index)]
    return df_downsampled, unused_majority

# df_downsampled, removed_samples = downsampling_lumA(metadata, 23)

def add_back_test(rna_df, removed_samples, X_test,y_test):
    add_back_features = rna_df.loc[rna_df.index.intersection(removed_samples.index)]
    add_back_features = add_back_features.sort_index()

    add_back_y = removed_samples.loc[removed_samples.index.intersection(rna_df.index)]
    add_back_y = add_back_y.sort_index()

    X_test = pd.concat([X_test, add_back_features])
    y_test = pd.concat([y_test, add_back_y['HRD-sum'].squeeze()])
    return X_test, y_test

# X_test,y_test = add_back_test(deconvo, removed_samples, X_test, y_test)

def sigmoid_transform(values, shift=0, scale=1):
    return 1 / (1 + np.exp(-scale * (values - shift)))
def binary_hrd(values, threshold):
    return (values >= threshold).astype(int)

def plot_test_train_pam50_dist(metadata, X_train, X_test):
    plt.figure(figsize=(10, 6))
    sns.histplot(data=metadata.loc[metadata.index.intersection(X_test.index)], x='HRD-sum', hue='event.PAM50', multiple='layer', kde=True)
    plt.title('Test Distribution')
    plt.xlabel('HRD-sum')
    plt.ylabel('Frequency')
    plt.show()

    plt.figure(figsize=(10, 6))
    sns.histplot(data=metadata.loc[metadata.index.intersection(X_train.index)], x='HRD-sum', hue='event.PAM50', multiple='layer', kde=True)
    plt.title('Train Distribution')
    plt.xlabel('HRD-sum')
    plt.ylabel('Frequency')
    plt.show()

    # pam50_counts = df_downsampled['event.PAM50'].value_counts()
    # print(pam50_counts)
    # print(df_downsampled.shape)
# plot_test_train_pam50_dist(metadata, X_train, X_test)



alphas = [0.01, 0.1, 0.25, 0.5, 1.0]
l1_ratios = [0.1, 0.5, 0.7, 0.9]
rna_seqs = [deconvo, tpm, fpkm]
downsample = [(True, False),(True,False), (False,False)]
downsample_thresholds = [x for x in range(10, 80, 5)]
softlabels = ["None", "Sigmoid", "Binary"]
softlabel_thresholds = [x for x in range(10, 80, 5)]
softlabel_gradients = np.arange(0, 1, 0.1)
normalization = ['StandardScaler','log2','None']



best_model = None
best_metrics = {'Mean Squared Error': float('inf'), 'R^2 Score': -float('inf')}
best_params = {}

# Grid search
param_distributions = {'alpha': alphas, 'l1_ratio': l1_ratios,
                            'rna-seq': rna_seqs, 'downsample': downsample,
                            'downsample_thresholds': downsample_thresholds,'softlabels': softlabels,
                            'softlabel_thresholds': softlabel_thresholds, 'softlabel_gradients': softlabel_gradients,
                            'normalization': normalization}


n_iter = 3000  # You can adjust this value based on your computational resources

# Generate random combinations of parameters
random_params = list(ParameterSampler(param_distributions, n_iter=n_iter, random_state=42))

def get_variable_name(var):
    conv = {deconvo:'deconvo', tpm:'tpm', fpkm:'fpkm'}
    return conv[var]

for iteration, params in enumerate(random_params):
    # print(params)

    print({k: (v.index.name if k == 'rna-seq' else v) for k, v in params.items()})
    
    # Split data into training and test sets and get intersecting indices
    features_df = params['rna-seq']

    metadata_truncated = metadata.loc[metadata.index.intersection(features_df.index)]
    # features_df = features_df.loc[features_df.index.intersection(metadata_truncated.index)]

    labels_df = metadata_truncated['HRD-sum']
    labels_df = labels_df.sort_index()
    features_df = features_df.sort_index()
    
    debug_features = features_df.copy()
    debug_labels = labels_df.copy()

    if params['normalization'] == 'StandardScaler':
        scaler = StandardScaler()
        features_df = pd.DataFrame(scaler.fit_transform(features_df), index=features_df.index, columns=features_df.columns)
    elif params['normalization'] == 'log2':
        features_df = np.log2(features_df + 1)

    # Apply soft labels
    if params['softlabels'] == "Sigmoid":
        labels_df = sigmoid_transform(labels_df, params['softlabel_thresholds'], params['softlabel_gradients'])
    elif params['softlabels'] == "Binary":
        labels_df = binary_hrd(labels_df, params['softlabel_thresholds'])


    # Downsample LumA samples
    if params['downsample'][0]:
        df_downsampled, removed_samples = downsampling_lumA(metadata_truncated, params['downsample_thresholds'])
        features_df = features_df.loc[features_df.index.intersection(df_downsampled.index)]
        labels_df = df_downsampled.loc[df_downsampled.index.intersection(features_df.index), 'HRD-sum']
    labels = labels_df.squeeze()

    #TEST TRAIN SPLIT
    X_train, X_test, y_train, y_test = train_test_split(features_df, labels, test_size=0.2, random_state=42)
    
    if params['downsample'][1]:
        X_test,y_test = add_back_test(features_df, removed_samples, X_test, y_test)
        
    
    model = ElasticNet(alpha=params['alpha'], l1_ratio=params['l1_ratio'], max_iter=1000, random_state=42)
    model.fit(X_train, y_train)
    
    # Predict on test set
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    
    # Update best model if current is better
    # if mse < best_metrics['Mean Squared Error']:
    if mse < best_metrics['Mean Squared Error'] or (mse == best_metrics['Mean Squared Error'] and r2 > best_metrics['R^2 Score']):
        best_model = model
        best_metrics = {'Mean Squared Error': mse, 'R^2 Score': r2}
        best_params = params
        best_X_train = X_train
        best_X_test = X_test
    if iteration % 3 == 0:
        joblib.dump(best_model, 'model.joblib')

# Scatter plot of predictions vs actual values for the best model
y_pred_best = best_model.predict(X_test)
plt.figure(figsize=(8, 6))
plt.scatter(y_test, y_pred_best, alpha=0.6, label='Predicted vs Actual')
plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'k--', label='Ideal Prediction')
plt.xlabel('Actual Values')
plt.ylabel('Predicted Values')
plt.title('Elastic Net Regression: Predicted vs Actual (Best Model)')
plt.legend(loc="upper left")
plt.show()
plot_test_train_pam50_dist(metadata, best_X_train, best_X_test)
print(f"Best Parameters: {best_params}")
print(f"Mean Squared Error: {best_metrics['Mean Squared Error']:.3f}")
print(f"R^2 Score: {best_metrics['R^2 Score']:.3f}")

# return best_model, best_metrics, best_params
elastic_model = best_model
metrics = best_metrics
params = best_params
joblib.dump(best_model, 'model.joblib')
