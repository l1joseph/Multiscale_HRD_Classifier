import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc
from scipy import stats

def f_score(confusion_matrix, weight=1):
    """
    Calculate F-score from confusion matrix
    """
    tn, fp, fn, tp = confusion_matrix.ravel()
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)
    f_score = (1 + weight**2) * (precision * recall) / ((weight**2 * precision) + recall)
    return pd.Series([recall, precision, f_score], index=['recall', 'precision', 'f_score'])

# Load TCGA assignments with BRCA-defect and ER status labels
ann_tcga = pd.read_pickle('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl')

# Process the data
ann_tcga = ann_tcga[ann_tcga['ER_status'].isin(['Negative', 'Positive'])]
ann_tcga['HR_defect'] = np.where(ann_tcga['BRCA_status'] != 'none', 'defective', 'proficient')
ann_tcga['HR_defect'] = pd.Categorical(ann_tcga['HR_defect'], categories=['defective', 'proficient'])

# Initialize dataframe to store results
hrd_thresholds = pd.DataFrame(columns=['p_HRD', 'recall', 'precision', 'F-score'])

# Calculate F-score for different thresholds
for i in np.arange(0, 1.01, 0.01):
    ann_tcga['HRD'] = np.where(ann_tcga['HRD_prob'] > i, 'HRD', 'HR-proficient')
    ann_tcga['HRD'] = pd.Categorical(ann_tcga['HRD'], categories=['HRD', 'HR-proficient'])
    
    confusion_matrix = pd.crosstab(ann_tcga['HRD'], ann_tcga['HR_defect'])
    f_scores = f_score(confusion_matrix)
    
    hrd_thresholds = hrd_thresholds.append(
        pd.Series([i] + list(f_scores), index=['p_HRD', 'recall', 'precision', 'F-score']),
        ignore_index=True
    )

# Plot results
plt.figure(figsize=(10, 6))
sns.lineplot(x='p_HRD', y='value', hue='variable', 
             data=pd.melt(hrd_thresholds, ['p_HRD']))
plt.axvline(x=hrd_thresholds.loc[hrd_thresholds['F-score'].idxmax(), 'p_HRD'], 
            color='red', linestyle='--')
plt.xlabel('p(HRD)')
plt.ylabel('Score')
plt.title('Performance Metrics vs HRD Probability Threshold')
plt.savefig('~/Projects/HRD_MutationalSignature/Figures/Supp_TCGAoptimalHRDthreshold.pdf')
plt.close()

# Plot AUC curve
fpr, tpr, thresholds = roc_curve(ann_tcga['HR_defect'] == 'defective', ann_tcga['HRD_prob'])
roc_auc = auc(fpr, tpr)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.savefig('~/Projects/HRD_MutationalSignature/Figures/Supp_TCGAHRDthresholdAUC.pdf')
plt.close()

# Repeat analysis for ER-positive and ER-negative subsets
for er_status in ['Positive', 'Negative']:
    ann_tcga_er = ann_tcga[ann_tcga['ER_status'] == er_status]
    
    hrd_thresholds_er = pd.DataFrame(columns=['p_HRD', 'recall', 'precision', 'F-score'])
    
    for i in np.arange(0, 1.01, 0.01):
        ann_tcga_er['HRD'] = np.where(ann_tcga_er['HRD_prob'] > i, 'HRD', 'HR-proficient')
        ann_tcga_er['HRD'] = pd.Categorical(ann_tcga_er['HRD'], categories=['HRD', 'HR-proficient'])
        
        confusion_matrix = pd.crosstab(ann_tcga_er['HRD'], ann_tcga_er['HR_defect'])
        f_scores = f_score(confusion_matrix)
        
        hrd_thresholds_er = hrd_thresholds_er.append(
            pd.Series([i] + list(f_scores), index=['p_HRD', 'recall', 'precision', 'F-score']),
            ignore_index=True
        )
    
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='p_HRD', y='value', hue='variable', 
                 data=pd.melt(hrd_thresholds_er, ['p_HRD']))
    plt.axvline(x=hrd_thresholds_er.loc[hrd_thresholds_er['F-score'].idxmax(), 'p_HRD'], 
                color='red', linestyle='--')
    plt.xlabel('p(HRD)')
    plt.ylabel('Score')
    plt.title(f'Performance Metrics vs HRD Probability Threshold (ER-{er_status})')
    plt.savefig(f'~/Projects/HRD_MutationalSignature/Figures/Supp_TCGAoptimalHRDthreshold_{er_status}.pdf')
    plt.close()