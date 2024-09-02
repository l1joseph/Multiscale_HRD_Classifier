import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc

def f_score(confusion_matrix, weight=1):
    tn, fp, fn, tp = confusion_matrix.ravel()
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)
    f_score = (1 + weight**2) * (precision * recall) / ((weight**2 * precision) + recall)
    return pd.Series({'recall': recall, 'precision': precision, 'F-score': f_score})

# Set working directory
os.chdir('~/Projects/HRD_MutationalSignature/Results')

# Load TCGA assignments with BRCA-defect and ER status labels
ann_tcga = pd.read_pickle('TCGA_HRDclassification_BRCAannotation.pkl')

# Calculate F-scores for different thresholds
hrd_thresholds = []
for i in np.arange(0, 1.01, 0.01):
    tcga_results = ann_tcga.copy()
    tcga_results['HR_defect'] = np.where(tcga_results['BRCA_status'] != 'none', 'defective', 'proficient')
    tcga_results['HRD'] = np.where(tcga_results['HRD_prob'] > i, 'HRD', 'HR-proficient')
    
    confusion_matrix = pd.crosstab(tcga_results['HRD'], tcga_results['HR_defect'])
    f_scores = f_score(confusion_matrix)
    hrd_thresholds.append([i] + f_scores.tolist())

hrd_thresholds = pd.DataFrame(hrd_thresholds, columns=['p_HRD', 'recall', 'precision', 'F-score'])

# Plot results
plt.figure(figsize=(10, 6))
for measure in ['recall', 'precision', 'F-score']:
    plt.plot(hrd_thresholds['p_HRD'], hrd_thresholds[measure], label=measure)
plt.axvline(x=hrd_thresholds.loc[hrd_thresholds['F-score'].idxmax(), 'p_HRD'], color='r', linestyle='--')
plt.xlabel('p(HRD)')
plt.ylabel('Score')
plt.legend()
plt.title('HRD Classification Thresholds')
plt.savefig('../Figures/Supp_TCGAoptimalHRDthreshold.pdf')
plt.close()

# Plot AUC curve
ann_tcga_auc = ann_tcga.dropna(subset=['BRCA_status'])
ann_tcga_auc['HR_geneDefective'] = (ann_tcga_auc['BRCA_status'] != 'none').astype(int)

fpr, tpr, _ = roc_curve(ann_tcga_auc['HR_geneDefective'], ann_tcga_auc['HRD_prob'])
roc_auc = auc(fpr, tpr)

plt.figure(figsize=(8, 8))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.savefig('../Figures/Supp_TCGAHRDthresholdAUC.pdf')
plt.close()

# Repeat analysis for ER-positive and ER-negative subsets
for er_sub in ['Positive', 'Negative']:
    ann_tcga_er = ann_tcga[ann_tcga['ER_status'] == er_sub]
    
    hrd_thresholds_er = []
    for i in np.arange(0, 1.01, 0.01):
        tcga_results = ann_tcga_er.copy()
        tcga_results['HR_defect'] = np.where(tcga_results['BRCA_status'] != 'none', 'defective', 'proficient')
        tcga_results['HRD'] = np.where(tcga_results['HRD_prob'] > i, 'HRD', 'HR-proficient')
        
        confusion_matrix = pd.crosstab(tcga_results['HRD'], tcga_results['HR_defect'])
        f_scores = f_score(confusion_matrix)
        hrd_thresholds_er.append([i] + f_scores.tolist())

    hrd_thresholds_er = pd.DataFrame(hrd_thresholds_er, columns=['p_HRD', 'recall', 'precision', 'F-score'])

    plt.figure(figsize=(10, 6))
    for measure in ['recall', 'precision', 'F-score']:
        plt.plot(hrd_thresholds_er['p_HRD'], hrd_thresholds_er[measure], label=measure)
    plt.axvline(x=hrd_thresholds_er.loc[hrd_thresholds_er['F-score'].idxmax(), 'p_HRD'], color='r', linestyle='--')
    plt.xlabel('p(HRD)')
    plt.ylabel('Score')
    plt.legend()
    plt.title(f'HRD Classification Thresholds (ER-{er_sub})')
    plt.savefig(f'../Figures/Supp_TCGAoptimalHRDthreshold_{er_sub}.pdf')
    plt.close()
