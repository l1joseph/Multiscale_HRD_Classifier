import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
import joblib
import os

# Set working directory
os.chdir("~/TranscriptionalSignatures/Revisions")

# Load in training data
Z_tumor_training = pd.read_pickle(
    "Data/TCGA_BRCA.ExprDeconvolution_270723_p0.79_training.pkl"
)
Z_tumor_training.index = Z_tumor_training.index.str[:12]

# Organize output variable
ann_tcga = pd.read_pickle("Data/TCGA_HRDclassification_BRCAannotation_HRD0.79.pkl")
ann_tcga = ann_tcga.set_index("Patient")
ann_tcga = ann_tcga[["HRD", "BRCA_status"]]
ann_tcga = ann_tcga[~ann_tcga["BRCA_status"].isna()]
ann_tcga = ann_tcga[~ann_tcga["BRCA_status"].isin(["PALB2", "RAD51C"])]

# Match output to training data
samples_intersect = Z_tumor_training.index.intersection(ann_tcga.index)

ann_tcga = ann_tcga.loc[samples_intersect]
ann_tcga["group"] = ann_tcga["BRCA_status"]
ann_tcga.loc[
    (ann_tcga["HRD"] == "HRD") & (ann_tcga["BRCA_status"] == "none"), "group"
] = "HRD_BRCA+"
ann_tcga.loc[ann_tcga["group"] == "none", "group"] = "HR-proficient"
ann_tcga["group"] = pd.Categorical(
    ann_tcga["group"], categories=["HR-proficient", "HRD_BRCA+", "BRCA1", "BRCA2"]
)

Z_tumor_training = Z_tumor_training.loc[samples_intersect]

# Separate out expression data from BRCA_defect status
input_x = Z_tumor_training
input_y = ann_tcga["group"]


# Calculate weightings
def weight_function(index, vector):
    return 1 / sum(vector == vector[index])


weights_input = np.array([weight_function(i, input_y) for i in range(len(input_y))])

# Conduct 100 iterations of 10-fold cross validation
iterations = range(1, 101)
coefs = []

for i in iterations:
    # Set seed with each iteration
    np.random.seed(123 * i)

    print(f"Iteration {i} of {iterations[-1]}...{pd.Timestamp.now()}")

    # Conduct grouped multinomial 10-fold cross validation
    kf = KFold(n_splits=10, shuffle=True, random_state=123 * i)

    fold_coefs = []
    for train_index, val_index in kf.split(input_x):
        X_train, X_val = input_x.iloc[train_index], input_x.iloc[val_index]
        y_train, y_val = input_y.iloc[train_index], input_y.iloc[val_index]
        weights_train = weights_input[train_index]

        # Standardize the features
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)

        # Fit the model
        model = LogisticRegression(
            penalty="elasticnet",
            solver="saga",
            l1_ratio=0.25,
            max_iter=1000,
            multi_class="multinomial",
            random_state=123 * i,
        )
        model.fit(X_train_scaled, y_train, sample_weight=weights_train)

        fold_coefs.append(model.coef_)

    coefs.append(np.mean(fold_coefs, axis=0))

# Save coefficients
joblib.dump(
    coefs, "Results/CVcoefficients_MultiElasticNet_alpha0.25_p0.79_iter1to100.joblib"
)
