import pandas as pd
import numpy as np
import glob
import os

# Load deconvoluted training data for reference
Z_tumor_training = pd.read_pickle(
    "~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.pkl"
)
Z_tumor_training.index = Z_tumor_training.index.str[:12]
Z_tumor_training = np.log2(Z_tumor_training + 1)
samples = Z_tumor_training.index
Z_tumor_training = (Z_tumor_training - Z_tumor_training.mean()) / Z_tumor_training.std()
Z_tumor_training = pd.DataFrame(Z_tumor_training, index=samples)

ann_tcga = pd.read_pickle(
    "~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.pkl"
)
ann_tcga = ann_tcga.set_index("Patient")
ann_tcga = ann_tcga[["HRD", "BRCA_status"]]
ann_tcga = ann_tcga[~ann_tcga["BRCA_status"].isna()]
ann_tcga = ann_tcga[~ann_tcga["BRCA_status"].isin(["PALB2", "RAD51C"])]

ann_tcga["group"] = ann_tcga["BRCA_status"]
ann_tcga.loc[
    (ann_tcga["HRD"] == "HRD") & (ann_tcga["BRCA_status"] == "none"), "group"
] = "HRD_BRCA+"
ann_tcga.loc[ann_tcga["group"] == "none", "group"] = "HR-proficient"
ann_tcga["group"] = pd.Categorical(
    ann_tcga["group"], categories=["HR-proficient", "HRD_BRCA+", "BRCA1", "BRCA2"]
)

samples_intersect = Z_tumor_training.index.intersection(ann_tcga.index)

ann_tcga_train = ann_tcga.loc[samples_intersect]
Z_tumor_training = Z_tumor_training.loc[samples_intersect]

# Move to CV coefficients and generate centroids
os.chdir("~/Projects/HRD_TranscriptionalSignature/CV_Coefficients/Run6/")
signature_centroid_list = {}

models = set("_".join(file.split("_")[1:3]) for file in glob.glob("*"))

for model in models:
    print(model)

    files_coefs = glob.glob(f"*{model}*")
    coefs_join = []
    for file in files_coefs:
        coefs = pd.read_pickle(file)
        coefs_join.extend(coefs)

    coef_mat = pd.DataFrame(
        [coef["HR-proficient"].iloc[:, 0] != 0 for coef in coefs_join]
    ).T
    coef_mat = coef_mat.iloc[1:]  # Remove intercept
    genes_include = coef_mat.index[coef_mat.sum(axis=1) == 1000].tolist()

    print(f"Number of genes in model {model}: {len(genes_include)}")

    ## Create centroids
    centroid_model = pd.DataFrame(
        {
            "HRD": Z_tumor_training.loc[
                ann_tcga_train["HRD"] == "HRD", genes_include
            ].mean(),
            "HR_proficient": Z_tumor_training.loc[
                ann_tcga_train["HRD"] == "HR-proficient", genes_include
            ].mean(),
            "BRCA1": Z_tumor_training.loc[
                ann_tcga_train["group"] == "BRCA1", genes_include
            ].mean(),
            "BRCA2": Z_tumor_training.loc[
                ann_tcga_train["group"] == "BRCA2", genes_include
            ].mean(),
            "HRD_BRCApos": Z_tumor_training.loc[
                ann_tcga_train["group"] == "HRD_BRCA+", genes_include
            ].mean(),
            "HR_BRCA_proficient": Z_tumor_training.loc[
                ann_tcga_train["group"] == "HR-proficient", genes_include
            ].mean(),
        }
    )

    signature_centroid_list[model] = centroid_model

pd.to_pickle(
    signature_centroid_list,
    "~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.pkl",
)
