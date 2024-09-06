import os
import pandas as pd
import numpy as np

# Set working directory
os.chdir("/home/zcqsdhj/TranscriptionalSignatures/Revisions/Data")

# Load the data
Z_tumor_training = pd.read_pickle(
    "TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.pkl"
)

# Remove columns (genes) with zero sum
z2 = Z_tumor_training.loc[:, (Z_tumor_training.sum() > 0)]

# Log2 transform the data
z3 = np.log2(z2 + 1)

# Scale each column (gene)
z4 = (z3 - z3.mean()) / z3.std()

# Assign the original index (sample names) to z4
Z_tumor_training = pd.DataFrame(z4, index=z3.index)

# Save the transformed data
Z_tumor_training.to_pickle("TCGA_BRCA.ExprDeconvolution_270723_p0.79_training.pkl")


# notes:
# n R, operations like scale() by default use n-1 in the denominator when calculating standard deviation.
# In Python, I used the numpy implementation which uses n.
# So if we find that this difference is important for our analysis, we might want to use ddof=1 in the std calculation.
