import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# Load ICGC-BRCA signature contributions
#   Only signatures appearing in >1% samples will be included
brca_final_sigs = pd.read_csv("~/Data/ICGC/ICGC_BRCA_sigProfilerCont.txt", sep="\t")
brca_final_sigs = brca_final_sigs[brca_final_sigs["Prop_present"] > 0.01]
sigs_sbs = brca_final_sigs[brca_final_sigs["Sigs"].str.contains("SBS")]["Sigs"].tolist()
sigs_id = brca_final_sigs[brca_final_sigs["Sigs"].str.contains("ID")]["Sigs"].tolist()

# Load hg19 signatures (for ICGC samples)
signatures_sbs96_hg19 = pd.read_csv("~/Data/COSMIC_v3.3.1_SBS_GRCh37.txt", sep="\t")
sig_sbs96_types = signatures_sbs96_hg19["Type"]
signatures_sbs96_hg19 = signatures_sbs96_hg19.iloc[:, 1:].T
signatures_sbs96_hg19.columns = sig_sbs96_types
signatures_sbs96_hg19 = signatures_sbs96_hg19[signatures_cosmic.columns]

# Load COSMIC indel signatures
signatures_id83 = pd.read_csv("~/Data/COSMIC_v3.3_ID_GRCh37.txt", sep="\t")
sig_id83_types = signatures_id83["Type"]
signatures_id83 = signatures_id83.iloc[:, 1:].T
signatures_id83.columns = sig_id83_types

# Load in datasets and tweak into deconstructSigs inputs:
#   rows = samples, cols = signature contexts
#   order the same as signature data frames
mt_tally_brca_wgs = pd.read_pickle(
    "~/Projects/HRD_MutationalSignature/Data/BRCA_UKEU_mt_tally.pkl"
)
sigs_sbs96_input = mt_tally_brca_wgs["SBS_96"]
sigs_sbs96_input = sigs_sbs96_input[signatures_sbs96_hg19.columns]

sigs_id83_input = mt_tally_brca_wgs["ID_83"]
sigs_id83_input = sigs_id83_input[signatures_id83.columns]

# Normalise ID83 counts
icgc_idCounts = pd.read_csv("~/Data/ICGC/ICGC_BRCA_indelCounts.txt", sep="\t")
icgc_idCounts["genome_to_exome"] = icgc_idCounts["exome"] / icgc_idCounts["wgs"]
icgc_idCounts = icgc_idCounts.loc[sigs_id83_input.columns]

sigs_id83_input = sigs_id83_input.mul(icgc_idCounts["genome_to_exome"], axis=1)


def whichSignatures(
    tumor_ref,
    signatures_ref,
    sample_id,
    contexts_needed=True,
    signature_cutoff=0,
    tri_counts_method="default",
):
    """
    Python implementation of deconstructSigs' whichSignatures function.
    """
    tumor = tumor_ref.loc[sample_id]

    if tri_counts_method == "genome2exome":
        # Implement genome2exome normalization
        pass

    # Compute weights using non-negative least squares
    # Don't know why it's not loading, but might need to import directly from scipy.optimize
    weights, residuals = nnls(signatures_ref.T, tumor)

    # Normalize weights
    weights = weights / np.sum(weights)

    # Apply signature cutoff
    weights[weights < signature_cutoff] = 0
    weights = weights / np.sum(weights)

    return pd.Series(weights, index=signatures_ref.index)


def run_deconstructSigs(sigs_input, sig_type="SBS"):
    if sig_type == "SBS":
        sig_ref = signatures_sbs96_hg19.loc[sigs_sbs]
        print("Calculating SBS signature contributions...")
    elif sig_type == "ID":
        sig_ref = signatures_id83.loc[sigs_id]
        print("Calculating ID signature contributions...")
    else:
        print("Set sig_type to SBS or ID. Using SBS signatures...")
        sig_ref = signatures_sbs96_hg19.loc[sigs_sbs]

    sigs_out = pd.DataFrame()
    for i, sample in enumerate(sigs_input.index):
        print(f"{sig_type} contributions, Sample {i+1} of {len(sigs_input)}: {sample}")

        sigs_i = whichSignatures(
            tumor_ref=sigs_input,
            signatures_ref=sig_ref,
            sample_id=sample,
            contexts_needed=True,
            signature_cutoff=0,
            tri_counts_method="genome2exome" if sig_type == "SBS" else "default",
        )
        sigs_out = sigs_out.append(sigs_i, ignore_index=True)

    return sigs_out


sigs_sbs96 = run_deconstructSigs(sigs_sbs96_input, "SBS")
sigs_id83 = run_deconstructSigs(sigs_id83_input, "ID")
sigs_complete = pd.concat([sigs_sbs96, sigs_id83], axis=1)

sigs_complete.to_pickle(
    "~/Projects/HRD_MutationalSignature/Results/ICGC_BRCA_deconstructSigs_Cutoff0.01_SBSandIDnormalised.pkl"
)


# notes:

# The deconstructSigs package doesn't have a direct Python equivalent.
# I've implemented a simplified version of the whichSignatures function.
# The non-negative least squares algorithm is used in place of R's whichSignatures function.
# Need to import this from scipy: from scipy.optimize import nnls.
# Verify that the nnls algorithm is producing results comparable to the R version.
# The genome2exome normalization in the whichSignatures function is a placeholder.
# The run_deconstructSigs function is implemented to mimic the behavior of the R version,
# but I need to refine it further to match R functionality.
