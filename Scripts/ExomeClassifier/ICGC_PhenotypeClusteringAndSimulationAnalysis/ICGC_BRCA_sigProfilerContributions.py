import pandas as pd
import numpy as np


# Load SBS data
sbs_pcawg = pd.read_csv("PCAWG_sigProfiler_SBS_signatures_in_samples.csv", index_col=1)
sbs_nopcawg = pd.read_csv(
    "nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv", index_col=1
)

sbs = pd.concat([sbs_pcawg, sbs_nopcawg])
sbs_brca = sbs[sbs["Cancer.Types"].str.contains("Breast")].iloc[:, 2:]

sbs_brca_props = pd.DataFrame(
    {"Sigs": sbs_brca.columns, "Prop_present": (sbs_brca > 0).mean()}
)

# Load ID data
id_data = pd.read_csv("PCAWG_SigProfiler_ID_signatures_in_samples.csv", index_col=1)
id_brca = id_data[id_data["Cancer.Types"].str.contains("Breast")].iloc[:, 2:]

id_brca_props = pd.DataFrame(
    {"Sigs": id_brca.columns, "Prop_present": (id_brca > 0).mean()}
)

# Collate and save
brca_sigs = pd.concat([sbs_brca_props, id_brca_props])
brca_sigs.to_csv("ICGC_BRCA_sigProfilerCont.txt", sep="\t", index=False)


# make sure files are being saved in correct directory.
