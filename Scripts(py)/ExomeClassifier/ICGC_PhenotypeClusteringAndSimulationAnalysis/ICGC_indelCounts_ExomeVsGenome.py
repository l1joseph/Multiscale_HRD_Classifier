import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Load libraries
# Note: We might need to install biopython and pysam for maftools and sigminer functionality
# Could also just run some of this in R, using the rpy2 package.

from Bio import SeqIO
import pysam
import maftools
import sigminer

# Load UK and EU BRCA data (US is exome, FR contains no indels)
# Extract WGS strategies, remove duplicated mutations

brca_uk = pd.read_csv(
    "~/Data/ICGC/simple_somatic_mutation.open.BRCA-UK.tsv.gz", sep="\t"
)
brca_uk_wgs = brca_uk[brca_uk["sequencing_strategy"] == "WGS"]
brca_uk_wgs = brca_uk_wgs.drop_duplicates(subset="icgc_mutation_id")

brca_eu = pd.read_csv(
    "~/Data/ICGC/simple_somatic_mutation.open.BRCA-EU.tsv.gz", sep="\t"
)
brca_eu_wgs = brca_eu[brca_eu["sequencing_strategy"] == "WGS"]
brca_eu_wgs = brca_eu_wgs.drop_duplicates(subset="icgc_mutation_id")

# Collate ICGC data and extract mutations in intron or intergenic regions

brca_wgs = pd.concat([brca_uk_wgs, brca_eu_wgs])
brca_wgs_input = pd.DataFrame(
    {
        "Tumor_Sample_Barcode": brca_wgs["icgc_donor_id"],
        "Hugo_Symbol": np.nan,
        "Chromosome": brca_wgs["chromosome"],
        "Start_position": brca_wgs["chromosome_start"],
        "End_position": brca_wgs["chromosome_end"],
        "Variant_Classification": brca_wgs["consequence_type"],
        "Variant_Type": brca_wgs["mutation_type"].map(
            {
                "single base substitution": "SNP",
                "insertion of <=200bp": "INS",
                "deletion of <=200bp": "DEL",
            }
        ),
        "Reference_Allele": brca_wgs["reference_genome_allele"],
        "Tumor_Seq_Allele2": brca_wgs["mutated_to_allele"],
    }
)

brca_exome_input = brca_wgs_input[
    brca_wgs_input["Variant_Classification"] != "intergenic_region"
]

# Read MAF files
# Note: We might need to implement a custom read_maf function if maftools alternative is not available in Python
brca_wgs_maf = maftools.read_maf(brca_wgs_input)
brca_exome_maf = maftools.read_maf(brca_exome_input)

# Calculate ID-83 counts for each MAF file
# Note: We might need to implement a custom sig_tally function if sigminer alternative is not available in Python
mt_tally_brca_wgs = sigminer.sig_tally(
    brca_wgs_maf, ref_genome="BSgenome.Hsapiens.UCSC.hg19", mode="ID", use_syn=True
)

mt_tally_brca_exome = sigminer.sig_tally(
    brca_exome_maf, ref_genome="BSgenome.Hsapiens.UCSC.hg19", mode="ID", use_syn=True
)

# Count total mutations of each ID-83 type
icgc_indelCounts = pd.DataFrame(
    {
        "wgs": mt_tally_brca_wgs["all_matrices"]["ID_83"].sum(),
        "exome": mt_tally_brca_exome["all_matrices"]["ID_83"].sum(),
    }
)

icgc_indelCounts.to_csv("ICGC_BRCA_indelCounts.txt", sep="\t")

icgc_indelCounts["ex_gen"] = icgc_indelCounts["exome"] / icgc_indelCounts["wgs"]

plt.figure(figsize=(10, 6))
plt.hist(icgc_indelCounts["ex_gen"], bins=30)
plt.title("Histogram of exome/genome ratio")
plt.xlabel("exome/genome ratio")
plt.ylabel("Frequency")
plt.savefig("exome_genome_ratio_histogram.png")
plt.close()

plt.figure(figsize=(10, 6))
plt.scatter(np.log10(icgc_indelCounts["wgs"]), np.log10(icgc_indelCounts["exome"]))
plt.xlabel("log10(wgs)")
plt.ylabel("log10(exome)")
plt.title("WGS vs Exome Indel Counts")

# Add regression line
slope, intercept, r_value, p_value, std_err = stats.linregress(
    np.log10(icgc_indelCounts["wgs"]), np.log10(icgc_indelCounts["exome"])
)
line = slope * np.log10(icgc_indelCounts["wgs"]) + intercept
plt.plot(
    np.log10(icgc_indelCounts["wgs"]),
    line,
    "r",
    label=f"y={slope:.2f}x+{intercept:.2f}",
)
plt.legend()

plt.savefig("wgs_vs_exome_scatter.png")
plt.close()


# notes:

# need to implement custom functions to replace read_maf and sig_tally.
# I haven't checked if BSgenome.Hsapiens.UCSC.hg19 reference genome is available,
# so we might need to download this.

