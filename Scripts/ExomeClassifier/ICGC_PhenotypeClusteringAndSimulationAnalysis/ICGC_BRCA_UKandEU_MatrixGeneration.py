import pandas as pd
import numpy as np
from scipy import stats
import pyranges as pr
from pyfaidx import Fasta

# Load libraries equivalent to maftools and sigminer, maybe using R snippet.
#or use the following code to load the libraries

def read_maf(maf_file):
    """Read MAF file and return a DataFrame"""
    return pd.read_csv(maf_file, sep='\t', low_memory=False)

def sig_tally(maf_df, ref_genome, mode='ALL', use_syn=True):
    """
    Generate mutation tally similar to sigminer's sig_tally function
    This is a simplified version and may need to be expanded based on specific requirements
    """
    # Load reference genome
    genome = Fasta(ref_genome)
    
    # Generate SBS96 context
    def get_context(chrom, pos, ref, alt):
        seq = genome[chrom][pos-2:pos+1].seq
        return f"{seq[0]}[{ref}>{alt}]{seq[2]}"
    
    maf_df['Context'] = maf_df.apply(lambda row: get_context(row['Chromosome'], 
                                                             row['Start_position'], 
                                                             row['Reference_Allele'], 
                                                             row['Tumor_Seq_Allele2']), axis=1)
    
    # Generate ID83 context (simplified version)
    def get_indel_context(ref, alt):
        if len(ref) > len(alt):
            return f"DEL:{len(ref) - len(alt)}:1"
        elif len(alt) > len(ref):
            return f"INS:{len(alt) - len(ref)}:1"
        else:
            return "OTHER"
    
    maf_df['IndelContext'] = maf_df.apply(lambda row: get_indel_context(row['Reference_Allele'], 
                                                                        row['Tumor_Seq_Allele2']), axis=1)
    
    # Tally mutations
    sbs_tally = maf_df['Context'].value_counts()
    id_tally = maf_df['IndelContext'].value_counts()
    
    return {'SBS_96': sbs_tally, 'ID_83': id_tally}

# Load UK and EU BRCA data
brca_uk = pd.read_csv('~/Data/ICGC/simple_somatic_mutation.open.BRCA-UK.tsv.gz', sep='\t')
brca_eu = pd.read_csv('~/Data/ICGC/simple_somatic_mutation.open.BRCA-EU.tsv.gz', sep='\t')

# Extract WGS strategies, remove duplicated mutations
brca_uk_wgs = brca_uk[brca_uk['sequencing_strategy'] == 'WGS'].drop_duplicates(subset='icgc_mutation_id')
brca_eu_wgs = brca_eu[brca_eu['sequencing_strategy'] == 'WGS'].drop_duplicates(subset='icgc_mutation_id')

# Collate ICGC data and organize into MAF-readable format
brca_wgs = pd.concat([brca_uk_wgs, brca_eu_wgs])

brca_wgs_input = pd.DataFrame({
    'Tumor_Sample_Barcode': brca_wgs['icgc_donor_id'],
    'Hugo_Symbol': np.nan,
    'Chromosome': brca_wgs['chromosome'],
    'Start_position': brca_wgs['chromosome_start'],
    'End_position': brca_wgs['chromosome_end'],
    'Variant_Classification': brca_wgs['consequence_type'],
    'Variant_Type': brca_wgs['mutation_type'].map({
        'single base substitution': 'SNP',
        'insertion of <=200bp': 'INS',
        'deletion of <=200bp': 'DEL'
    }),
    'Reference_Allele': brca_wgs['reference_genome_allele'],
    'Tumor_Seq_Allele2': brca_wgs['mutated_to_allele']
})

# Create MAF object
brca_maf = read_maf(brca_wgs_input)

# Run mutation tally
mt_tally_brca_wgs = sig_tally(
    brca_maf,
    ref_genome='path/to/BSgenome.Hsapiens.UCSC.hg19.fasta',
    mode='ALL',
    use_syn=True
)

# Save the results
import pickle
with open('~/Projects/HRD_MutationalSignature/Data/BRCA_UKEU_mt_tally.pkl', 'wb') as f:
    pickle.dump(mt_tally_brca_wgs, f)