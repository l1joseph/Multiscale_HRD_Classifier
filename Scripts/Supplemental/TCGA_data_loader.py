import requests
import json
import pandas as pd
import numpy as np

def load_tcga_data(project, data_category, data_type, workflow_type, barcodes):
    """
    Load TCGA data using the GDC API.
    
    :param project: TCGA project (e.g., 'TCGA-BRCA')
    :param data_category: Data category (e.g., 'Transcriptome Profiling')
    :param data_type: Data type (e.g., 'Gene Expression Quantification')
    :param workflow_type: Workflow type (e.g., 'STAR - Counts')
    :param barcodes: List of sample barcodes
    :return: Processed expression data
    """
    
    # Construct the query
    filters = {
        "op": "and",
        "content":[
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "files.data_category", "value": [data_category]}},
            {"op": "in", "content": {"field": "files.data_type", "value": [data_type]}},
            {"op": "in", "content": {"field": "files.analysis.workflow_type", "value": [workflow_type]}},
            {"op": "in", "content": {"field": "cases.samples.portions.analytes.aliquots.submitter_id", "value": barcodes}}
        ]
    }

    params = {
        "filters": json.dumps(filters),
        "fields": "file_id",
        "format": "JSON",
        "size": "100"
    }

    # Send the query
    response = requests.get("https://api.gdc.cancer.gov/files", params=params)
    file_ids = [file['file_id'] for file in response.json()['data']['hits']]

    # Download the files
    params = {"ids": file_ids}
    response = requests.post("https://api.gdc.cancer.gov/data", data=json.dumps(params), headers={"Content-Type": "application/json"})

    # Process the data
    expression_data = []
    for file_content in response.iter_content(chunk_size=None):
        df = pd.read_csv(file_content, sep='\t', index_col=0)
        expression_data.append(df)

    # Combine all files
    expr_test = pd.concat(expression_data, axis=1)

    # Filter for primary tumor samples
    expr_test = expr_test.loc[:, expr_test.columns.str.contains('-01A-')]

    # Remove duplicate gene names and NA values
    expr_test = expr_test.loc[~expr_test.index.duplicated(keep='first')]
    expr_test = expr_test.dropna()

    return expr_test

# Usage:
Z_tumor_testing = pd.read_pickle('~/Data/TCGA/TCGA_BRCA.ExprDeconvolution_050823_p0.79_testing.pkl')
barcodes = Z_tumor_testing.index.str[:12].tolist()

expr_test = load_tcga_data(
    project='TCGA-BRCA',
    data_category='Transcriptome Profiling',
    data_type='Gene Expression Quantification',
    workflow_type='STAR - Counts',
    barcodes=barcodes
)

# Further processing as needed
expr_tumor_testing = expr_test.T  # transpose to match the R script
expr_tumor_testing.index = expr_tumor_testing.index.str[:12]  # keep only the first 12 characters of the barcode




# notes:


# This Python function load_tcga_data does the following:

# Constructs a query to the GDC API based on the provided parameters.
# Sends the query to get the file IDs.
# Downloads the files using the obtained file IDs.
# Processes each file and combines them into a single DataFrame.
# Filters for primary tumor samples and removes duplicates and NA values.

# To use this function, you would replace the lines that load expr_tumor_testing with this new function call.
# We might want to add error handling and potentially a caching mechanism to avoid re-downloading data unnecessarily.
# Also this might be subject to rate limiting or other restrictions from  GDC API.