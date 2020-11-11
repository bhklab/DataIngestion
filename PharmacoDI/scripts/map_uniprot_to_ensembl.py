import os
import pandas as pd
import numpy as np
import requests
import json
from multiprocessing import Pool, cpu_count

metadata_path = os.path.join("data", "metadata", "Annotations")

# Import Drugbank data
drugbank_file = "drugbank_targets_has_ref_has_uniprot.csv"
drugbank_df = pd.read_csv(os.path.join(metadata_path, drugbank_file))

# Get UniProtKB IDs
uniprot_ids = pd.Series(
    pd.unique(drugbank_df['polypeptide.external.identifiers.UniProtKB']))

# Split UniProtKB into groups of 1000 entries
queries = [uniprot_ids.iloc[i:i+1000]
           for i in np.arange(0, len(uniprot_ids), 1000)]


def mapUniprotToEnsembl(idList):
    # Function to make GET request and process response
    
    # Make API call
    params = {
        'from': 'ID',
        'to': 'ENSEMBL_ID',
        'format': 'tab',
        'query': " ".join(idList)
    }
    r = requests.get('https://www.uniprot.org/uploadlists/', params=params)
    
    # Check that request was successful
    if r.status_code != 200:
        return f'ERROR: {r.status_code}'
    
    # Split text into rows (one row per ID) and build df
    # Exclude first row (header)
    gene_id_df = pd.DataFrame(r.text.split("\n")[1:])
    
    # Split into two columns, UniprotId and EnsemblId
    gene_id_df = gene_id_df[0].str.split(pat="\t", expand=True)
    
    # Drop empty rows and rename columns
    gene_id_df.dropna(inplace=True)
    gene_id_df.rename(columns={0:'UniProtId', 1: 'EnsemblGeneId'}, inplace=True)

    return gene_id_df


# Start pool to make API calls in parallel
pool = Pool(cpu_count() - 1)
gene_id_df_list = pool.map(mapUniprotToEnsembl, queries)
pool.close()

# Combine all dfs into one and reset index
gene_id_df = pd.concat(gene_id_df_list)
gene_id_df.reset_index(drop=True, inplace=True)
gene_id_df.index = gene_id_df.index + 1

# Write to csv file
mapping_file = "uniprot_ensembl_mappings.csv"
gene_id_df.to_csv(os.path.join(metadata_path, mapping_file), index=True)
