import os
import PharmacoDI as di
import pandas as pd
import numpy as np
import urllib
import requests


def get_target_annotations(pset, annot_dir):
    # Read in drug target annotations and gene annotations
    drug_targets = pd.read_csv(os.path.join(annot_dir, 'drugbank_drug_targets_all.csv'))
    rnaseq_df = pset.get("molecularProfiles").get("Kallisto_0.46.1.rnaseq").get("elementMetadata")

    # Map genes to drugbank drug ids
    genes_to_drugs = pd.merge(drug_targets.loc[:, ['Name', 'Gene Name', 'Drug IDs']],
                              rnaseq_df.loc[:, ['gene_name', 'gene_id']],
                              left_on='Gene Name', right_on='gene_name')
    genes_to_drugs['Drug IDs'] = [str.split(ids, '; ') for ids in genes_to_drugs['Drug IDs'].values]
    genes_to_drugs = genes_to_drugs.explode('Drug IDs')

    file_path = os.path.join(annot_dir, 'drugbank_drug_to_gene_mapplings.csv')
    if not os.isfile(file_path):
        pd.write_csv(genes_to_drugs, file_path)


def query_uniprot_mapping_api(ids):
    url = 'https://www.uniprot.org/uploadlists/'

    # Build the query
    query = ' '.join([str(id) for id in ids])

    params = {
        'from': 'ACC+ID',
        'to': 'ENSEMBL_ID',
        'format': 'tab',
        'query': query
    }

    query_string = str(urllib.parse.urlencode(params))
    req = requests.get(f'{url}{query_string}')
