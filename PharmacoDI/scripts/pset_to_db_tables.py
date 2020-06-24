import pandas as pd
import numpy as np
import re

def pset_to_db_tables(pset, save_dir, api_url= "https://www.orcestra.ca/api/psets/canonical"):
    """
    Take in a Python dictionary of a PSet and convert it to database tables for PharmacoDB, with a .csv file for
    each table saved to `save_dir`.

    :param pset: [dict] a nested dictionary containing data from a PSet, as returned by the `convert_pset_to_py`
        function.
    :param save_dir: [string] path to the directory where the database table .csv files should be written
    :param api_url: [string] URL to fetch available PSets from Orecstra, defaults to current URL
    :return: [None] writes .csv files to `save_dir`
    """
    canonical_names = pd.read_json(api_url).name
    name = re.sub('_', '.*', pset.get('annotation').get('name')[0])

    # ---- Primary tables ----
    dataset = pd.DataFrame.from_dict({
        ## TODO:: Determine dataset id from availablePSets df
        "id": [np.where(canonical_names.str.match(name))[0][0] + 1], # Wrap single values in list to make 1 row DF
        "name": [pset.get('annotation').get('name')[0]],
    })

    tissue = pd.DataFrame({
        'id': np.arange(len(pset.get('cell')['tissueid'].unique())) + 1,
        'name': pset.get('cell')['tissueid'].unique(),
    })

    cell = pd.DataFrame({
        "id": 1,
        "name": [],
        "tissue_id": [],
    })

    compound = pd.DataFrame({
        "id": 1,
        "name": []
    })

    gene = pd.DataFrame({
        "id": 1,
        "name": []
    })

    target = pd.DataFrame({
        "id": 1,
        "name": [],
        "gene_id": []
    })


    # ---- Annotation tables ----
    compound_annotation = pd.DataFrame({
        "drug_id": 1,
        "name": []
    })

    compound_target = pd.DataFrame({
        "id": 1,
        "name": []
    })

    gene_annotations = pd.DataFrame({
        "id": 1,
        "name": []
    })

    cellosaurus = pd.DataFrame({
        "id": 1,
        "name": []
    })

    dataset_statistics = pd.DataFrame({
        "id": 1,
        "name": []
    })


    # ---- Derived tables ----



    # ---- Join tables ----

    # ---- Write to .csv ----
