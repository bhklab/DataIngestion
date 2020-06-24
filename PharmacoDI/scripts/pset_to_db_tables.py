import pandas as pd
import numpy as np

def pset_to_db_tables(pset, save_dir, api_url= "https://www.orcestra.ca/api/psets/available"):
    """
    Take in a Python dictionary of a PSet and convert it to database tables for PharmacoDB, with a .csv file for
    each table saved to `save_dir`.

    :param pset: [dict] a nested dictionary containing data from a PSet, as returned by the `convert_pset_to_py`
        function.
    :param save_dir: [string] path to the directory where the database table .csv files should be written
    :param api_url: [string] URL to fetch available PSets from Orecstra, defaults to current URL
    :return: [None] writes .csv files to `save_dir`
    """
    available_psets = pd.read_json(api_url)
    canonical_psets = available_psets.loc[available_psets.canonical, 'URL']

    # ---- Primary tables ----
    dataset = pd.DataFrame({
        ## TODO:: Determine dataset id from availablePSets df
        "id": np.where(pset.get('annotation').get('name')[0] in canonical_psets),
        "name": pset.get('annotation').get('name')[0]
    })

    tissue = pd.DataFrame({
        'id': 1,
        'name': [],
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
