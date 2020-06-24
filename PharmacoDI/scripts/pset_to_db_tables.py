import pandas as pd
import numpy as np

def pset_to_db_tables(pset, save_dir):

    # ---- Primary tables ----
    dataset = pd.DataFrame({
        ## TODO:: Determine dataset id from availablePSets df
        "id": 1,
        "name": pset.get('annotation').get('name')
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
        "name": []
    })


    # ---- Annotation tables ----
    compound_annotation = pd.DataFrame({
        "id": 1,
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
