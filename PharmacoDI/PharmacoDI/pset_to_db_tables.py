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

    ## ---- dataset
    dataset = pd.DataFrame.from_dict({
        "id": [np.where(canonical_names.str.match(name))[0][0] + 1], # Wrap single values in list to make 1 row DF
        "name": [pset.get('annotation').get('name')[0]],
    })

    ## ---- tissue
    tissue = pd.DataFrame({
        'id': np.arange(len(pset.get('cell')['tissueid'].unique())) + 1,
        'name': pset.get('cell')['tissueid'].unique(),
    })

    cell_info = pset.get('cell')
    tissue_id_map = dict(zip(tissue.name, tissue.id))
    cell_info['tissue_id'] = [tissue_id_map[tissue_id] for tissue_id in cell_info.tissueid.to_numpy()] # .to_numpy() should speed up iteration
    cell = cell_info[['cellid', 'tissue_id']].copy()
    cell.insert(0, 'id', np.arange(1, len(cell_info.cellid) + 1))
    cell.columns = ['id', 'name', 'tissue_id']

    ## ---- compound
    compound = pd.DataFrame({
        "id": 1,
        "name": []
    })

    ## gene
    gene = pd.DataFrame({
        "id": 1,
        "name": []
    })

    ## target
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

    ## ---- cellosaurus
    cellosaurus = pd.DataFrame({
        'id': np.arange(1, len(cell.id) + 1), ## FIXME: Do can't we just use cell_id as PK?
        'cell_id': cell.id,
        'identifier': cell_info['COSMIC.identifier'],
        'accession': cell_info['Cellosaurus.Accession.id'],
        'as': np.repeat(None, len(cell.id)),
        'sy': np.repeat(None, len(cell.id)),
        'dr': np.repeat(None, len(cell.id)),
        'rx': np.repeat(None, len(cell.id)),
        'ww': np.repeat(None, len(cell.id)),
        'cc': np.repeat(None, len(cell.id)),
        'st': np.repeat(None, len(cell.id)),
        'di': np.repeat(None, len(cell.id)),
        'ox': np.repeat(None, len(cell.id)),
        'hi': np.repeat(None, len(cell.id)),
        'oi': np.repeat(None, len(cell.id)),
        'sx': np.repeat(None, len(cell.id)),
        'ca': np.repeat(None, len(cell.id)),
    })


    # ---- Join tables ----

    # ---- Write to .csv ----
