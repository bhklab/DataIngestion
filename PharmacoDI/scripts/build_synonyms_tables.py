import os
import re
import glob
import numpy as np
import pandas as pd
import dask.dataframe as dd
from scripts.join_pset_tables import *

read_file_path = os.path.join('data', 'procdata')
write_file_path = os.path.join('data', 'latest')
metadata_path = os.path.join('data', 'metadata')

cell_meta_file = "cell_annotation_all.csv"
tissue_meta_file = "cell_annotation_all.csv"
drug_meta_file = "drugs_with_ids.csv"


def get_metadata(file_name, metadata_path):
    # Find correct metadata annotations CSV file
    annotations_file = glob.glob(
        os.path.join(metadata_path, file_name))[0]
    if annotations_file is None:
        raise ValueError(
            f'No metadata file named {file_name} could be found in {metadata_path}')
    
    # Read csv file and return df
    return pd.read_csv(annotations_file, index_col=[0])


# --- SYNONYMS TABLES --------------------------------------------------------------------------

# TODO - 3 rows; 954 rows
def build_cell_synonyms_df(cell_df, cell_meta_file, metadata_path, write_file_path):
    # Get metadata file
    cell_metadata = get_metadata(cell_meta_file, metadata_path)

    # Find all columns relevant to cellid
    pattern = re.compile('cellid')
    cell_columns = cell_metadata[[
        col for col in cell_metadata.columns if pattern.search(col)]]

    # Convert cell columns to Dask dataframe to enable merging
    cell_columns = dd.from_pandas(cell_columns, chunksize=500000)

    # Get all unique synonyms and join with cell_df
    cell_synonym_df = melt_and_join(cell_columns, 'unique.cellid', cell_df)
    cell_synonym_df = cell_synonym_df.rename(columns={'id': 'cell_id', 'value': 'cell_name'})

    # Add blank col for dataset_id (TODO)
    cell_synonym_df['dataset_id'] = np.nan

    # Reindex and write to disk
    cell_synonym_df = reindex_table(cell_synonym_df)
    write_table(cell_synonym_df, 'cell_synonym', write_file_path)
    

def build_tissue_synonyms_df(tissue_df, tissue_meta_file, metadata_path, write_file_path):
    # Get metadata file
    tissue_metadata = get_metadata(tissue_meta_file, metadata_path)

    # Find all columns relevant to tissueid
    pattern = re.compile('tissueid')
    tissue_cols = tissue_metadata[[
        col for col in tissue_metadata.columns if pattern.search(col)]]

    # Convert tissue columns to Dask dataframe to enable merging
    tissue_cols = dd.from_pandas(tissue_cols, chunksize=500000)

    # Get all unique synonyms and join with tissue_df
    tissue_synonym_df = melt_and_join(tissue_cols, 'unique.tissueid', tissue_df)
    tissue_synonym_df = tissue_synonym_df.rename(columns={'id': 'tissue_id', 'value': 'tissue_name'})

    # Add blank col for dataset_id (TODO)
    tissue_synonym_df['dataset_id'] = np.nan

    # Reindex and write to disk
    tissue_synonym_df = reindex_table(tissue_synonym_df)
    write_table(tissue_synonym_df, 'tissue_synonym', write_file_path)


def build_drug_synonyms_df(drug_df, drug_meta_file, metadata_path, write_file_path):
    # Get metadata file
    drug_metadata = get_metadata(drug_meta_file, metadata_path)

    # Find all columns relevant to drugid
    # Right now only FDA col is dropped, but may be more metadata in the future
    pattern = re.compile('drugid')
    drug_cols= drug_metadata[[
        col for col in drug_metadata.columns if pattern.search(col)]]

    # Convert drug columns to Dask dataframe to enable merging
    drug_cols = dd.from_pandas(drug_cols, chunksize=500000)

    # Get all unique synonyms and join with drugs_df
    drug_synonym_df = melt_and_join(drug_cols, 'unique.drugid', drug_df)
    drug_synonym_df = drug_synonym_df.rename(columns={'id': 'drug_id', 'value': 'drug_name'})

    # Add blank col for dataset_id (TODO)
    drug_synonym_df['dataset_id'] = np.nan

    # Reindex and write to disk
    drug_synonym_df = reindex_table(drug_synonym_df)
    write_table(drug_synonym_df, 'drug_synonym', write_file_path)


# Helper function for getting all synonyms related to a certain df
def melt_and_join(meta_df, unique_id, join_df):
    """
    @param meta_df: [`Dask DataFrame`] The DataFrame containing all the synonyms (metadata)
    @param unique_id: [`string`] The name of the column in the metadata containing the unique IDs
    @param join_df: [`Dask DataFrame`] THe DataFrame containing the primary keys that will be used as
        foreign keys in the new synonyms df

    @return [`DataFrame`] The synonys dataframe, with a PK, FK based on join_df, and all unique synonyms
    """
    join_df['id'] = join_df.index

    # Convert wide meta_df to long table
    # Drop 'variable' col (leave only unique ID and synonyms), drop duplicates
    synonyms = dd.melt(meta_df, id_vars=[unique_id])[
        [unique_id, 'value']].drop_duplicates()

    # Drop all rows where value is NA
    synonyms = synonyms[synonyms['value'].notnull()]

    # Join with join_df based on unique_id
    synonyms = dd.merge(synonyms, join_df, left_on=unique_id,
                        right_on='name', how='left')[['id', 'value']]
    if synonyms[synonyms['id'].isna()].shape[0].compute() > 0:
        print(f'ERROR - some rows did not join to a {unique_id}!')
    
    synonyms['id'] = synonyms['id'].astype(pd.Int64Dtype())

    return synonyms

