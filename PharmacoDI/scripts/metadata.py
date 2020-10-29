import os
import re
import glob
import numpy as np
import pandas as pd


# --- SYNONYMS TABLES --------------------------------------------------------------------------

annotations_file_path = os.path.join('data', 'metadata', 'Annotations')


def get_annotations(file_name, annotations_file_path):
    # Find correct metadata annotations CSV file
    annotations_file = glob.glob(
        os.path.join(annotations_file_path, file_name))[0]
    if annotations_file is None:
        raise ValueError(
            f'No metadata file named {file_name} could be found in {annotations_file_path}')

    # Read csv file and return df
    return pd.read_csv(annotations_file)


# TODO - 3 rows; 954 rows
def cell_synonyms(cells_df, file_name, file_path):
    # Get metadata file
    file_name = 'cell_annotation_all.csv'
    cell_names_df = get_annotations(file_name, file_path)

    # Find all columns relevant to cellid
    pattern = re.compile('cellid')
    cell_columns = cell_names_df[[
        col for col in cell_names_df.columns if pattern.search(col)]]

    # Get all unique synonyms and join with cells_df
    cell_synonyms = melt_and_join(cell_columns, 'unique.cellid', cells_df)

    # Rename columns
    cell_synonyms.columns = ['cell_id', 'cell_name']

    # Add primary key and blank col for dataset_id (TODO)
    cell_synonyms['id'] = cell_synonyms.index + 1
    cell_synonyms['dataset_id'] = np.nan

    # Reorder columns to match ERD
    return cell_synonyms[['id', 'cell_id', 'dataset_id', 'cell_name']]


def tissue_synonyms(tissues_df, file_name, file_path):
    # Get metadata file
    file_name = 'cell_annotation_all.csv'
    tissues_metadata = get_annotations(file_name, file_path)

    # Find all columns relevant to tissueid
    pattern = re.compile('tissueid')
    tissue_cols = tissues_metadata[[
        col for col in tissues_metadata.columns if pattern.search(col)]]

    # Get all unique synonyms and join with tissues_df
    tissue_synonyms = melt_and_join(tissue_cols, 'unique.tissueid', tissues_df)

    # Rename columns
    tissue_synonyms.columns = ['tissue_id', 'tissue_name']

    # Add primary key and blank col for dataset_id (TODO)
    tissue_synonyms['id'] = tissue_synonyms.index + 1
    tissue_synonyms['dataset_id'] = np.nan

    # Reorder columns to match ERD
    return tissue_synonyms[['id', 'tissue_id', 'dataset_id', 'tissue_name']]


def drug_synonyms(drugs_df, file_name, file_path):
    # Get metadata file
    file_name = 'drugs_with_ids.csv'
    drugs_metadata = get_annotations(file_name, file_path)

    # Find all columns relevant to drugid
    # Right now only FDA col is dropped, but may be more metadata in the future
    pattern = re.compile('drugid')
    drugs_metadata = drugs_metadata[[
        col for col in drugs_metadata.columns if pattern.search(col)]]

    # Get all unique synonyms and join with drugs_df
    drug_synonyms = melt_and_join(drugs_metadata, 'unique.drugid', drugs_df)

    # Rename columns
    drug_synonyms.columns = ['tissue_id', 'tissue_name']

    # Add primary key and blank col for dataset_id (TODO)
    drug_synonyms['id'] = drug_synonyms.index + 1
    drug_synonyms['dataset_id'] = np.nan

    # Reorder columns to match ERD
    return drug_synonyms[['id', 'tissue_id', 'dataset_id', 'tissue_name']]


# Helper function for getting all synonyms related to a certain df
def melt_and_join(meta_df, unique_id, join_df):
    """
    @param meta_df: [`DataFrame`] The dataframe containing all the synonyms (metadata)
    @param unique_id: [`string`] The name of the column in the metadata containing the unique IDs
    @param join_df: [`DataFrame`] THe dataframe containing the primary keys that will be used as
        foreign keys in the new synonyms df

    @return [`DataFrame`] The synonys dataframe, with a PK, FK based on join_df, and all unique synonyms
    """
    # Convert wide meta_df to long table
    # Drop 'variable' col (leave only unique ID and synonyms), drop duplicates
    synonyms = pd.melt(meta_df, id_vars=[unique_id])[
        [unique_id, 'value']].drop_duplicates()

    # Drop all rows where value is NA
    synonyms = synonyms[synonyms['value'].notna()]

    # Join with join_df based on unique_id
    synonyms = pd.merge(synonyms, join_df, left_on=unique_id,
                        right_on='name', how='left')[['id', 'value']]

    return synonyms