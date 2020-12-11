import os
import re
import glob
import numpy as np
import pandas as pd
import dask.dataframe as dd
from join_pset_tables import *

read_file_path = os.path.join('data', 'procdata')
write_file_path = os.path.join('data', 'latest')
metadata_path = os.path.join('data', 'metadata')
annotations_path = os.path.join(metadata_path, 'Annotations')


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

drug_syns_file = "drugs_with_ids.csv"
cell_syns_file = "cell_annotation_all.csv"


# TODO - 3 rows; 954 rows
def cell_synonyms(cell_df, cell_syns_file, annotations_path):
    # Get metadata file
    cell_names_df = get_metadata(cell_syns_file, annotations_path)

    # Find all columns relevant to cellid
    pattern = re.compile('cellid')
    cell_columns = cell_names_df[[
        col for col in cell_names_df.columns if pattern.search(col)]]

    # Get all unique synonyms and join with cells_df
    cell_synonyms = melt_and_join(cell_columns, 'unique.cellid', cell_df)

    # Rename columns
    cell_synonyms.columns = ['cell_id', 'cell_name']

    # Add primary key and blank col for dataset_id (TODO)
    cell_synonyms['id'] = cell_synonyms.index + 1
    cell_synonyms['dataset_id'] = np.nan

    # Reorder columns to match ERD
    return cell_synonyms[['id', 'cell_id', 'dataset_id', 'cell_name']]


def tissue_synonyms(tissues_df, tissue_syns_file, annotations_path):
    # Get metadata file
    tissues_metadata = get_metadata(tissue_syns_file, annotations_path)

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


def drug_synonyms(drugs_df, drug_syns_file, annotations_path):
    # Get metadata file
    drugs_metadata = get_metadata(drug_syns_file, annotations_path)

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


# --- TARGET TABLES --------------------------------------------------------------------------

drugbank_file = "drugbank_targets_has_ref_has_uniprot.csv"
uniprot_ensemble_file = "uniprot_ensembl_mappings.csv"


def build_metadata_dfs(join_dfs, uniprot_ensemble_file, metdata_path):
    gene_mappings_df = get_metadata(uniprot_ensemble_file, metadata_path)
    target_file_path = os.path.join(metadata_path, 'target')
    target_df = build_target_df(target_file_path, gene_mappings_df, join_dfs['gene'])
    #drug_target_df = build_drug_target_df(target_df, join_dfs['drug'])
        #TODO: generalize drug-target mappings

    # TODO: how to pass this many parameters
    cell_syn_df = cell_synonyms(join_dfs['cell'], cell_syns_file, metadata_path)
    tissue_syn_df = tissue_synonyms(join_dfs['tissue'], cell_syns_file, metadata_path)
    drug_syn_df = drug_synonyms(join_dfs['drug'], drug_syns_file, metadata_path)


#TODO - generalize so it goes through all targets (ChEMBL, etc.)
def build_target_df(target_file_path, gene_mappings_df, gene_df):
    """
    Map all Drugbank (TODO - generalize) targets to EnsemblGene IDs
    based on UniProt IDs (TODO - join with gene table?).
    
    @param drugbank_targets_df: [`DataFrame`] A table of all Drugbank targets
    @param gene_mappings_df: [`DataFrame`] A table mapping UniProt IDs to EnsemlGene IDs
    @return: [`DataFrame`] A table of target to gene mappings.
    """
    # Get all target files (Drugbank, ChemBL, etc.)
    target_files = glob.glob(target_file_path)
    if len(target_files) == 0:
        raise ValueError(
            f'No target directory named {target_file_path} could be found!')

    # Load all target tables
    dfs = [] 
    for target_file in os.listdir(target_files):
        df = (get_metadata(target_file, target_file_path)) #change so that it uses dask to load
        dfs.append(df[['name', 'polypeptide.external.identifiers.UniProtKB']])

    # Append all target tables to each other
    target_df = concat_tables(dfs)
    
    # Rename columns for easier merging
    target_df.rename(columns={'polypeptide.external.identifiers.UniProtKB': 'UniProtId'}, inplace=True)
    
    # Join on UniProt ID to get EnsemblGene IDs (1300 rows)
    target_df = pd.merge(target_df, gene_mappings_df, on='UniProtId', how='inner')
    target_df.drop('UniProtId', axis='columns', inplace=True)
    target_df.rename(columns={'EnsemblGeneId': 'gene_id'}, inplace=True)

    # Join with gene table
    target_df = safe_merge(target_df, gene_df, 'gene_id')
    
    return target_df


def build_drug_target_df(drugbank_target_df, target_df, drug_df):
    # Get all drug-target mappings from Drugbank
    drug_target_df = drugbank_target_df[['name', 'drugName']].drop_duplicates()
    drug_target_df.rename(columns={'name': 'target_id', 'drugName': 'drug_id'}, inplace=True)

    # Build foreign keys by joining with drug table and target table
    drug_target_df = safe_merge(drug_target_df, drug_df, 'drug_id')
    drug_target_df = safe_merge(drug_target_df, target_df, 'target_id')

    # Re-index
    drug_target_df = reindex_table(drug_target_df)
    
    return drug_target_df
