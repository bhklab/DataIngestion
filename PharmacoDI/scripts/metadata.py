import os
import re
import glob
import numpy as np
import pandas as pd

metadata_path = os.path.join('data', 'metadata', 'Annotations')


def get_annotations(file_name, metadata_path):
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
tissue_syns_file = "cell_annotation_all.csv"


# TODO - 3 rows; 954 rows
def cell_synonyms(cells_df, cell_syns_file, metadata_path):
    # Get metadata file
    cell_names_df = get_annotations(cell_syns_file, metadata_path)

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


def tissue_synonyms(tissues_df, tissue_syns_file, metadata_path):
    # Get metadata file
    tissues_metadata = get_annotations(tissue_syns_file, metadata_path)

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


def drug_synonyms(drugs_df, drug_syns_file, metadata_path):
    # Get metadata file
    drugs_metadata = get_annotations(drug_syns_file, metadata_path)

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

#TODO - generalize so it goes through all targets (ChEMBL, etc.)
# TODO - clean this up, too many parameters?
def build_target_dfs(drugbank_file, uniprot_ensemble_file, metadata_path, metadata_dfs):
    # Get metadata
    drugbank_targets_df = get_annotations(drugbank_file, metadata_path)
    gene_mappings_df = get_annotations(uniprot_ensemble_file, metadata_path)

    # Build targets table
    #TODO - get gene table
    metadata_dfs['target'] = build_target_df(drugbank_targets_df, gene_mappings_df)

    # Build drug_targets table
    #TODO - get drug table
    metadata_dfs['drug_targets'] = build_drug_targets_df(drugbank_targets_df, metadata_dfs['target'], None)


#TODO - generalize so it goes through all targets (ChEMBL, etc.)
def build_target_df(drugbank_targets_df, gene_mappings_df):
    """
    Map all Drugbank (TODO - generalize) targets to EnsemblGene IDs
    based on UniProt IDs (TODO - join with gene table?).
    
    @param drugbank_targets_df: [`DataFrame`] A table of all Drugbank targets
    @param gene_mappings_df: [`DataFrame`] A table mapping UniProt IDs to EnsemlGene IDs
    @return: [`DataFrame`] A table of target to gene mappings.
    """
    
    # Get all targets & UniProtIds from Drugbank
    target_df = drugbank_targets_df[['name', 'polypeptide.external.identifiers.UniProtKB']].drop_duplicates()
    
    # Rename columns for easier merging
    target_df.rename(columns={'polypeptide.external.identifiers.UniProtKB': 'UniProtId'}, inplace=True)
    
    # Join on UniProt ID to get EnsemblGene IDs
    target_df = pd.merge(target_df, gene_mappings_df, on='UniProtId', how='inner')
    
    # Drop UniProtId and rename gene column
    target_df.drop('UniProtId', axis='columns', inplace=True)
    target_df.rename(columns={'EnsemblGeneId': 'gene_id'}, inplace=True)
    
    # Re-index to start at 1
    target_df.index = target_df.index + 1
    
    return target_df


#TODO - generalize, and make cleaner(?)
def build_drug_targets_df(drugbank_targets_df, target_df, drug_df):

    # Get all drug-target mappings from Drugbank
    drug_targets_df = drugbank_targets_df[['name', 'drugName']].drop_duplicates()

    # Target df with index and target name only
    # find a better way to do this? or just make sure everything has an 'id' column?
    targets = target_df[['name']]
    targets['id'] = targets.index

    # Join with target_df
    drug_targets_df = pd.merge(drug_targets_df, targets, on='name', how='inner')

    # Drop target name column and rename columns
    drug_targets_df.drop('name', axis='columns', inplace=True)
    drug_targets_df.rename(columns={'drugName': 'name', 'id': 'target_id'}, inplace=True)

    # Join with drug_df
    drug_targets_df = pd.merge(drug_targets_df, drug_df, on='name', how='inner')

    # Drop drug name column and rename columns
    drug_targets_df.drop('name', axis='columns', inplace=True)
    drug_targets_df.rename(columns={'id': 'drug_id'}, inplace=True)

    # Reset index
    drug_targets_df.reset_index(drop=True, inplace=True)
    drug_targets_df.index = drug_targets_df.index + 1

    return drug_targets_df