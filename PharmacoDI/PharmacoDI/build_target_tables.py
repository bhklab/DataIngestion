import os
import requests
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import datatable as dt
from PharmacoDI.combine_pset_tables import index_and_write, rename_and_key, join_tables
from PharmacoDI.get_chembl_drug_targets import parallelize


drugbank_file = os.path.join(
    "data", "metadata", "drugbank_targets_has_ref_has_uniprot.csv")
chembl_file = os.path.join('data', 'metadata', 'chembl_drug_targets.csv')
output_dir = os.path.join("data", "demo")

# TODO: need error handling for missing directories
# TODO: need to abstract some functions


def build_target_tables(drugbank_file, chembl_file, output_dir):
    """
    Build the target and drug target tables using data from Drugbank
    and ChEMBL.

    @param drugbank_file: [`string`] The full file path to Drugbank targets
    @param chembl_file: [`string`] The full file path to ChEMBL targets
    @param output_dir: [`string`] The directory of all final PharmacoDB tables
    @return: None
    """
    # Get Drugbank data
    drugbank_df = pd.read_csv(drugbank_file)
    drugbank_df.rename(columns={'polypeptide.external.identifiers.UniProtKB': 'uniprot_id',
                                'drugName': 'drug_name'}, inplace=True)

    # Get ChEMBL data
    chembl_df = pd.read_csv(chembl_file, index_col=0)
    chembl_df.rename(columns={'pref_name': 'name',
                              'accession': 'uniprot_id'}, inplace=True)

    target_df = build_target_table(chembl_df, drugbank_df, output_dir)
    build_drug_target_table(chembl_df, drugbank_df, target_df, output_dir)


def build_target_table(chembl_df, drugbank_df, output_dir):
    """
    Using data from the Drugbank and ChEMBL drug target files and
    the UniProt API, build the target table.

    @param chembl_df: [`pd.DataFrame`] The ChEMBL drug target table
    @param drugbank_df: [`pd.DataFrame`] The DrugBank drug target table
    @param output_dir: [`string`] The file path to write the final target table
    @return: [`pd.DataFrame`] The target table
    """
    # Combine ChEMBL and Drugbank tables to make target table
    target_df = pd.concat([chembl_df[['name', 'uniprot_id']].copy(),
                           drugbank_df[['name', 'uniprot_id']].copy()])
    target_df.drop_duplicates(inplace=True)
    uniprot_ids = pd.Series(pd.unique(target_df['uniprot_id']))

    # Retrieve Uniprot-ENSEMBL gene ID mappings
    uniprot_ensembl_mappings = pd.concat(
        parallelize(uniprot_ids, map_uniprot_to_ensembl, 1000))
    uniprot_ensembl_mappings.drop_duplicates(inplace=True)

    # Join target table with gene table based on uniprot-ensembl mappings
    target_df = pd.merge(target_df, uniprot_ensembl_mappings,
                         on='uniprot_id', how='left')
    # TODO: 3 uniprot ids don't join to a gene id (O90777, P0A953, P0ABF6)

    # Load gene table from output_dir
    gene_file = os.path.join(output_dir, 'gene.csv')
    if not os.path.exists(gene_file):
        raise FileNotFoundError(f"ERROR: There is no gene file in {output_dir}!")
    gene_df = pd.read_csv(gene_file)
    gene_df.rename(columns={'name': 'gene_id'}, inplace=True)

    # Join target table with gene table based on gene name
    target_df = pd.merge(target_df, gene_df, on='gene_id', how='left')
    # TODO: 99 rows (excluding gene_id=NA) don't map to any genes in gene table
    target_df.drop(columns=['uniprot_id', 'gene_id'], inplace=True)
    target_df.rename(columns={'id': 'gene_id'}, inplace=True)
    # idk how there are duplicates but there are
    target_df.drop_duplicates(inplace=True)

    target_df = index_and_write(dt.Frame(target_df), 'target', output_dir)
    return target_df # convert to pandas


def build_drug_target_table(chembl_df, drugbank_df, target_df, output_dir):
    """
    Using data from the Drugbank and ChEMBL drug target files and 
    the target table, build the drug target table.

    @param chembl_df: [`pd.DataFrame`] The ChEMBL drug target table
    @param drugbank_df: [`pd.DataFrame`] The DrugBank drug target table
    @param target_df: [`pd.DataFrame`] The target table
    @param output_dir: [`string`] The file path with all final PharmacoDB tables
    @return: [`pd.DataFrame`] The drug target table
    """
    # Load drug synonym table from output_dir
    drug_synonym_file = os.path.join(output_dir, 'drug_synonym.csv')
    if not os.path.exists(drug_synonym_file):
        raise FileNotFoundError(f"ERROR: There is no drug synonym file in {output_dir}!")
    drug_syn_df = pd.read_csv(drug_synonym_file, dtype={'drug_id': 'int32'})
    
    # Join drugbank df with drug table (TODO: are we really using drug name to map?)
    drugbank_df = pd.merge(drugbank_df, drug_syn_df, on='drug_name')
    # TODO: from 7521 down to only 122 rows :/

    # Combine ChEMBL and Drugbank tables to make drug target table
    drug_target_df = pd.concat([chembl_df[['name', 'drug_id']].copy(),
                                drugbank_df[['name', 'drug_id']].copy()])
    drug_target_df.drop_duplicates(inplace=True)

    # Join with target table
    target_df = rename_and_key(target_df, 'target_id')
    drug_target_df = join_tables(dt.Frame(drug_target_df), target_df, 'target_id')
    # Drop rows with no target_id, drop duplicates
    drug_target_df = drug_target_df[dt.f.target_id >= 1, ['drug_id', 'target_id']]
    drug_target_df = drug_target_df[0, :, dt.by(drug_target_df.names)]

    drug_target_df = index_and_write(drug_target_df, 'drug_target', output_dir)
    return drug_target_df


def map_uniprot_to_ensembl(uniprot_ids):
    """
    Use the UniProt API to retrieve the ENSEMBL gene IDs
    corresponding to the UniProt IDS.

    @param uniprot_ids: [`list(string)`] A list of UniProt IDs.
    @return: [`pd.DataFrame`] A table mapping UniProt IDs to ENSEMBL gene IDs.
    """
    # Make API call
    params = {
        'from': 'ID',
        'to': 'ENSEMBL_ID',
        'format': 'tab',
        'query': " ".join(uniprot_ids)
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
    gene_id_df.rename(columns={0: 'uniprot_id', 1: 'gene_id'}, inplace=True)

    return gene_id_df

