import os
import pandas as pd
import numpy as np
from PharmacoDI.get_chembl_drug_targets import parallelize
from multiprocessing import Pool, cpu_count
import requests

drugbank_file = os.path.join(
    "data", "metadata", "drugbank_targets_has_ref_has_uniprot.csv")
chembl_file = os.path.join('data', 'metadata', 'chembl_drug_targets.csv')
gene_dir = os.path.join("data", "latest", "gene")
drug_synonym_dir = os.path.join("data", "latest", "drug_synonym")
target_file = os.path.join("data", "latest", "target", "target.csv")
drug_target_file = os.path.join(
    "data", "latest", "drug_target", "drug_target.csv")

# TODO: need error handling for missing directories
# TODO: need to update docstrings
# TODO: need to abstract some functions


def bulid_target_tables(drugbank_file, chembl_file, gene_dir, drug_synonym_dir,
                        target_file, drug_target_file):
    """
    Build the target and drug target tables using data from Drugbank
    and ChEMBL.

    :param: drugbank_file The full file path to Drugbank targets
    :param: chembl_file The full file path to ChEMBL targets
    :param: gene_dir The directory with the gene table
    :param: drug_synonym_dir The directory with the drug synonyms table
    """
    # Get Drugbank data
    drugbank_df = pd.read_csv(drugbank_file)
    drugbank_df.rename(
        columns={'polypeptide.external.identifiers.UniProtKB': 'uniprot_id', 'drugName': 'drug_name'}, inplace=True)

    # Get ChEMBL data
    chembl_df = pd.read_csv(chembl_file, index_col=0)
    chembl_df.rename(columns={'pref_name': 'name',
                              'accession': 'uniprot_id'}, inplace=True)

    target_df = build_target_table(
        chembl_df, drugbank_df, gene_dir, target_file)
    build_drug_target_table(chembl_df, drugbank_df, target_df,
                            drug_synonym_dir, drug_target_file)


def build_target_table(chembl_df, drugbank_df, gene_dir, target_file):
    """
    Using data from the Drugbank and ChEMBL drug target files and
    the UniProt API, build the target table.
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
    gene_df = load_gene_table(gene_dir)
    target_df = pd.merge(target_df, gene_df, on='gene_id', how='left')
    # TODO: 99 rows (excluding gene_id=NA) don't map to any genes in gene table
    target_df.drop(columns=['uniprot_id', 'gene_id'], inplace=True)
    target_df.rename(columns={'id': 'gene_id'}, inplace=True)

    # idk how there are duplicates but there are
    target_df.drop_duplicates(inplace=True)
    target_df.reset_index(drop=True)
    target_df.index += 1
    target_df.index = target_df.index.rename('id')
    target_df.to_csv(target_file)

    return target_df


def build_drug_target_table(chembl_df, drugbank_df, target_df,
                            drug_synonym_dir, drug_target_file):
    """
    Using data from the Drugbank and ChEMBL drug target files and 
    the target table, build the drug target table.
    """
    # Join drugbank df with drug table (TODO: are we really using drug name to map?)
    drug_syn_df = load_drug_synonym_table(drug_synonym_dir)
    drugbank_df = pd.merge(drugbank_df, drug_syn_df, on='drug_name')
    # TODO: from 7521 down to only 122 rows :/
    # retaining NA drug_ids only improves the join a bit, up to 158 rows

    # Combine ChEMBL and Drugbank tables to make drug target table
    drug_target_df = pd.concat([chembl_df[['name', 'drug_id']].copy(),
                                drugbank_df[['name', 'drug_id']].copy()])
    drug_target_df.drop_duplicates(inplace=True)

    # Join with target table
    target_df['id'] = target_df.index
    drug_target_df = pd.merge(drug_target_df, target_df, on='name')
    drug_target_df = drug_target_df[['drug_id', 'id']].copy().drop_duplicates()
    drug_target_df.rename(columns={'id': 'target_id'}, inplace=True)

    drug_target_df.reset_index(drop=True)
    drug_target_df.index += 1
    drug_target_df.index = drug_target_df.index.rename('id')
    drug_target_df.to_csv(drug_target_file)

    return drug_target_df


def map_uniprot_to_ensembl(uniprot_ids):
    """
    Use the UniProt API to retrieve the ENSEMBL gene IDs
    correspondong to the UniProt IDS.

    @param uniprot_ids: [`list(string)`] A list of UniProt IDs.
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


# These will be abstracted
def load_gene_table(gene_dir):
    if not os.path.exists(gene_dir):
        print(f"ERROR: the gene directory {gene_dir} doesn't exist!")
    gene_df = pd.DataFrame()
    for f in os.listdir(gene_dir):
        df = pd.read_csv(os.path.join(gene_dir, f),
                         dtype={'id': pd.Int64Dtype()})
        gene_df = gene_df.append(df)
    gene_df.rename(columns={'name': 'gene_id'}, inplace=True)
    return gene_df


def load_drug_synonym_table(drug_synonym_dir):
    if not os.path.exists(drug_synonym_dir):
        print(
            f"ERROR: the drug synonym directory {drug_synonym_dir} doesn't exist!")
    drug_syn_df = pd.DataFrame()
    for f in os.listdir(drug_synonym_dir):
        df = pd.read_csv(os.path.join(drug_synonym_dir, f), index_col='id')
        drug_syn_df = drug_syn_df.append(df)
    # TODO: temporary; there shouldn't be NAs in final synonym table
    drug_syn_df = drug_syn_df.query('drug_id.notna()').copy()
    drug_syn_df = drug_syn_df.astype({'drug_id': 'int32'})
    return drug_syn_df
