import os
import glob
import pandas as pd
import numpy as np
import dask.dataframe as dd

pset_name = 'GDSC_v1'
# the directory where you want to store the tables (processed data)
file_path = 'procdata'
gene_sig_file_path = os.path.join(
    'data', 'rawdata', 'gene_signatures', 'pearson_perm_res')  # the directory with data for gene_drugs

# the maximum number of rows in a pandas df that can be written to a csv in a reasonable amount of time
dask_threshold = 1000000


def build_pset_tables(pset_dict, pset_name, file_path):
    """
    Build all tables for this

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    pset_dfs = {}

    # Build primary tables
    pset_dfs['dataset'] = pd.DataFrame({'id': pset_name, 'name': pset_name})
    tissues_df, drugs_df, genes_df = build_primary_tables(
        pset_dict)
    pset_dfs.update({'tissue': tissues_df, 'drug': drugs_df, 'gene': genes_df})

    # Build annotation tables
    gene_annotations_df, drug_annotations_df = build_annotation_dfs(pset_dict)
    pset_dfs.update({'gene_annotations': gene_annotations_df,
                     'drug_annotations': drug_annotations_df})

    # Build secondary tables
    pset_dfs['cell'] = build_cells_df(pset_dict, tissues_df)
    pset_dfs['target'] = build_targets_df(pset_dict, genes_df)
    pset_dfs['experiments'] = build_experiments_df(
        pset_dict, pset_dfs['cell'], pset_name)
    pset_dfs['clinical_trials'] = build_clinical_trials_df(pset_dict, drugs_df)
    pset_dfs['drug_targets'] = build_drug_targets_df(pset_dict)
    pset_dfs['dose_responses'] = build_dose_response_df(
        pset_dict, pset_dfs['experiments'])
    pset_dfs['profiles'] = build_profiles_df(pset_dict)

    # Build gene drugs table
    pset_dfs['gene_drugs'] = build_gene_drugs_df(gene_sig_file_path, pset_name)

    # Build summary/stats tables
    pset_dfs['datasets_cells'] = build_datasets_cells_df(
        pset_dict, pset_dfs['cell'], pset_name)
    pset_dfs['mol_cells'] = build_mol_cells_df(
        pset_dict, pset_dfs['datasets_cells'], pset_dfs['gene_drugs'])
    pset_dfs['dataset_statistics'] = build_dataset_stats_df(
        pset_dict, pset_dfs, pset_name)

    # Write all tables to csv
    write_dfs_to_csv(pset_dfs, pset_name, file_path)


def write_dfs_to_csv(pset_dfs, pset_name, df_dir):
    """
    @param pset_dfs: [`dict`] A dictionary of all the pset DataFrames you want to write to csv
    @param pset_name: [`string`] The name of the pset, used to create a subdirectory for the tables in this pset
    @param df_dir: [`string`] The name of the directory to hold all the tables
    @return [`None`]
    """
    file_path = os.path.join(df_dir, pset_name)

    # TODO - test out and get a better dask_threshold if needed
    for df in pset_dfs.keys():
        if len(df.index) < dask_threshold:
            # Use pandas to convert df to csv
            df.to_csv(os.path.join(file_path, df,
                                   f'{pset_name}_{df}.csv'), index=False)
        else:
            # Convert pandas df into dask df TODO - check how many partitions it makes and adjust if necessary
            dask_df = dd.from_pandas(pset_dfs[df])
            # Write dask_df to csv
            dd.to_csv(dask_df, os.path.join(
                file_path, df, f'{pset_name}_{df}-*.csv'))


# --- PRIMARY TABLES --------------------------------------------------------------------------

def build_primary_tables(pset_dict):
    """
    Build the tissues, drug, gene DataFrames.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    @return: [(`DataFrame`, `DataFrame`, `DataFrame`)]
    """
    tissues_df = tissue_table(pset_dict)
    drugs_df = drug_table(pset_dict)
    genes_df = gene_table(pset_dict)

    return tissues_df, drugs_df, genes_df


# TODO - check that you're getting the correct ID
def gene_table(pset_dict):
    """
    Build the gene dataframe.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    genes = pd.Series(pd.unique(
        pset_dict['molecularProfiles']['rna']['rowData']['EnsemblGeneId']))

    return pd.DataFrame({'id': genes, 'name': genes})


# TODO - check that you should get it from 'cell' rather than 'curation'
def tissue_table(pset_dict):
    """
    Build the tissue dataframe.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    tissues = pd.Series(pd.unique(pset_dict['cell']['tissueid']))
    return pd.DataFrame({'id': tissues, 'name': tissues})


# Make drugs df; TODO - confirm whether to use drugid (?) or DRUG_NAME (targets)
def drug_table(pset_dict):
    """
    Build the drug dataframe.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    drugs = pd.Series(pd.unique(pset_dict['drug']['drugid']))
    return pd.DataFrame({'id': drugs, 'name': drugs})


# --- ANNOTATION TABLES -----------------------------------------------------------------------

def build_annotation_dfs(pset_dict):
    """
    build out drug_annotations and gene_annotations dfs

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    @return: [(`DataFrame`, `DataFrame`)]
    """
    gene_annotations_df = build_gene_annotations_df(pset_dict)
    drug_annotations_df = build_drug_annotations_df(pset_dict)

    return gene_annotations_df, drug_annotations_df


def build_gene_annotations_df(pset_dict):
    """
    add documentation

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    # Make gene_annotations df
    gene_annotations_df = pset_dict['molecularProfiles']['rna']['rowData'][[
        'EnsemblGeneId', 'Symbol']]
    gene_annotations_df.columns = ['gene_id', 'symbol']
    gene_annotations_df.loc[:, ('gene_seq_start')] = np.nan  # TODO - fill in
    # TODO - fill in & fix copy warning
    gene_annotations_df.loc[:, ('gene_seq_end')] = np.nan

    return gene_annotations_df


def build_drug_annotations_df(pset_dict):
    """
    add documentation

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    # Make drug_annotations df
    drug_annotations_df = pset_dict['drug'][[
        'rownames', 'smiles', 'inchikey', 'cid', 'FDA']]
    drug_annotations_df.columns = [
        'drug_id', 'smiles', 'inchikey', 'pubchem', 'fda_status']

    return drug_annotations_df


# --- OTHER TABLES ----------------------------------------------------------------------------

def build_targets_df(pset_dict, genes_df):
    """
    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    # Make the targets df -- NEED TO ADD GENE ID TODO - how to relate target to genes??
    targets = pd.Series(pd.unique(pset_dict['drug']['TARGET']))
    targets_df = pd.DataFrame({'id': targets, 'name': targets})
    targets_df['gene_id'] = np.nan  # TODO - fill in
    # target_df = pd.merge(pset_dict['drug'][['TARGET']], gene_df)[
    #     ['TARGET', 'id']]
    # target_df.columns = ['name', 'gene_id']

    return targets_df


# TODO - confirm that you're using the correct cell id
def build_cells_df(pset_dict, tissues_df):
    """
    ~ description ~ (build out cell and target dataframes)

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    @param tissue_df: [`DataFrame`]
    @return: [(`DataFrame`, `DataFrame`)]
    """
    cells_df = pset_dict['cell'][['cellid', 'tissueid']]
    cells_df.columns = ['name', 'tissue_id']
    cells_df['id'] = cells_df.loc[:, ('name')]  # TODO - check copy warning

    return cells_df[['id', 'name', 'tissue_id']]


def build_clinical_trials_df(pset_dict, drugs_df):
    """
    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    return pd.DataFrame()

# ---- Build dose response table
# Chris' implementation
# TODO:: Do Python functions pass my reference or copy by default?
# If pass by reference may need to make a copy before using inplace=TRUE argument
# to prevent modifying the original experiments table
# NOTE: database tables should always use singular names
def build_dose_response_df(pset_dict, experiment_df):
    """
    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    # Get dose and response info from pset
    dose = pset_dict['sensitivity']['raw.Dose']
    response = pset_dict['sensitivity']['raw.Viability']

    # rename columns so they can be coerced to int later
    rename_dict = {f'dose.{n}': str(n) for n in np.arange(1, dose.shape[0])}
    dose.rename(columns=rename_dict, inplace=True)
    response.rename(columns=rename_dict, inplace=True)

    # reshape the DataFrames using melt and pivot to go from 'wide' to 'long' or back, respectively
    # these are much faster than using Python loops because all of the looping is done in the Pandas  C++ code
    dose = dose.melt(id_vars='.exp_id', value_name='dose',
                     var_name='dose_id').dropna()
    dose['dose_id'] = dose.dose_id.astype('int')
    response = response.melt(
        id_vars='.exp_id', value_name='response', var_name='dose_id').dropna()
    response['dose_id'] = response.dose_id.astype('int')

    # set indexes for faster joins (~3x)
    dose.set_index(['.exp_id', 'dose_id'], inplace=True)
    response.set_index(['.exp_id', 'dose_id'], inplace=True)

    # because I set indexes on both table I don't need to specify on
    dose_response_df = pd.merge(
        dose, response, left_index=True, right_index=True).reset_index()

    dose_response_df.rename(columns={'.exp_id': 'exp_id'}, inplace=True)
    dose_response_df.set_index('exp_id', inplace=True)
    experiment_df.set_index('exp_id', inplace=True)

    dose_response_df = pd.merge(dose_response_df, experiment_df[['id']],
                                right_index=True, left_index=True).reset_index()
    dose_response_df.drop('exp_id', 1, inplace=True)
    dose_response_df.rename(columns={'id': 'experiment_id'}, inplace=True)
    dose_response_df.index = dose_response_df.index + 1

    return dose_response_df


def build_drug_targets_df(pset_dict):
    """
    add documentation

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    drug_targets_df = pset_dict['drug'][['drugid', 'TARGET']]
    drug_targets_df.columns = ['drug_id', 'target_id']
    drug_targets_df['id'] = drug_targets_df.index + 1

    return drug_targets_df[['id', 'drug_id', 'target_id']]


def build_experiments_df(pset_dict, cells_df, dataset_id):
    """
    add description

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    # Extract relelvant experiment columns
    experiments_df = pset_dict['sensitivity']['info'][[
        'exp_id', 'cellid', 'drugid']]

    # Rename columns
    experiments_df.columns = ['id', 'cell_id', 'drug_id']

    # Add datset_id column TODO - copy warning
    experiments_df['dataset_id'] = dataset_id

    # Add tissue_id column by joining with cells_df
    experiments_df = pd.merge(
        experiments_df, cells_df[['name', 'tissue_id']], left_on='cell_id', right_on='name',
        how='left')[['id', 'cell_id', 'drug_id', 'dataset_id', 'tissue_id']]

    return experiments_df


def build_profiles_df(pset_dict):
    """
    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    """
    # Get profiles info
    profiles_df = pset_dict['sensitivity']['profiles']

    # Get experiment_id by joining with sensitivity info df
    profiles_df = pd.merge(profiles_df, pset_dict['sensitivity']['info'][[
                           '.rownames', 'exp_id']], on='.rownames', how='left')

    # Drop .rownames column after join
    profiles_df.drop('.rownames', axis='columns', inplace=True)

    # Rename columns TODO - double check tomorrow that these are correctly renamed
    profiles_df.columns = ['AAC', 'IC50',
                           'HS', 'Einf', 'EC50', 'experiment_id']

    # Add DSS columns - TODO get these values from somewhere?
    profiles_df['DSS1'] = np.nan
    profiles_df['DSS2'] = np.nan
    profiles_df['DSS3'] = np.nan

    return profiles_df[['experiment_id', 'HS', 'Einf', 'EC50', 'AAC', 'IC50', 'DSS1', 'DSS2', 'DSS3']]


# --- GENE_DRUGS TABLE --------------------------------------------------------------------------

def read_gene_sig(pset_name, file_path):
    # Find correct pset gene signature CSV file
    pset_file = glob.glob(f'{os.path.join(file_path, pset_name)}.csv')[0]
    if pset_file is None:
        raise ValueError(
            f'No PSet gene signatures file named {pset_name} could be found in {file_path}')

    # Read csv file and return df
    return pd.read_csv(pset_file)


def build_gene_drugs_df(gene_sig_file_path, pset_name):
    """
    add documentation
    """
    # Get gene_sig_df from gene_sig_file
    gene_sig_df = read_gene_sig(pset_name, gene_sig_file_path)

    # Extract relevant columns
    gene_drugs_df = gene_sig_df[[
        'gene_id', 'drug', 'estimate', 'n', 'pvalue', 'df', 'fdr', 'tissue', 'mDataType']]
    # TODO - cannot find 'se', 'sens_stat' -- is one of these 'significant'???
    # Chris: You will determine significance based on the fdr (false discovery rate) at alpha = 0.05, it will be TRUE or FALSE (or 1 or 0)

    # Chris: 'se' - leave NA/Null/None for now, it will be added as a column to the gene signatures the next time we run them.
    gene_drugs_df['se'] = np.nan
    # Chris: 'sens_stat' - I will add this to the function for extracting per PSet gene signatures - for now it is always 'AAC' (Area above dose-response curve)
    gene_drugs_df['sens_stat'] = 'AAC'
    # TODO - cannot find 'drug_like_molecule', 'in_clinical_trials'
    # Chris: Have renamed it to tested_in_human_trials, it will indicate a 1 if it has ever been tested in a human clinical trial (even if it failed)
    # Chris: Source for this data will be clinicaltrails.gov
    # TODO - check out API, leave NA for now
    gene_drugs_df['tested_in_human_trials'] = np.nan
    gene_drugs_df['in_clinical_trials'] = np.nan

    # Rename foreign key columns
    gene_drugs_df.rename(columns={'drug': 'drug_id', 'tissue': 'tissue_id'})

    # Don't need joins anymore because just need names not indices??
    """ # Join with genes_df
    gene_drugs_df = pd.merge(genes_df, gene_drugs_df,
                             left_on='name', right_on='gene_id', how='right')
    # Drop gene name columns after merge & rename 'id' column
    gene_drugs_df.drop(['name', 'gene_id'], axis='columns', inplace=True)
    gene_drugs_df.rename(columns={"id": "gene_id"}, inplace=True)

    # Join with drugs_df
    gene_drugs_df = pd.merge(drugs_df, gene_drugs_df,
                             left_on='name', right_on='drug', how='right')
    # TODO - this merges better on drugid than DRUG_NAME
    # Drop drug name columns after merge & rename 'id' column
    gene_drugs_df.drop(['name', 'drug'], axis='columns', inplace=True)
    gene_drugs_df.rename(columns={"id": "drug_id"}, inplace=True)

    # Join with tissues_df
    gene_drugs_df = pd.merge(gene_drugs_df, tissues_df,
                             left_on='tissue', right_on='name', how='left')
    # Drop tissue name columns after merge & rename 'id' column
    gene_drugs_df.drop(['name', 'tissue'], axis='columns', inplace=True)
    gene_drugs_df.rename(columns={"id": "tissue_id"}, inplace=True) """

    # Add dataset id
    gene_drugs_df['dataset_id'] = pset_name

    # Add missing columns (TODO - get this data)
    gene_drugs_df['tstat'] = np.nan
    gene_drugs_df['fstat'] = np.nan
    gene_drugs_df['FWER_genes'] = np.nan
    gene_drugs_df['FWER_drugs'] = np.nan
    gene_drugs_df['FWER_all'] = np.nan
    gene_drugs_df['BF_p_all'] = np.nan
    gene_drugs_df['meta_res'] = np.nan

    # Reorder columns
    return gene_drugs_df[['gene_id', 'drug_id', 'estimate', 'se', 'n', 'tstat', 'fstat',
                          'pvalue', 'df', 'fdr', 'FWER_genes', 'FWER_drugs', 'FWER_all',
                          'BF_p_all', 'meta_res', 'dataset_id', 'sens_stat', 'tissue_id',
                          'mDataType', 'tested_in_human_trials', 'in_clinical_trials']]


# --- STATS/SUMMARY TABLES -----------------------------------------------------------------------

def build_datasets_cells_df(pset_dict, cells_df, dataset_id):
    """
    Builds a join table summarizing the cell lines that are in this dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    @param cells_df: [`DataFrame`] The cell table for this dataset
    @param dataset_id: [`string`] The name/id of the datset
    @return: [`DataFrame`] The join table with all cell lines in this dataset
    """
    datasets_cells_df = pd.DataFrame(
        {'dataset_id': dataset_id, 'cell_id': cells_df['id']})
    datasets_cells_df['id'] = datasets_cells_df.index + 1

    return datasets_cells_df[['id', 'dataset_id', 'cell_id']]


def build_mol_cells_df(pset_dict, datasets_cells_df, gene_drugs_df):
    """
    Builds a table that summarizes which molecular data types have calculated stats in
    this dataset for each cell line, and the number of profiles of each molecular data type.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    @param datasets_cells_df: [`DataFrame`] A table containing all the cells in this  
        dataset and the datset name/id.
    @return: [`DataFrame`] 
    """
    mol_cells_df = pd.DataFrame(
        columns=['id', 'cell_id', 'dataset_id', 'mDataType', 'num_prof'])
    molecularTypes = pd.unique(gene_drugs_df['mDataType'])

    for mDataType in molecularTypes:
        # Get the number of times each cellid appears in colData for that mDataType
        num_profiles = pset_dict['molecularProfiles'][mDataType]['colData']['cellid'].value_counts(
        )

        # Join with datasets cells on cellid
        df = pd.merge(datasets_cells_df, num_profiles,
                      left_on='cell_id', right_on=num_profiles.index, how='left')

        # Rename num_profiles column
        df.rename(columns={'cellid': 'num_prof'}, inplace=True)

        # Set mDataType column to the current molecular type
        df['mDataType'] = mDataType

        # Append to mol_cells_df
        mol_cells_df = mol_cells_df.append(df)

    # Replace any NaN in the num_profiles column with 0
    mask = mol_cells_df.query('num_prof.isna()').index
    mol_cells_df.loc[mask, 'num_prof'] = 0

    return mol_cells_df


def build_dataset_stats_df(pset_dict, pset_dfs, pset_name):
    """
    Summarizes how many cell lines, tissues, drugs, and experiments are contained
    within the dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    @param pset_name: [`string`]
    @return: [(`DataFrame`, `DataFrame`)]
    """
    return pd.DataFrame({
        'id': 1,
        'dataset_id': pset_name,
        'cell_lines': len(pset_dfs['cell'].index),
        'tissues': len(pset_dfs['tissue'].index),
        'drugs': len(pset_dfs['drug'].index),
        'experiments': len(pset_dfs['experiments'].index)
    })
