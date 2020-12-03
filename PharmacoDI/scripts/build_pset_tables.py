import os
import glob
import pandas as pd
import numpy as np
import dask.dataframe as dd

pset_name = 'GDSC_v1'
# the directory where you want to store the tables (processed data)
file_path = 'procdata'
metadata_path = os.path.join("data", "metadata", "Annotations")
gene_sig_file_path = os.path.join(
    'data', 'rawdata', 'gene_signatures')  # the directory with data for gene_drugs

# the maximum number of rows in a pandas df that can be written to a csv in a reasonable amount of time
dask_threshold = 1000000


def build_pset_tables(pset_dict, pset_name, file_path):
    """
    Build all tables for a dataset and write them to a directory of all processed data.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the pset
    @param pset_name: [`string`] The name of the dataset
    @param file_path: [`string`] The file path to the directory containing processed data
    @return: [`None`]
    """
    pset_dfs = {}

    # Build primary tables
    print('Building primary tables...')
    # pd.DataFrame({'id': pset_name, 'name': pset_name})
    pset_dfs['dataset'] = pd.Series(pset_name, name='name')
    build_primary_tables(pset_dict, pset_dfs)

    # Build annotation tables
    print('Building annotation tables...')
    build_annotation_dfs(pset_dict, pset_dfs)

    # Build secondary tables
    print('Building secondary tables...')
    pset_dfs['cell'] = build_cells_df(pset_dict, pset_dfs['tissue'])
    pset_dfs['experiments'] = build_experiments_df(
        pset_dict, pset_dfs['cell'], pset_name)
    # pset_dfs['clinical_trials'] = build_clinical_trials_df(pset_dict, pset_dfs['drug']) # TODO - API
    pset_dfs['dose_responses'] = build_dose_response_df(
        pset_dict, pset_dfs['experiments'])
    pset_dfs['profiles'] = build_profiles_df(pset_dict)

    # Build gene drugs table
    if os.path.exists(os.path.join(gene_sig_file_path, pset_name)):
        print('Building gene drugs table...')
        pset_dfs['gene_drugs'] = build_gene_drugs_df(
            gene_sig_file_path, pset_name)

    # Build summary/stats tables
    print('Building summary/stats tables...')
    pset_dfs['datasets_cells'] = build_datasets_cells_df(
        pset_dict, pset_dfs['cell'], pset_name)
    if 'gene_drugs' in pset_dfs:
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
    # Make sure directory for this PSet exists
    if not os.path.exists(file_path):
        os.mkdir(file_path)

    for df_name in pset_dfs.keys():
        print(f'Writing {df_name} table to csv...')
        df = pset_dfs[df_name]
        
        df.index.rename('id', inplace=True)

        df_path = os.path.join(file_path, df_name)
        # Make sure subdirectory for this table exists
        if not os.path.exists(df_path):
            os.mkdir(df_path)
        else:
            # Clear all old files so there aren't overlaps/errors
            for f in os.listdir(df_path):
                os.remove(os.path.join(df_path, f))

        if len(df.index) < dask_threshold:
            # Use pandas to convert df to csv
            df.to_csv(os.path.join(df_path,
                f'{pset_name}_{df_name}.csv'), index=True)
        else:
            # Convert pandas df into dask df TODO - adjust chunksize if needed
            dask_df = dd.from_pandas(df, chunksize=500000)
            # Write dask_df to csv
            dd.to_csv(dask_df, os.path.join(
                df_path, f'{pset_name}_{df_name}-*.csv'), index=True)


# --- PRIMARY TABLES --------------------------------------------------------------------------

def build_primary_tables(pset_dict, pset_dfs):
    """
    Build the tissue, drug, and gene tables for a dataset and add
    them to the dictionary containing all tables for that dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_dfs: [`dict`] A dictionary of all the tables corresponding to a PSet
    @return: [`None`]
    """
    pset_dfs['tissue'] = build_tissue_table(pset_dict)
    pset_dfs['drug'] = build_drug_table(pset_dict)
    # Don't make gene table if there are no molecular profiles (TODO - check with chris)
    if 'molecularProfiles' in pset_dict:
        pset_dfs['gene'] = build_gene_table(pset_dict)


def build_gene_table(pset_dict):
    """
    Build a table containing all genes in a dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] The gene table
    """
    genes = pd.Series([], name='name')
    for mDataType in pset_dict['molecularProfiles']:
        if 'EnsemblGeneId' in pset_dict['molecularProfiles'][mDataType]['rowData']:
            genes = genes.append(pd.Series(pd.unique(
                pset_dict['molecularProfiles'][mDataType]['rowData']['EnsemblGeneId']), name='name'))

    genes.drop_duplicates(inplace=True)
    return genes  # pd.DataFrame({'id': genes, 'name': genes})


def build_tissue_table(pset_dict):
    """
    Build a table containing all tissues in a dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] The tissue table
    """
    tissues = pd.Series(pd.unique(pset_dict['cell']['tissueid']), name='name')
    return tissues  # pd.DataFrame({'id': tissues, 'name': tissues})


# Make drugs df; TODO - confirm whether to use drugid (?) or DRUG_NAME (targets)
def build_drug_table(pset_dict):
    """
    Build a table containing all drugs in a dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] The drug table
    """
    drugs = pd.Series(pd.unique(pset_dict['drug']['drugid']), name='name')
    return drugs  # pd.DataFrame({'id': drugs, 'name': drugs})


# --- ANNOTATION TABLES -----------------------------------------------------------------------

def build_annotation_dfs(pset_dict, pset_dfs):
    """
    Build the drug_annotations and gene_annotations tables for a dataset and add
    them to the dictionary containing all tables for that dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_dfs: [`dict`] A dictionary of all the tables corresponding to a PSet
    @return: [`None`]
    """
    # Only make gene_annotations table if gene table already exists
    if 'gene' in pset_dfs:
        pset_dfs['gene_annotations'] = build_gene_annotations_df(pset_dict)

    pset_dfs['drug_annotations'] = build_drug_annotations_df(pset_dict)


def build_gene_annotations_df(pset_dict):
    """
    Build a table mapping each gene in a dataset to its gene annotations.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] A table of all gene annotations, mapped to genes
    """
    gene_annotations_df = pd.DataFrame(
        columns=['gene_id', 'symbol', 'gene_seq_start', 'gene_seq_end'])

    for mDataType in pset_dict['molecularProfiles']:
        df = pset_dict['molecularProfiles'][mDataType]['rowData'].copy()
        # Get gene annotation columns
        cols = ['.features']
        if 'Symbol' in df.columns:
            cols.append('Symbol')

        df = df[cols]
        df.rename(columns={'.features': 'gene_id', 'Symbol': 'symbol'}, inplace=True)
        gene_annotations_df = gene_annotations_df.append(df)

    gene_annotations_df.drop_duplicates(subset=['gene_id'], inplace=True)

    return gene_annotations_df


def build_drug_annotations_df(pset_dict):
    """
    Build a table mapping each drug in a dataset to its drug annotations.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] A table of all drug annotations, mapped to drugs
    """
    # Make drug_annotations df
    drug_annotations_df = pset_dict['drug'][[
        'rownames', 'smiles', 'inchikey', 'cid', 'FDA']].copy()
    drug_annotations_df.rename(
        columns={'rownames': 'drug_id', 'cid': 'pubchem', 'FDA': 'fda_status'}, inplace=True)

    return drug_annotations_df


# --- OTHER TABLES ----------------------------------------------------------------------------

def load_drugbank_mappings(gene_file, file_path):
    """
    Read all Drugbank gene ID mappings from the directory file_path.

    @param gene_file: [`string`] The name of the Drugbank file
    @param file_path: [`string`] The directory that holds all gene mapping files
    @return: [`DataFrame`] A dataframe containing all Drugbank gene mappings
    """
    # Find correct CSV file
    drugbank_file = glob.glob(os.path.join(file_path, gene_file))
    if len(drugbank_file) == 0:
        raise ValueError(
            f'No Drugbank file named {gene_file} could be found in {file_path}')

    # Read csv file and return df
    return pd.read_csv(drugbank_file[0])


# TODO - confirm that you're using the correct cell id
def build_cells_df(pset_dict, tissues_df):
    """
    Build a table containing all the cells in a dataset, mapped to their tissues.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param tissues_df: [`DataFrame`] A table of all tissues in the dataset
    @return: [`DataFrame`] A table of all cell lines, mapped to tissues
    """
    cells_df = pset_dict['cell'][['cellid', 'tissueid']].copy()
    cells_df.rename(columns={'cellid': 'name',
                             'tissueid': 'tissue_id'}, inplace=True)
    # cells_df['id'] = cells_df.loc[:, 'name']  # TODO - check copy warning
    # maybe try cells_df.loc[:, 'id'] = cells_df[:, 'name']

    return cells_df  # [['id', 'name', 'tissue_id']]

# TODO


def build_clinical_trials_df(pset_dict, drugs_df):
    """
    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
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
    Build a table that, for each experiment in a dataset, lists the drug that was
    tested, the doses in which that drug was administered, and the viability responses 
    corresponding to all the doses.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] A table with all the dose-response mappings for each experiment
    """
    # Get dose and response info from pset
    dose = pset_dict['sensitivity']['raw.Dose']
    response = pset_dict['sensitivity']['raw.Viability']

    # rename columns so they can be coerced to int later
    dose_pattern = dose.columns[1][:-1]
    rename_dict = {f'{dose_pattern}{n}': str(
        n) for n in np.arange(1, dose.shape[0])}

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

    dose_response_df.rename(columns={'.exp_id': 'experiment_id'}, inplace=True)

    # Not necessary since merging will occur after we join all PSets
    #dose_response_df.set_index('exp_id', inplace=True)
    #experiment_df.set_index('exp_id', inplace=True)

    # dose_response_df = pd.merge(dose_response_df, experiment_df[['id']],
    #                            right_index=True, left_index=True).reset_index()
    #dose_response_df.drop('exp_id', 1, inplace=True)
    #dose_response_df.rename(columns={'id': 'experiment_id'}, inplace=True)
    #dose_response_df.index = dose_response_df.index + 1

    return dose_response_df


def build_experiments_df(pset_dict, cells_df, dataset_id):
    """
    Build a table with all the experiments in a dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param cells_df: [`DataFrame`] A table of all the cells in the PSet and their tissues
    @param dataset_id: [`string`] The name of the PSet 
    @return: [`DataFrame`] A table containing all experiments in the dataset
    """
    # Extract relelvant experiment columns
    experiments_df = pset_dict['sensitivity']['info'][[
        '.rownames', 'cellid', 'drugid']].copy()

    # Rename columns
    experiments_df.rename(
        columns={'.rownames': 'experiment_id', 'cellid': 'cell_id', 'drugid': 'drug_id'}, inplace=True)

    # Add datset_id column TODO - copy warning
    experiments_df['dataset_id'] = dataset_id

    # Add tissue_id column by joining with cells_df
    experiments_df = pd.merge(
        experiments_df, cells_df[['name', 'tissue_id']], left_on='cell_id', right_on='name',
        how='left')[['experiment_id', 'cell_id', 'drug_id', 'dataset_id', 'tissue_id']]

    experiments_df.rename(columns={'experiment_id': 'name'}, inplace=True)

    return experiments_df


def build_profiles_df(pset_dict):
    """
    TODO: ask Chris

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] A table containing all statistics for each profile in the PSet (?)
    """
    # Get profiles info
    if 'E_inf' in pset_dict['sensitivity']['profiles'].columns:
        profiles_df = pset_dict['sensitivity']['profiles'][[
            '.rownames', 'aac_recomputed', 'ic50_recomputed', 'HS', 'E_inf', 'EC50']].copy()
        profiles_df.rename(columns={'.rownames': 'experiment_id', 'aac_recomputed': 'AAC',
                                    'ic50_recomputed': 'IC50', 'E_inf': 'Einf'}, inplace=True)
    else:
        profiles_df = pset_dict['sensitivity']['profiles'][[
            '.rownames', 'aac_recomputed', 'ic50_recomputed', 'slope_recomputed', 'einf', 'ec50']].copy()
        profiles_df.rename(columns={'.rownames': 'experiment_id', 'aac_recomputed': 'AAC', 'slope_recomputed': 'HS',
                                    'ic50_recomputed': 'IC50', 'einf': 'Einf', 'ec50': 'EC50'}, inplace=True)

    # Add DSS columns - TODO get these values when they are eventually computed
    profiles_df['DSS1'] = np.nan
    profiles_df['DSS2'] = np.nan
    profiles_df['DSS3'] = np.nan

    return profiles_df[['experiment_id', 'HS', 'Einf', 'EC50', 'AAC', 'IC50', 'DSS1', 'DSS2', 'DSS3']]


# --- GENE_DRUGS TABLE --------------------------------------------------------------------------

def read_gene_sig(pset_name, file_path):
    """
    Read all gene signatures for a PSet (to be used in gene_drugs table) from the directory file_path.

    @param pset_name: [`string`] The name of the PSet
    @param file_path: [`string`] The directory that holds all gene signature files
    @return: [`DataFrame`] A dataframe containing all gene signatures for the PSet
    """
    # Find correct pset gene signature CSV file
    pset_file = glob.glob(
        f'{os.path.join(file_path, pset_name, pset_name)}_gene_sig.csv')
    if len(pset_file) == 0:
        raise ValueError(
            f'No PSet gene signatures file named {pset_name} could be found in {file_path}')

    # Read csv file and return df
    return pd.read_csv(pset_file[0])


def build_gene_drugs_df(gene_sig_file_path, pset_name):
    """
    TODO - ask Chris to explain this table again

    @param gene_sig_file_path: [`string`] The file path to the directory containing the gene
        signatures for each PSet
    @param pset_name: [`string`] The name of the PSet
    @return: [`DataFrame`] The gene_drugs table for this PSet, containing all stats (?)
    """
    # Get gene_sig_df from gene_sig_file
    gene_sig_df = read_gene_sig(pset_name, gene_sig_file_path)

    # Extract relevant columns
    gene_drugs_df = gene_sig_df[[
        'gene', 'drug', 'estimate', 'n', 'pvalue', 'df', 'fdr', 'tissue', 'mDataType']].copy()
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
    gene_drugs_df.rename(
        columns={'gene': 'gene_id', 'drug': 'drug_id', 'tissue': 'tissue_id'}, inplace=True)

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

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param cells_df: [`DataFrame`] The cell table for this dataset
    @param dataset_id: [`string`] The name/id of the datset
    @return: [`DataFrame`] The join table with all cell lines in this dataset
    """
    datasets_cells_df = pd.DataFrame(
        {'dataset_id': dataset_id, 'cell_id': cells_df['name']})
    datasets_cells_df['id'] = datasets_cells_df.index + 1

    # datasets_cells_df[['id', 'dataset_id', 'cell_id']]
    return datasets_cells_df[['dataset_id', 'cell_id']]


def build_mol_cells_df(pset_dict, datasets_cells_df, gene_drugs_df):
    """
    Builds a table that summarizes the number of profiles, per cell line, per molecular data
    type, in this dataset. (Only considers molecular data types for which there are sens stats?)

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param datasets_cells_df: [`DataFrame`] A table containing all the cells in this  
        dataset and the datset name/id.
    @param gene_drugs_df: [`DataFrame`] A table containing (..?)
    @return: [`DataFrame`] The table with the number of profiles for each cell line, for each
        molecular data type
    """
    mol_cells_df = pd.DataFrame(
        columns=['id', 'cell_id', 'dataset_id', 'mDataType', 'num_prof'])
    molecularTypes = pd.unique(gene_drugs_df['mDataType'])

    for mDataType in molecularTypes:
        if 'molecularProfiles' in pset_dict:
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
        else:
            # If PSet contains no molecular profiles, set num_prof to 0
            # for all celll lines and all molecular data types
            df = datasets_cells_df.copy()
            df['mDataType'] = mDataType
            df['num_prof'] = 0

        # Append to mol_cells_df
        mol_cells_df = mol_cells_df.append(df)

    # Replace any NaN in the num_profiles column with 0
    mask = mol_cells_df.query('num_prof.isna()').index
    mol_cells_df.loc[mask, 'num_prof'] = 0
    mol_cells_df['num_prof'].astype('int32', copy=False)

    return mol_cells_df


def build_dataset_stats_df(pset_dict, pset_dfs, pset_name):
    """
    Summarizes how many cell lines, tissues, drugs, and experiments are contained
    within the dataset.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_dfs: [`dict`] A dictionary of all the data frames associated with the PSet
    @param pset_name: [`string`] The name of the PSet
    @return: [`DataFrame`] A one-row table with the summary stats for this PSet
    """
    return pd.DataFrame({
        'dataset_id': [pset_name],
        'cell_lines': [len(pset_dfs['cell'].index)],
        'tissues': [len(pset_dfs['tissue'].index)],
        'drugs': [len(pset_dfs['drug'].index)],
        'experiments': [len(pset_dfs['experiments'].index)]
    })
