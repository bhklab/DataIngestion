import glob
import os
import re
import numpy as np
import pandas as pd
import swifter  # Library to parallelize apply statements automagically

pset_name = 'GDSC_v1'
file_path = os.path.join('data', 'rawdata')
slot_names = ['curation', 'drug', 'molecularProfiles',
              'sensitivity', 'annotation', 'cell']


def read_pset(pset_name, file_path, slot_names=['curation', 'drug', 'molecularProfiles', 'sensitivity', 'annotation', 'cell']):
    """
    Read in all the data associated with a PharmacoSet from the .csv files exported by the writeToCsv method from rPharmacoDI.

    @param pset_name: [`string`] Name of the PharmacoSet object as it appears in the directory name for files exported using
        rPharmacoDIs writeToCsv function.
    @param file_path: [`string`] Path to the directory where PharmacoSet data is stored.
    @param slot_names: [`list`] A list of PharmacoSet slots to read in. Defaults to all slots in a PharmacoSet.

    @return: [`DataFrame`] A DataFrame with the columns slot, containing the name of the slot data is from, one or more subitem columns, indicating
        the subitem in a slot a data came from, file path, the path data was read from and data, which contains the object associated with the
        selected slot and subitem.
    """
    # Use regex to find the appropriate files
    pset_dir = glob.glob(f'{os.path.join(file_path, pset_name)}_PSet')[0]
    if pset_dir is None:
        raise ValueError(
            f'No PSet directory named {pset_name} could be found in {file_path}')

    # List all files for the select PSet, then split on $ to make a DataFrame
    pset_files = pd.Series(os.listdir(pset_dir))
    pset_files_df = pset_files.str.split('$', expand=True)

    # Build the file paths to read in data for each row of the DataFrame
    pset_files_df['file_paths'] = [os.path.join(
        pset_dir, file_name) for file_name in pset_files]

    # Rename columns
    pset_files_df.columns = [
        'slot', *[f'subitems{i}' for i in range(1, pset_files_df.shape[1] - 1)], 'file_paths']

    # Read in PSet data
    pset_files_df['data'] = pset_files_df['file_paths'].swifter.apply(
        read_pset_file)

    # Drop file_paths column by reference
    pset_files_df.drop('file_paths', axis='columns', inplace=True)

    # Process id columns to use the proper slot names
    pset_files_df.iloc[:, 0:-1] = pset_files_df.iloc[:, 0:-
                                                     1].apply(lambda col: col.str.replace('.*@|.csv.gz$|.txt', ''))

    return pset_files_df


# ---- Helper methods for read_csv
def read_pset_file(file_path):
    """Deal with text files which can't be read in with pd.read_csv"""
    if '.csv.gz' in file_path:
        return pd.read_csv(file_path)
    elif '.txt' in file_path:
        with open(file_path) as f:
            text = [line for line in f]
        return text
    else:
        raise ValueError(
            f'Unsupported file type passed to this function from: {file_path}')


def pset_df_to_nested_dict(df):
    """
    Recrusively turn unique values in the first column of a DataFrame into dictionary keys until only 1 column remains.

    @param df: [`DataFrame`] With one or more ID columns before a data column.
    @return: [`dict`] A nested dict with levels equal to the number of columns of the original DataFrame minus one
    """
    for key_check in pd.unique(df[df.columns[0]]):
        if key_check is not None:
            # Recursively nest the first column as key and remaning columns as values
            return {key: pset_df_to_nested_dict(df.loc[df[df.columns[0]] == key, df.columns[1:]]) if
                    df[1:].shape[1] > 2 else df.loc[df[df.columns[0]] == key, 'data'].values[0] for
                    key in pd.unique(df[df.columns[0]])}
        else:
            # Return the data column if there no no key
            return df.iloc[:, -1].values[0]


def build_pset_tables(pset_dict, pset_name):
    """
    Build all tables for this....
    """

    datasets_df = pd.DataFrame({'id': pset_name, 'name': pset_name})

    tissues_df, drugs_df, genes_df = build_primary_tables(
        pset_dict)

    cells_df = build_cells_df(pset_dict, tissues_df)
    targets_df = build_targets_df(pset_dict, genes_df)
    drug_annotations_df, gene_annotations_df = build_annotation_dfs(
        pset_dict, drugs_df, genes_df)

    datasets_cells_df = build_datasets_cells_df(
        pset_dict, cells_df, datasets_df)
    mol_cels_df = build_mol_cells_df(pset_dict, datasets_cells_df)

    clinical_trials_df = build_clinical_trials_df(pset_dict, drugs_df)
    experiments_df = build_experiments_df(
        pset_dict, cells_df, drugs_df, datasets_df, tissues_df)
    drug_targets_df = build_drug_targets_df(pset_dict, drugs_df, targets_df)

    dose_responses_df = build_dose_response_df(pset_dict, experiments_df)
    profiles_df = build_profiles_df(pset_dict, experiments_df)
    dataset_statistics_df = build_dataset_stats_df(pset_dict, datasets_df)

    gene_drugs_df = build_gene_drugs_table()

    # TODO - how to not repeat the same code over and over? how to iterate through all tables
    tissues_df.to_csv(f'{pset_name}_Tissues.csv.gz', header=True,
                      index=False, compression='infer')


# --- PRIMARY TABLES --------------------------------------------------------------------------

def build_primary_tables(pset_dict):
    """
    Build the tissues, drug, gene DataFrames.

    @param pset_dict: [`dict`]
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
    """
    genes = pd.Series(pd.unique(
        pset_dict['molecularProfiles']['rna']['rowData']['EnsemblGeneId']))

    return pd.DataFrame({'id': genes, 'name': genes})


# TODO - check that you should get it from 'cell' rather than 'curation'
def tissue_table(pset_dict):
    """
    Build the tissue dataframe.
    """
    tissues = pd.Series(pd.unique(pset_dict['cell']['tissueid']))
    return pd.DataFrame({'id': tissues, 'name': tissues})


# Make drugs df; TODO - confirm whether to use drugid (?) or DRUG_NAME (targets)
def drug_table(pset_dict):
    """
    Build the drug dataframe.
    """
    drugs = pd.Series(pd.unique(pset_dict['drug']['drugid']))
    return pd.DataFrame({'id': drugs, 'name': drugs})

# --- OTHER TABLES -----------------------------------------------------------------------------

def build_targets_df(pset_dict, genes_df):
    # Make the targets df -- NEED TO ADD GENE ID TODO - how to relate target to genes??
    targets = pd.Series(pd.unique(pset_dict['drug']['TARGET']))
    targets_df = pd.DataFrame({'id': targets, 'name': targets})
    targets_df['gene_id'] = np.nan #TODO - fill in
    # target_df = pd.merge(pset_dict['drug'][['TARGET']], gene_df)[
    #     ['TARGET', 'id']]
    # target_df.columns = ['name', 'gene_id']

    return targets_df


# TODO - confirm that you're using the correct cell id
def build_cells_df(pset_dict, tissues_df):
    """
    ~ description ~ (build out cell and target dataframes)

    @param pset_dict: [`dict`]
    @param tissue_df: [`DataFrame`]
    @return: [(`DataFrame`, `DataFrame`)]
    """
    cells_df = pset_dict['cell'][['cellid', 'tissueid']]
    cells_df.columns = ['name', 'tissue_id']
    cells_df['id'] = cells_df.loc[:, ('name')] # TODO - check copy warning

    return cells_df[['id', 'name', 'tissue_id']]


def build_annotation_dfs(pset_dict):
    """
    build out drug_annotations and gene_annotations dfs

    @param pset_dict: [`dict`]
    @return: [(`DataFrame`, `DataFrame`)]
    """
    # Make gene_annotations df
    gene_annotations_df = pset_dict['molecularProfiles']['rna']['rowData'][['EnsemblGeneId', 'Symbol']]
    gene_annotations_df.columns = ['gene_id', 'symbol']
    gene_annotations_df.loc[:, ('gene_seq_start')] = np.nan  #TODO - fill in
    gene_annotations_df.loc[:, ('gene_seq_end')] = np.nan    #TODO - fill in & fix copy warning

    # Make drug_annotations df
    drug_annotations_df = pset_dict['drug'][['rownames', 'smiles', 'inchikey', 'cid', 'FDA']]
    drug_annotations_df.columns = ['drug_id', 'smiles', 'inchikey', 'pubchem', 'fda_status']

    return drug_annotations_df, gene_annotations_df


def build_datasets_cells_df(pset_dict, cells_df, dataset_id):
    """
    add documentation
    """
    datasets_cells_df = pd.DataFrame({'dataset_id': dataset_id, 'cell_id': cells_df['id']})
    datasets_cells_df['id'] = datasets_cells_df.index + 1

    return datasets_cells_df[['id', 'dataset_id', 'cell_id']]


def build_mol_cells_df(pset_dict, datasets_cells_df):
    # Get the number of times each cellid appears in colData
    num_profiles = pset_dict['molecularProfiles']['rna']['colData']['cellid'].value_counts()

    # Join with datasets cells on cellid
    mol_cells_df = pd.merge(datasets_cells_df, num_profiles, left_on='cell_id', right_on=num_profiles.index, how='left')
    mol_cells_df.rename(columns={'cellid': 'num_prof'}, inplace=True)

    # Replace any NaN in the num_profiles column with 0
    mask = mol_cells_df.query('num_prof.isna()').index
    mol_cells_df.loc[mask, 'num_prof'] = 0

    # Set mDataType TODO - find mDataType???
    mol_cells_df['mDataType'] = np.nan

    return mol_cells_df[['id', 'cell_id', 'dataset_id', 'mDataType', 'num_prof']]


def build_clinical_trials_df(pset_dict, drugs_df):
    return None

# ---- Build dose response table
# Chris' implementation
# TODO:: Do Python functions pass my reference or copy by default?
# If pass by reference may need to make a copy before using inplace=TRUE argument
# to prevent modifying the original experiments table
# NOTE: database tables should always use singular names
def build_dose_response_df(pset_dict, experiment_df):
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


def build_drug_targets_df(pset_dict, drug_df, target_df):
    """
    add documentation
    """
    drug_targets_df = pset_dict['drug'][['drugid', 'TARGET']]
    drug_targets_df.columns = ['drug_id', 'target_id']
    drug_targets_df['id'] = drug_targets_df.index + 1

    return drug_targets_df[['id', 'drug_id', 'target_id']]


def build_experiment_df(pset_dict, cell_df, drug_df, dataset_id):
    # Extract relelvant experiment columns
    experiments_df = pset_dict['sensitivity']['info'][[
        'exp_id', 'cellid', 'drugid']]

    # Join with cell_df
    experiments_df = pd.merge(experiments_df, cell_df, left_on='cellid', right_on='name',
                              how='left')[['exp_id', 'id', 'drugid', 'tissue_id']]
    # Rename cell_id column (FK)
    experiments_df.rename(columns={"id": "cell_id"}, inplace=True)

    # Join with drug_df
    experiments_df = pd.merge(experiments_df, drug_df, left_on='drugid', right_on='name',
                              how='left')[['exp_id', 'cell_id', 'id', 'tissue_id']]
    # Rename drug_id column (FK)
    experiments_df.rename(columns={"id": "drug_id"}, inplace=True)

    # Add dataset_id
    experiments_df['dataset_id'] = dataset_id

    # TODO - get rid of exp_id; check if col order matters

    return experiments_df


def build_profiles_df(pset_dict):


# --- GENE_DRUGS TABLE --------------------------------------------------------------------------

gene_sig_file_path = os.path.join(
    file_path, 'gene_signatures', 'pearson_perm_res')


def read_gene_sig(pset_name, file_path):
    # Find correct pset gene signature CSV file
    pset_file = glob.glob(f'{os.path.join(file_path, pset_name)}.csv')[0]
    if pset_file is None:
        raise ValueError(
            f'No PSet gene signatures file named {pset_name} could be found in {file_path}')

    # Read csv file and return df
    return pd.read_csv(pset_file)


def build_gene_drugs_df(gene_sig_df, genes_df, drugs_df, tissues_df, dataset_id):
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

    # Add all foreign keys (gene_id, drug_id, tissue_id, dataset_id)

    # Join with genes_df
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
    gene_drugs_df.rename(columns={"id": "tissue_id"}, inplace=True)

    # Add dataset id
    gene_drugs_df['dataset_id'] = dataset_id

    # Reorder columns (omit red cols and PK)
    return gene_drugs_df[['gene_id', 'drug_id', 'estimate', 'se', 'n', 'pvalue',
                          'df', 'fdr', 'dataset_id', 'sens_stat', 'tissue_id',
                          'mDataType', 'tested_in_human_trials', 'in_clinical_trials']]


# --- END GENE_DRUGS TABLE --------------------------------------------------------------------------

if __name__ == "__main__":
    # pset_names = ['GDSC_v1', 'GDSC_v2', 'gCSI',
    #               'FIMM', 'CTRPv2', 'CCLE', 'GRAY', 'UHNBreast']

    pset_names = ['GDSC_v1']

    # Convert pset files to nested dictionaries
    pset_dicts = []
    for pset in pset_names:
        pset_df = read_pset(pset_name, file_path)
        pset_dict = pset_df_to_nested_dict(pset_df)
        pset_dicts.append(pset_dict)


# Where is clinical_trials info stored?
   # Chris:
    # There is a ruby script to fetch this information
    # Code is at: https://github.com/bhklab/PharmacoDB-web/blob/master/lib/tasks/update_clinical_trials.rake
    # I know nothing about Ruby so I hope there is documentation :)
# oncotrees, cellosaurus, dataset_statistics, profiles, mol_cells
    # Chris:
    # oncotress - we don't have the data yet
    # dataset_statistics -
    # profiles - will be found in sensitivityProfiles, may not have all statistics yet
