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



# Table-making functions

def make_pset_dfs(pset_dict):

    # primary tables

    # make tissues df
    tissues = pd.Series(pd.unique(pset_dict['cell']['tissueid']))
    tissues_df = pd.DataFrame({'id': tissues.index, 'name': tissues})

    # make drugs df -- should i be using drugid or DRUG_NAME?? -- looks like DRUG_NAME has better mapping
    drugs = pd.Series(pd.unique(pset_dict['drug']['DRUG_NAME']))
    drugs_df = pd.DataFrame({'id': drugs.index, 'name': drugs})

    datasets = None #name of pset
    genes = None #molecularProfiles - RNA - ensembleGeneId

    #oncotrees
    #cell_synonyms
    #drug_synonyms
    #tissue_synonyms

    # non-primary tables

    # make the cells df
    cells_df = pd.merge(pset_dict['cell'], tissues_df, left_on='tissueid', right_on='name', how='left')[['cellid', 'name']]
    cells_df.columns = ['name', 'tissue_id']
    cells_df['id'] = cells_df.index

    # make the drug_annotations df
    drug_annotations_df = pset_dict['drug'][['rownames', 'smiles', 'inchikey', 'cid', 'FDA']]
    drug_annotations_df.columns = ['name', 'smiles', 'inchikey', 'pubchem', 'fda_status']
    drug_annotations_df = pd.merge(drugs_df, drug_annotations_df, on='name', how='right')
    # drop name column once you've merged on it
    drug_annotations_df.drop('name', axis='columns', inplace=True)

    # make the targets df -- NEED TO ADD GENE ID
    targets = pd.Series(pd.unique(pset_dict['drug']['TARGET']))
    targets_df = pd.DataFrame({'id': targets.index, 'name': targets})

    # make the drug_targets df
    # TODO - NEEDS TO BE MODIFIED TO JOIN WITH TARGET TABLE
    drug_targets_df = pset_dict['drug'][['DRUG_NAME', 'TARGET']]
    drug_targets_df = pd.merge(drugs_df, drug_targets_df, left_on='name', right_on='DRUG_NAME', how='right')
    # PROBLEM - some drug names don't map?? Ask Chris
    # drop name columns after merge
    drug_targets_df.drop(['name', 'DRUG_NAME'], axis='columns', inplace=True)
    # rename columns and add primary key
    drug_targets_df.columns = ['drug_id', 'target_id']
    drug_targets_df['id'] = drug_targets_df.index

    # make experiments df
    experiments_df = pset_dict['sensitivity']['info'][['exp_id', 'cellid', 'drugid']]
    experiments_df = pd.merge(experiments_df, cells_df[['name', 'tissue_id']], left_on='cellid', right_on='name', how='left')
    # drop name column after merge
    experiments_df.drop('name', axis='columns', inplace=True)
    # rename columns
    experiments_df.columns = ['exp_id', 'cell_id', 'drug_id', 'tissue_id']
    # TODO - need to add dataset_id, and eventually get rid of exp_id
    # add primary key
    experiments_df['id'] = experiments_df.index

    # make the dose_responses df
    #dose_responses = #pset_dict['sensitivity']['raw.Dose']
    #dose_responses = #pset_dict['sensitivity']['raw.Viability']
    pd.merge(experiments_df[['id', 'exp_id']], pset_dict['sensitivity']['raw.Dose'], left_on='exp_id', right_on='.exp_id', how='right')
    doses_df = pset_dict['sensitivity']['raw.Dose']
    responses_df = pset_dict['sensitivity']['raw.Viability']

    for i in doses_df.index:
        experiment = doses_df['.exp_id'][i]
        #for i in range(1, doses_df.shape[1]):

    #[f'dose.{i}' for i in range(1, doses_df.shape[1])]

    pd.DataFrame({'exp_id': doses_df['.exp_id'][0],
                  'dose': doses_df.iloc[0][1:]})

                  #doses_df[[f'dose.{i}' for i in range(1, doses_df.shape[1])]][0]

    #check joins with .notna() on cols
