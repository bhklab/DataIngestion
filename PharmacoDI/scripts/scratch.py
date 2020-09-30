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



gene_sig_file_path = os.path.join(file_path, 'gene_signatures', 'pearson_perm_res')

# Read gene signature
def read_gene_sig(pset_name, file_path):
    # Find correct pset gene signature CSV file 
    pset_file = glob.glob(f'{os.path.join(file_path, pset_name)}.csv')[0]
    if pset_file is None:
        raise ValueError(
            f'No PSet gene signatures file named {pset_name} could be found in {file_path}')

    # Read csv file and return df
    return pd.read_csv(pset_file)



annotations_file_path = os.path.join('data', 'metadata', 'Annotations')

# get metadata files
def get_annotations(file_name, annotations_path):
    # Find correct metadata annotations CSV file
    annotations_file = glob.glob({os.path.join(annotations_path, file_name)})[0]
    if annotations_file is None:
        raise ValueError(
            f'No metadata file named {file_name} could be found in {annotations_path}')
    
    # Read csv file and return df
    return pd.read_csv(annotations_file)


# Make cell synonyms df
def cell_synonyms(cells_df, file_name, file_path):
    #
    file_name = 'cell_annotation_all.csv'
    cell_names_df = get_annotations(file_name, file_path)

    #
    pattern = re.compile('cellid')
    cell_columns = cell_names_df[[col for col in cell_names_df.columns if pattern.search(col)]]

    # no
    cell_syn_df = pd.DataFrame(columns=['cell_id', 'cell_name'])

    for row in cell_columns.iterrows():
        (index, cell_row) = row
        cell_id = cells_df[cells_df['name'] == cell_names_df['X'][index]]['id']
        synonyms = pd.Series(pd.unique(cell_row)).dropna()
        cell_syn_df = cell_syn_df.append(
            pd.DataFrame({ 'cell_id': cell_id, 'cell_name': synonyms}), 
            ignore_index=True)

    # using X
    # do I ever use curation?
    # do I need datasetid? it's in red..
    #lots unmatched -- fix

# Make gene drugs df
def gene_drugs_df(gene_sig_df, genes_df, drugs_df, tissues_df):
    # Extract relevant columns
    gene_drugs_df = gene_sig_df[['gene_id', 'drug', 'estimate', 'n', 'pvalue', 'df', 'fdr', 'tissue', 'mDataType']]
    #TODO - cannot find 'se', 'sens_stat' -- is one of these 'significant'???
    #TODO - cannot find 'level' ('lower'? 'upper'?)
    #TODO - cannot find 'drug_like_molecule', 'in_clinical_trials'
    
    # Join with genes_df
    gene_drugs_df = pd.merge(genes_df, gene_drugs_df, left_on='name', right_on='gene_id', how='right')
    # Drop gene name columns after merge & rename 'id' column
    gene_drugs_df.drop(['name', 'gene_id'], axis='columns', inplace=True)
    gene_drugs_df = gene_drugs_df.rename(columns={"id": "gene_id"})

    # Join with drugs_df
    gene_drugs_df = pd.merge(drugs_df, gene_drugs_df, left_on='name', right_on='drug', how='right')
        # TODO - this merges better on drugid than DRUG_NAME
    # Drop drug name columns after merge & rename 'id' column
    gene_drugs_df.drop(['name', 'drug'], axis='columns', inplace=True)
    gene_drugs_df = gene_drugs_df.rename(columns={"id": "drug_id"})

    # Join with tissues_df
    gene_drugs_df = pd.merge(gene_drugs_df, tissues_df, left_on='tissue', right_on='name', how='left')
        # TODO - test this merge
    # Drop tissue name columns after merge & rename 'id' column
    gene_drugs_df.drop(['name', 'tissue'], axis='columns', inplace=True)
    gene_drugs_df = gene_drugs_df.rename(columns={"id": "tissue_id"})

    # dataset id..

    return gene_drugs_df

# TODO - Make dose resonses df
def dose_responses_df(pset_dict):
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
                  'dose': doses_df.iloc[0][1:],
                  'response': responses_df.iloc[0][1:]})

                  #doses_df[[f'dose.{i}' for i in range(1, doses_df.shape[1])]][0]

#TODO - datasets df

def make_pset_dfs(pset_dict, dataset_id):

    # make tissues df
    tissues = pd.Series(pd.unique(pset_dict['cell']['tissueid']))
    tissues_df = pd.DataFrame({'id': tissues.index, 'name': tissues})

    # make drugs df -- should i be using drugid or DRUG_NAME?? -- looks like DRUG_NAME has better mapping
    drugs = pd.Series(pd.unique(pset_dict['drug']['DRUG_NAME']))
    drugs_df = pd.DataFrame({'id': drugs.index, 'name': drugs})

    # makes genes df
    genes = pd.Series(pd.unique(pset_dict['molecularProfiles']['rna']['rowData']['EnsemblGeneId'])) # check this
    genes_df = pd.DataFrame({'id': genes.index, 'name': genes})

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


    #check joins with .notna() and .isna() on cols

    dose_responses_df = dose_responses_df(pset_dict)

    #pset_name??? gene_sig_file_path?
    gene_sig_df = read_gene_sig(pset_name, gene_sig_file_path)
    gene_drugs_df = gene_drugs_df(gene_sig_df, genes_df, drugs_df, tissues_df)
