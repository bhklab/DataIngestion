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


gene_sig_file_path = os.path.join(
    file_path, 'gene_signatures', 'pearson_perm_res')


# ---- Read gene signature


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
    annotations_file = glob.glob(
        {os.path.join(annotations_path, file_name)})[0]
    if annotations_file is None:
        raise ValueError(
            f'No metadata file named {file_name} could be found in {annotations_path}')

    # Read csv file and return df
    return pd.read_csv(annotations_file)


# add useful comments

# Make cell synonyms df
def cell_synonyms(cells_df, file_name, file_path):
    # Get metadata file
    file_name = 'cell_annotation_all.csv'
    cell_names_df = get_annotations(file_name, file_path)

    # Find all columns relevant to cellid
    pattern = re.compile('cellid')
    cell_columns = cell_names_df[[
        col for col in cell_names_df.columns if pattern.search(col)]]

    # Initialize DataFrame for cell synonyms
    cell_syn_df = pd.DataFrame(columns=['cell_id', 'cell_name'])

    for (index, cell_row) in cell_columns.iterrows():
        cell_id = cells_df[cells_df['name'] == cell_names_df['X'][index]]['id']
        synonyms = pd.Series(pd.unique(cell_row)).dropna()
        cell_syn_df = cell_syn_df.append(
            pd.DataFrame({'cell_id': cell_id, 'cell_name': synonyms}),
            ignore_index=True)

    # Add primary key
    cell_syn_df['id'] = cell_syn_df.index

    return cell_syn_df

    # using X
    # do I ever use curation?
    # do I need datasetid? it's in red..
    # lots unmatched -- fix


# very similar to cell_synonyms and drug_synonyms fxn.. should i condense, and how?
def tissue_synonyms(tissues_df, file_name, file_path):
    # Get metadata file
    file_name = 'cell_annotation_all.csv'
    tissues_metadata = get_annotations(file_name, file_path)

    # Find all columns relevant to tissueid
    pattern = re.compile('tissueid')
    tissue_cols = tissues_metadata[[
        col for col in tissues_metadata.columns if pattern.search(col)]]

    # Initialize DataFrame for tissue synonyms
    tissues_syn_df = pd.DataFrame(columns=['tissue_id', 'tissue_name'])

    for (index, tissue_row) in tissue_cols.iterrows():
        tissue_id = tissues_df[tissues_df['name'] ==
                               tissue_cols['unique.tissueid'][index]]['id']
        synonyms = pd.Series(pd.unique(tissue_row)).dropna()
        tissues_syn_df.append(
            pd.DataFrame({'tissue_id': tissue_id, 'tissue_name': synonyms}),
            ignore_index=True)

    # Add primary key
    tissues_syn_df['id'] = tissues_syn_df.index

    return tissues_syn_df

    # add datasetid???


# Make gene drugs df
def build_gene_drugs_df(gene_sig_df, genes_df, drugs_df, tissues_df):
    # Extract relevant columns
    gene_drugs_df = gene_sig_df[[
        'gene_id', 'drug', 'estimate', 'n', 'pvalue', 'df', 'fdr', 'tissue', 'mDataType']]
    # TODO - cannot find 'se', 'sens_stat' -- is one of these 'significant'???
    # TODO - cannot find 'level' ('lower'? 'upper'?)
    # TODO - cannot find 'drug_like_molecule', 'in_clinical_trials'

    # Join with genes_df
    gene_drugs_df = pd.merge(genes_df, gene_drugs_df,
                             left_on='name', right_on='gene_id', how='right')
    # Drop gene name columns after merge & rename 'id' column
    gene_drugs_df.drop(['name', 'gene_id'], axis='columns', inplace=True)
    gene_drugs_df = gene_drugs_df.rename(columns={"id": "gene_id"})

    # Join with drugs_df
    gene_drugs_df = pd.merge(drugs_df, gene_drugs_df,
                             left_on='name', right_on='drug', how='right')
    # TODO - this merges better on drugid than DRUG_NAME
    # Drop drug name columns after merge & rename 'id' column
    gene_drugs_df.drop(['name', 'drug'], axis='columns', inplace=True)
    gene_drugs_df = gene_drugs_df.rename(columns={"id": "drug_id"})

    # Join with tissues_df
    gene_drugs_df = pd.merge(gene_drugs_df, tissues_df,
                             left_on='tissue', right_on='name', how='left')
    # TODO - test this merge
    # Drop tissue name columns after merge & rename 'id' column
    gene_drugs_df.drop(['name', 'tissue'], axis='columns', inplace=True)
    gene_drugs_df = gene_drugs_df.rename(columns={"id": "tissue_id"})

    # dataset id..

    return gene_drugs_df


def build_datasets_cells_df(pset_dict, cell_df, dataset_id):
    datasets_cells_df = pd.merge(pset_dict['cell'][[
                                 'cellid']], cell_df, left_on='cellid', right_on='name', how='left')[['id']]
    datasets_cells_df.rename(columns={"id": "cell_id"})
    datasets_cells_df['dataset_id'] = dataset_id
    return datasets_cells_df


# TODO - needs to be faster
def build_dose_responses_df(pset_dict, experiment_df):
    # Get dose and response info from pset
    doses = pset_dict['sensitivity']['raw.Dose']
    responses = pset_dict['sensitivity']['raw.Viability']

    dose_responses_df = pd.DataFrame(columns=['exp_id', 'dose', 'response'])

    for (index, dose_row) in doses.iterrows():
        exp = dose_row['.exp_id']
        doses_df = pd.DataFrame(
            {'exp_id': exp, 'dose_num': dose_row.index[1:], 'dose': dose_row[1:]})

        response_row = responses[responses['.exp_id'] == exp].transpose()
        response_row = response_row.drop('.exp_id')
        responses_df = pd.DataFrame(
            {'dose_num': response_row.index, 'response': response_row[response_row.columns[0]]})

        dose_resp = pd.merge(doses_df, responses_df, on='dose_num')
        dose_resp.drop('dose_num', axis='columns', inplace=True)

        dose_responses_df = dose_responses_df.append(dose_resp)

    dose_responses_df = pd.merge(dose_responses_df, experiment_df, on='exp_id', how='left')[
        ['id', 'dose', 'response']]
    dose_responses_df.rename(columns={"id": "experiment_id"})

    return dose_responses_df

# ---- Build dose response table
# Chris' implementation
def build_dose_response_df(pset_dict, experiment_df):  # NOTE: database tables should always use singular names
    # Get dose and response info from pset
    dose = pset_dict['sensitivity']['raw.Dose']
    response = pset_dict['sensitivity']['raw.Viability']

    # rename columns so they can be coerced to int later
    rename_dict = {f'dose.{n}': str(n) for n in np.arange(1, dose.shape[0])}
    dose.rename(columns=rename_dict, inplace=True)
    response.rename(columns=rename_dict, inplace=True)

    # reshape the DataFrames using melt and pivot to go from 'wide' to 'long' or back, respectively
    # these are much faster than using Python loops because all of the looping is done in the Pandas C++ code
    dose = dose.melt(id_vars='.exp_id', value_name='dose', var_name='dose_id').dropna()
    dose['dose_id'] = dose.dose_id.astype('int')
    response = response.melt(id_vars='.exp_id', value_name='response', var_name='dose_id').dropna()
    response['dose_id'] = response.dose_id.astype('int')

    # set indexes for faster joins (~3x)
    dose.set_index(['.exp_id', 'dose_id'], inplace=True)
    response.set_index(['.exp_id', 'dose_id'], inplace=True)

    # because I set indexes on both table I don't need to specify on
    dose_response_df = pd.merge(dose, response, left_index=True, right_index=True).reset_index()

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
    # TODO - this join works better with DRUG_NAME rather than drugid
    drug_targets_df = pset_dict['drug'][['drugid', 'TARGET']]
    drug_targets_df = pd.merge(drug_df, drug_targets_df, left_on='name',
                               right_on='drugid', how='right')[['id', 'TARGET']]
    # Rename drug id column (FK)
    drug_targets_df.rename(columns={"id": "drug_id"})

    # Join with target df
    drug_targets_df = pd.merge(drug_targets_df, target_df, left_on='TARGET',
                               right_on='name', how='left')[['drug_id', 'id']]
    # Rename target id column (FK)
    drug_targets_df.rename(columns={"id": "target_id"})

    return drug_targets_df


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
    ## FIXME: was missing inplace argument, otherwise df.rename returns a copy
    experiments_df.rename(columns={"id": "drug_id"}, inplace=True)

    # Add dataset_id
    experiments_df['dataset_id'] = dataset_id

    # TODO - get rid of exp_id; check if col order matters

    return experiments_df


def build_annotation_dfs(pset_dict, gene_df, drug_df):
    """
    build out drug_annotations and gene_annotations dfs

    @param pset_dict: [`dict`]
    @param gene_df: [`DataFrame`]
    @param drug_df: [`DataFrame`]
    @return: [(`DataFrame`, `DataFrame`)]
    """
    # Make gene_annotations df
    # TODO - add gene_seq_start, gene_seq_end
    gene_annotations_df = pd.DataFrame(columns=['gene_id', 'symbol',
                                                'gene_seq_start', 'gene_seq_end'])
    # pset_dict['molecularProfiles']['rna']['rowData'][['Symbol']]

    # Make drug_annotations df
    drug_annotations_df = pset_dict['drug'][['rownames', 'smiles', 'inchikey',
                                             'cid', 'FDA']]
    drug_annotations_df.columns = [
        'name', 'smiles', 'inchikey', 'pubchem', 'fda_status']
    # Create foreign key to drug_df
    drug_annotations_df = pd.merge(
        drug_df, drug_annotations_df, on='name', how='right')
    # Drop name column once you've merged on it
    drug_annotations_df.drop('name', axis='columns', inplace=True)

    return gene_annotations_df, drug_annotations_df


def build_cell_target_dfs(pset_dict, tissue_df, gene_df):
    """
    ~ description ~ (build out cell and target dataframes)

    @param pset_dict: [`dict`]
    @param tissue_df: [`DataFrame`]
    @param gene_df: [`DataFrame`]
    @return: [(`DataFrame`, `DataFrame`)]
    """
    # Make the cells df TODO
    cell_df = pd.merge(pset_dict['cell'], tissue_df, left_on='tissueid', right_on='name',
                       how='left')[['cellid', 'id']]
    cell_df.columns = ['name', 'tissue_id']

    # Make the targets df -- NEED TO ADD GENE ID TODO - how to relate target to genes??
    # target_df = pd.merge(pset_dict['drug'][['TARGET']], gene_df)[
    #     ['TARGET', 'id']]
    # target_df.columns = ['name', 'gene_id']

    return cell_df #, target_df


def build_primary_tables(pset_dict):
    """
    Build the tissues, drug, cell, dataset DataFrames, to be used to construct the 
    experiments DataFrame.

    @param pset_dict: [`dict`]
    @return: [(`DataFrame`, `DataFrame`, `DataFrame`)]
    """
    # Make tissues df; TODO - check that you should get it from 'cell' rather than 'curation'
    tissues = pd.Series(pd.unique(pset_dict['cell']['tissueid']))
    tissue_df = pd.DataFrame({'id': tissues.index, 'name': tissues})

    # Make drugs df; TODO - confirm whether to use drugid (?) or DRUG_NAME (targets)
    drugs = pd.Series(pd.unique(pset_dict['drug']['drugid']))
    drug_df = pd.DataFrame({'id': drugs.index, 'name': drugs})

    # Makes genes df; TODO - check that you're getting the correct ID
    genes = pd.Series(pd.unique(
        pset_dict['molecularProfiles']['rna']['rowData']['EnsemblGeneId']))
    gene_df = pd.DataFrame({'id': genes.index, 'name': genes})

    return tissue_df, drug_df, gene_df


if __name__ == "__main__":
    # pset_names = ['GDSC_v1', 'GDSC_v2', 'gCSI',
    #               'FIMM', 'CTRPv2', 'CCLE', 'GRAY', 'UHNBreast']

    pset_names = ['GDSC_v1']

    # TODO - turn this all into more functions, and condense into more concise code

    # Convert pset files to nested dictionaries
    pset_dicts = []
    for pset in pset_names:
        pset_df = read_pset(pset_name, file_path)
        pset_dict = pset_df_to_nested_dict(pset_df)
        pset_dicts.append(pset_dict)

    # Initialize primary DataFrames
    dataset_df = pd.DataFrame(columns=['id', 'name'])
    tissue_df = pd.DataFrame(columns=['id', 'name'])
    drug_df = pd.DataFrame(columns=['id', 'name'])
    gene_df = pd.DataFrame(columns=['id', 'name'])

    # Build out primary DataFrames
    for i in range(len(pset_dicts)):
        dataset_df = dataset_df.append(
            {'name': pset_names[i]}, ignore_index=True)
        tissues, drugs, genes = build_primary_tables(pset_dicts[i])
        tissue_df = tissue_df.append(tissues, ignore_index=True)
        drug_df = drug_df.append(drugs, ignore_index=True)
        gene_df = gene_df.append(genes, ignore_index=True)

    # Remove duplicate rows & reindex ID columns
    for df in [tissue_df, drug_df, gene_df]:
        df.drop_duplicates()
        df['id'] = df.index
    dataset_df['id'] = dataset_df.index

    # Initialize cell, target & annotation DataFrames
    cell_df = pd.DataFrame(columns=['id', 'name', 'tissue_id'])
    target_df = pd.DataFrame(columns=['id', 'name', 'gene_id'])
    gene_annotations_df = pd.DataFrame(
        columns=['gene_id', 'symbol', 'gene_seq_start', 'gene_seq_end'])
    drug_annotations_df = pd.DataFrame(
        columns=['drug_id', 'smiles', 'inchikey', 'pubchem', 'fda_status'])

    # Use tissue, drug, gene dfs to build out cell, target and annotation DataFrames
    for i in range(len(pset_dicts)):
        dataset_id = dataset_df[dataset_df['name'] ==
                                pset_names[i]]['id']  # TODO - check this
        cells, targets = build_cell_target_dfs(pset_dicts[i],
                                               tissue_df, gene_df)
        cell_df = cell_df.append(cells, ignore_index=True)
        target_df = target_df.append(targets, ignore_index=True)
        gene_annotations, drug_annotations = build_annotation_dfs(
            pset_dict, gene_df, drug_df)
        gene_annotations_df = gene_annotations_df.append(
            gene_annotations, ignore_index=True)
        drug_annotations_df = drug_annotations_df.append(
            drug_annotations, ignore_index=True)

    # Remove duplicates from cell and target dfs
    cell_df.drop_duplicates()
    target_df.drop_duplicates()

    # Update primary keys
    for df in [cell_df, target_df, gene_annotations_df, drug_annotations_df]:
        df['id'] = df.index + 1  # NOTE: database table indexing should start at 1; we use zero for special cases
                                 # such as missing PKs (in which case you just match to PK 0)

    # Initialize datasets_cells, experiments, mol_cells DataFrames
    datasets_cells_df = pd.DataFrame(columns=['id', 'dataset_id', 'cell_id'])
    experiments_df = pd.DataFrame(
        columns=['id', 'cell_id', 'drug_id', 'dataset_id', 'tissue_id'])
    mol_cells_df = pd.DataFrame(
        columns=['id', 'cell_id', 'dataset_id', 'mDataType, num_prof'])

    # Use cell df to build datasets_cells, experiments, mol_cells dfs
    for i in range(len(pset_dicts)):
        dataset_id = dataset_df[dataset_df['name'] ==
                                pset_names[i]]['id']  # TODO - check this
        datasets_cells_df = datasets_cells_df.append(build_datasets_cells_df(pset_dicts[i], cell_df, dataset_id))
        experiments_df = experiments_df.append(
            build_experiment_df(pset_dicts[i], cell_df, drug_df, dataset_id))
        # TODO - mol_cells dfs

    # Add primary key
    for df in [datasets_cells_df, experiments_df, mol_cells_df]:
        df['id'] = df.index

    # dose_response_df = pd.DataFrame(
    #    columns=['id', 'experiment_id', 'dose', 'response'])
    # dose_response_df = dose_response_df.append(
    #        dose_responses, ignore_index=True)

    # Initialize drug_targets DataFrame
    drug_targets_df = pd.DataFrame(columns=['id', 'drug_id', 'target_id'])

    # Use drug_df and target_df to build drug_targets df
    for pset_dict in pset_dicts:
        drug_targets_df = drug_targets_df.append(
            build_drug_targets_df(pset_dict, drug_df, target_df))

    # Add primary key
    drug_targets_df['id'] = drug_targets_df.index

    # Build out gene_drugs df

    # Build out synonyms DataFrames from metadata


# More comments -- confused about datasets_cells table? what is it for?
# Where is clinical_trials info stored?
# oncotrees, cellosaurus, dataset_statistics, profiles, mol_cells
