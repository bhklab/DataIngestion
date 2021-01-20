import pandas as pd
import numpy as np
from PharmacoDI.build_primary_pset_tables import build_cell_df

def build_experiment_tables(pset_dict, pset_name, cell_df=None):
    """
    Build the experiment, dose response, and profile tables for a PSet 
    and return them in a dictionary, with table names as the keys.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet
    @param cell_df: [`DataFrame`] A table of all the cells in the PSet and their tissues
    @return: [`dict`] A dictionary of experiment-related tables
    """
    pset_dfs = {}
    pset_dfs['experiment'] = build_experiment_df(pset_dict, pset_name, cell_df)
    pset_dfs['dose_response'] = build_dose_response_df(pset_dict, pset_name)
    pset_dfs['profile'] = build_profile_df(pset_dict)
    return pset_dfs


def build_experiment_df(pset_dict, pset_name, cell_df=None):
    """
    Build a table with all the experiments of a PSet.

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @param pset_name: [`string`] The name of the PSet 
    @param cell_df: [`DataFrame`] A table of all the cells in the PSet and their tissues
    @return: [`DataFrame`] A table containing all experiments in the dataset
    """
    # Build cell_df if not build already
    if not cell_df:
        cell_df = build_cell_df(pset_dict)

    # Extract relelvant experiment columns
    experiment_df = pset_dict['sensitivity']['info'][[
        '.rownames', 'cellid', 'drugid']].copy()

    # Rename columns
    experiment_df.rename(
        columns={'.rownames': 'experiment_id', 'cellid': 'cell_id', 'drugid': 'drug_id'}, inplace=True)

    # Add datset_id column
    experiment_df['dataset_id'] = pset_name

    # Add tissue_id column by joining with cells_df
    experiment_df = pd.merge(experiment_df, cell_df[['name', 'tissue_id']],
                             left_on='cell_id', right_on='name', how='left')
    experiment_df = experiment_df[[
        'experiment_id', 'cell_id', 'drug_id', 'dataset_id', 'tissue_id']]

    experiment_df.rename(columns={'experiment_id': 'name'}, inplace=True)

    return experiment_df


# TODO:: Do Python functions pass my reference or copy by default?
# If pass by reference may need to make a copy before using inplace=TRUE argument
# to prevent modifying the original experiments table
def build_dose_response_df(pset_dict, pset_name):
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

    # Rename columns so they can be coerced to int later
    dose_pattern = dose.columns[1][:-1]
    rename_dict = {f'{dose_pattern}{n}': str(
        n) for n in np.arange(1, dose.shape[0])}
    dose.rename(columns=rename_dict, inplace=True)
    response.rename(columns=rename_dict, inplace=True)

    # Reshape the DataFrames using melt to go from 'wide' to 'long'
    dose = dose.melt(id_vars='.exp_id', value_name='dose',
                     var_name='dose_id').dropna()
    dose['dose_id'] = dose.dose_id.astype('int')
    response = response.melt(
        id_vars='.exp_id', value_name='response', var_name='dose_id').dropna()
    response['dose_id'] = response.dose_id.astype('int')

    # Set indices for faster joins (~3x)
    dose.set_index(['.exp_id', 'dose_id'], inplace=True)
    response.set_index(['.exp_id', 'dose_id'], inplace=True)

    # Join on index
    dose_response_df = pd.merge(
        dose, response, left_index=True, right_index=True).reset_index()
    dose_response_df.rename(columns={'.exp_id': 'experiment_id'}, inplace=True)
    dose_response_df.drop(columns=['dose_id'], inplace=True)

    # Add dataset_id for joins
    dose_response_df['dataset_id'] = pset_name

    return dose_response_df
    

def build_profile_df(pset_dict):
    """
    TODO: ask Chris

    @param pset_dict: [`dict`] A nested dictionary containing all tables in the PSet
    @return: [`DataFrame`] A table containing all statistics for each profile in the PSet (?)
    """
    # Get profiles info
    if 'E_inf' in pset_dict['sensitivity']['profiles'].columns:
        profile_df = pset_dict['sensitivity']['profiles'][[
            '.rownames', 'aac_recomputed', 'ic50_recomputed', 'HS', 'E_inf', 'EC50']].copy()
        profile_df.rename(columns={'.rownames': 'experiment_id', 'aac_recomputed': 'AAC',
                                   'ic50_recomputed': 'IC50', 'E_inf': 'Einf'}, inplace=True)
    else:
        profile_df = pset_dict['sensitivity']['profiles'][[
            '.rownames', 'aac_recomputed', 'ic50_recomputed', 'slope_recomputed', 'einf', 'ec50']].copy()
        profile_df.rename(columns={'.rownames': 'experiment_id', 'aac_recomputed': 'AAC', 'slope_recomputed': 'HS',
                                   'ic50_recomputed': 'IC50', 'einf': 'Einf', 'ec50': 'EC50'}, inplace=True)

    # Add DSS columns - TODO get these values when they are eventually computed
    profile_df['DSS1'] = np.nan
    profile_df['DSS2'] = np.nan
    profile_df['DSS3'] = np.nan

    return profile_df[['experiment_id', 'HS', 'Einf', 'EC50', 'AAC', 'IC50', 'DSS1', 'DSS2', 'DSS3']]

