import glob
import os
import re
import numpy as np
import pandas as pd
from datatable import dt, fread, iread, join, by, rbind, cbind, f

data_dir = os.path.join('data', 'procdata')
output_dir = os.path.join('data', 'demo')

# TODO: change printed errors to actual errors

def load_table(name, data_dir):
    """
    Load all PSet tables with name into a datatable and reindex the rows.

    @name: [`string`] The name of the table
    @data_dir: [`string`] File path to the directory with all PSet tables
    @return: [`datatable.Frame`] A datatable containing all rows from all PSets
    """
    # Get all files
    files = glob.glob(os.path.join(data_dir, '**', f'*{name}.csv'))
    # Filter so that file path are '{data_dir}/{pset}/{pset}_{name}.csv'
    files = [file_name for file_name in files if re.search(data_dir + r'/(\w+)/\1_' + name + '.csv$', file_name)]
    # Read and concatenate tables
    df = rbind(*iread(files, sep=','))
    # Drop duplicates
    # (groups by all columns and selects only the first row from each group)
    df = df[0, :, by(df.names)]

    return df


def rename_and_key(df, join_col, og_col='name'):
    """
    Prepare df to be joined with other tables by renaming the column
    on which it will be joined and by keying it.

    @join_col: [`string`] The name of the join column in other tables
                            (ex. 'tissue_id', 'cell_id', etc.)
    @og_col: [`string`] The name of the join column in the join table
    """
    # Rename primary key to match foreign key name (necessary for joins)
    df.names = {og_col: join_col}
    # Only select necessary rows
    df = df[:, ['id', join_col]]
    # Set the key
    df.key = join_col
    return df # Not necessary? df passed by reference


def join_tables(df1, df2, join_col):
    """
    Join df2 and df1 based on join_col.

    @df1: [`datatable.Frame`] The datatable with the foreign key
    @df2: [`datatable.Frame`] The join table (ex. tissue datatable)
    @join_col: [`string`] The name of the columns on which the tables
                            will be joined (ex. 'tissue_id')
    """
    if (join_col not in df1.names) or (join_col not in df2.names):
        print(f'{join_col} is missing from one or both of the datatables passed!',
              'Make sure you have prepared df2 using rename_and_key().')
        return None

    # Join tables, then rename the join col and drop it
    df = df1[:, :, join(df2)]
    df.names = {join_col: 'drop', 'id': join_col}
    del df[:, 'drop']
    return df


def index_and_write(df, name, output_dir):
    # Index datatable
    df = cbind(dt.Frame(id=np.arange(df.nrows) + 1), df)
    # Write to .csv
    df.to_csv(os.path.join(output_dir, f'{name}.csv'))
    return df


def load_join_write(name, data_dir, output_dir, foreign_keys=[], join_dfs=None):
    df = load_table(name, data_dir)
    if foreign_keys and not join_tables:
        print(f'ERROR: The {name} table has foreign keys {foreign_keys}'
                'but you have not passed any join_tables.')
        return None

    for fk in foreign_keys:
        df = join_tables(df, join_dfs[fk], fk+'_id')
    
    df = index_and_write(df, name, output_dir)
    
    return df


# TODO - similar names to build_pset_tables, should I change it to be clearer?
def build_primary_tables(data_dir, output_dir):
    """
    Build all the primary tables, i.e., tables that require no joins,
    and return them in a dictionary.

    @param data_dir: [`string`] The file path to read the PSet tables
    @param output_dir: [`string`] The file path to write the final tables
    @return: [`dict(string: datatable.Frame)`] A dictionary of all the primary
                                                tables, with names as keys
    """
    # Load, concatenate, and write primary tables to disk
    tissue_df = load_join_write('tissue', data_dir, output_dir)
    drug_df = load_join_write('drug', data_dir, output_dir)
    gene_df = load_join_write('gene', data_dir, output_dir)
    dataset_df = load_join_write('dataset', data_dir, output_dir)

    # Transform tables to be used for joins
    dfs = {}
    dfs['tissue'] = rename_and_key(tissue_df, 'tissue_id')
    dfs['drug'] = rename_and_key(drug_df, 'drug_id')
    dfs['gene'] = rename_and_key(gene_df, 'gene_id')
    dfs['dataset'] = rename_and_key(dataset_df, 'dataset_id')
    return dfs


def build_secondary_tables(join_dfs, data_dir, output_dir):
    """
    Build all secondary tables, i.e., all tables that have foreign keys corresponding
    to primary keys of primary tables. The function reads PSet tables from 
    data_dir, concatenates and joins them with tables from join_dfs, and 
    writes them to output_dir.

    @param join_dfs: [`dict(string: datatable.Frame)`] A dictionary of all the primary
                                                    tables, with names as keys
    @param data_dir: [`string`] The file path to read the PSet tables
    @param output_dir: [`string`] The file path to write the final tables
    @return: [`dict(string: datatable.Frame)`] The updated dictionary of join tables
    """
    # Build cell table and add to join_dfs dictionary
    cell_df = load_join_write('cell', data_dir, output_dir, ['tissue'], join_dfs)
    join_dfs['cell'] = rename_and_key(cell_df, 'cell_id')

    # Build annotation tables
    load_join_write('drug_annotation', data_dir, output_dir, ['drug'], join_dfs)
    load_join_write('gene_annotation', data_dir, output_dir, ['gene'], join_dfs)

    # Build all other secondary tables
    load_join_write('dataset_cell', data_dir, output_dir, ['dataset', 'cell'], join_dfs)
    load_join_write('mol_cell', data_dir, output_dir, ['cell', 'dataset'], join_dfs)
    # mol_cells has Kallisto. not sure why. from CTRPv2
    load_join_write('dataset_statistics', data_dir, output_dir, ['dataset'], join_dfs)
    load_join_write('gene_drug', data_dir, output_dir, ['gene', 'drug', 'dataset', 'tissue'], join_dfs)

    return join_dfs


def build_experiment_tables(join_dfs, data_dir, output_dir):
    """
    Load and process experiment table, then use it to build the dose response
    and profile tables. Drop the 'name' column from the experiment table before
    writing to a CSV.

    @param join_dfs: [`dict(string: datatable.Frame)`]
    @param data_dir: [`string`] The file path to the PSet tables
    @param output_dir: [`string`] The file path to the final tables
    @return: [`None`]
    """
    # Load all experiments from PSets
    experiment_df = load_join_write('experiment', data_dir, output_dir, ['cell', 'drug', 'dataset', 'tissue'], join_dfs) 
    
    #concat_tables(load_pset_tables('experiment', data_dir))

    # Join experiment table and with primary tables
    #for fk in ['cell', 'drug', 'dataset', 'tissue']:
    #    join_df = join_dfs[fk]
    #    experiment_df = safe_merge(experiment_df, join_df, f'{fk}_id')

    # Reindex experiment table
    #experiment_df = reindex_table(experiment_df)
    join_dfs['experiment'] = rename_and_key(experiment_df, 'experiment_id')
    # Load and concatenate dose response and profile tables, join with experiment table
    #join_dict = {'experiment': experiment_df}
    load_join_write('dose_response', data_dir, output_dir, ['experiment'], join_dfs)
    load_join_write('profile', data_dir, output_dir, ['experiment'], join_dfs)

    # Drop experiment name from table after joins and write to disk
    del experiment_df[:, 'name']
    experiment_df.to_csv(os.path.join(output_dir, 'experiment.csv'))

    return join_dfs
#1019495 rows

# TODO - better name
def build_final_tables(data_dir, output_dir):
    join_dfs = build_primary_tables(data_dir, output_dir)
    join_dfs = build_secondary_tables(join_dfs, data_dir, output_dir)
    #build_experiment_tables(join_dfs, data_dir, output_dir)
    return join_dfs


# Duplicate experiment_ids
"""
   | name                      count
-- + ------------------------  -----
 0 | drugid_AZD6244_AU565          2
 1 | drugid_AZD6244_HCC1187        2
 2 | drugid_AZD6244_HCC1806        2
 3 | drugid_AZD6244_HCC1954        2
 4 | drugid_AZD6244_MCF7           2
 5 | drugid_Erlotinib_HCC70        2
 6 | drugid_Nilotinib_AU565        2
 7 | drugid_Nilotinib_HCC1187      2
 8 | drugid_Nilotinib_HCC1806      2
 9 | drugid_Nilotinib_HCC1954      2
10 | drugid_Nilotinib_MCF7         2
11 | drugid_Sorafenib_HCC70        2
"""