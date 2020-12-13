import glob
import os
import re
import numpy as np
import pandas as pd
import dask.dataframe as dd


read_file_path = os.path.join('data', 'procdata')
write_file_path = os.path.join('data', 'latest')
psets = ['CTRPv2', 'FIMM', 'gCSI', 'GDSC_v1',
         'GDSC_v2', 'GRAY', 'UHNBreast', 'CCLE']


def load_pset_tables(name, file_path, psets=['CTRPv2', 'FIMM', 'gCSI', 'GDSC_v1', 'GDSC_v2', 'GRAY', 'UHNBreast', 'CCLE']):
    """
    Given the name of a table and a file_path to all PSet tables,
    read all the tables into DataFrames and return a list.

    :param: name string - Name of table
    :param: file_path string - File path to all processed data (all PSet tables)
    :return: list[DataFrame] - List of Dask DataFrames, one from each PSet
    """
    dfs = []

    for pset in psets:
        # Check that this pset has this table
        pset_path = os.path.join(file_path, pset, name)
        if os.path.exists(pset_path):
            if len(os.listdir(pset_path)) == 1:
                # Written as a single CSV file with pandas read_csv
                dfs.append(dd.read_csv(f'{pset_path}/{pset}_{name}.csv'))
            else:
                # Written as many CSVs with Dask
                dfs.append(dd.read_csv(f'{pset_path}/{pset}_{name}-*.csv'))

    return dfs


def concat_tables(df_list, chunksize=500000):
    """
    Concatenate all DataFrames in df_list. Also drops duplicate rows
    and resets index to start at 1 and increment each row.

    :param: df_list list[DataFrame] List of Dask DataFrames to be concatenated
    :return: DataFrame A single DataFrame containing all rows of the dfs in df_list
    """
    # Dask concat all dfs in df_list
    df = dd.concat(df_list)

    # Drop id column (old index)
    if 'id' in df.columns:
        df = df.drop('id', axis=1)

    # Drop duplicate rows
    df = df.drop_duplicates()

    # Repartition into partitions, each with chunksize rows
    num_rows = df.shape[0].compute()  # This may take a while
    df = df.repartition(npartitions=(num_rows // chunksize + 1))

    return df


def reindex_table(df):
    """
    Reindex table so that its primary key follows SQL standards,
    i.e., starts at 1 and doesn't restart with each partition.
    Also, name the index 'id'.

    @param df: [`DataFrame`] A Dask DataFrame
    @return: [`DataFrame`] The same DataFrame, but reindexed
    """
    # Reindex
    if df.npartitions > 1:
        n = df.npartitions
        df = df.repartition(npartitions=1)
        df = df.reset_index(drop=True)
        df = df.repartition(npartitions=n)
    else:
        df = df.reset_index(drop=True)

    # Increment index by 1 and rename
    df.index += 1
    df.index = df.index.rename('id')

    return df


def write_table(df, name, file_path):
    """
    Write table to file_path
    """
    df_path = os.path.join(file_path, name)
    # Check that directory exists
    if not os.path.exists(df_path):
        os.mkdir(df_path)
    else:
        # Clear all old files
        for f in os.listdir(df_path):
            os.remove(os.path.join(df_path, f))

    # Write to CSV
    print(f'Writing {name} table to disk...')
    dd.to_csv(df, os.path.join(df_path, f'{name}-*.csv'))


def load_concat_write(name, read_file_path, write_file_path):
    """
    Given the name of the table, and read and write file path, read all
    "name" tables (one for each PSet) from the read_file_path, concatenate
    them, reindex them, and then write them to the write_file_path.

    @param name: [`string`] The name of the table
    @param read_file_path: [`string`] The file path containing the PSet tables
    @param write_file_path: [`string`] The file path containing the final tables
    @return: [`DataFrame`] The final DataFrame that was written to disk
    """
    df_list = load_pset_tables(name, read_file_path)
    df = concat_tables(df_list)
    df = reindex_table(df)
    write_table(df, name, write_file_path)
    return df


# I shouldn't do this because i'm making it run longer?
def load_join_table(name, file_path):
    df_files = glob.glob(os.path.join(file_path, name))
    if len(df_files) == 0:
        print('ERROR cannot find the join table')
        #raise ValueError

    return dd.read_csv(f'{os.path.join(file_path, name, name)}-*.csv')


def safe_merge(df1, df2, fk_name, right_on='name'):
    """
    Given two Dask DataFrames, performs a left join to add a foreign key from
    df2 to df1. By default, the function will join on 'name' in df2 (ex. referring
    to the gene name, drug name, etc.) and on fk_name in df1.

    Once the join occurs, the joining columns will be dropped, and only the forein
    ID column will remain, which will be renamed fk_name. safe_merge also prints a 
    warning if not all rows in df1 were matched to a record in df2.

    @param df1: [`DataFrame`] The DataFrame acquiring a foreign key
    @param df2: [`DataFrame`] The DataFrame providing the foreign key (its id column)
    @param fk_name: [`string`] The name of the column in df1 that will be joined on;
                                also the name of the foreign key once the join is complete
    @param right_on: [`string`] Optional: set to 'name' by default, but can be changed 
                                in special cases
    @return: [`DataFrame`] The joined DataFrame (df1 but with a foreign key corresponding
                            to the primary key of df2)
    """
    # Make copy of relevant cols of df2 so you don't have to drop all the other ones after merge
    df2 = df2[['id', right_on]].copy()
    # Rename df2 FK column to avoid naming conflicts (since df1 will probably have 'name' col too)
    df2 = df2.rename(columns={right_on: f'{fk_name}_{right_on}'})
    right_on = f'{fk_name}_{right_on}'

    # Join DataFrames
    join_df = dd.merge(df1, df2, left_on=fk_name,
                       right_on=right_on, how='left')

    # Remove join columns
    join_df = join_df.drop(columns=[fk_name, right_on])
    # Rename df2 id col to fk_name
    join_df = join_df.rename(columns={'id': fk_name})
    # Cast FK column as int
    join_df[fk_name] = join_df[fk_name].astype('int')

    # Check if any rows did not join properly with dataset
    # TODO - consider using validate param instead; consider throwing an error
    if join_df[join_df[fk_name].isna()].shape[0].compute() > 0:
        print(f'ERROR - some rows did not join to a {fk_name}!')

    return join_df


def load_join_write(name, fks, join_dfs, read_file_path, write_file_path):
    """
    Given the name to a table and a list of its foreign keys, load
    all PSet tables of that name from read_file_path, concatenate 
    them, build all the foreign keys by joining to the necessary 
    tables, and write the final DataFrame to write_file_path.

    @param name: [`string`] The name of the table
    @param fks: [`list(string)`] A list of tables to be joined with
    @param join_dfs: [`dict(string: DataFrame)`] A dictionary of tables
                                    to be joined to to build foreign keys
    @param read_file_path: [`string`] The file path to the PSet tables
    @param write_file_path: [`string`] The file path to the final tables
    @return: [`DataFrame`] The final DataFrame that was written to disk
    """
    # Load and concatenate all 'name' tables for all PSets
    df_list = load_pset_tables(name, read_file_path)
    df = concat_tables(df_list)

    # Build all foreign keys by joining with primary tables
    for fk in fks:
        join_df = join_dfs[fk]
        df = safe_merge(df, join_df, f'{fk}_id')

    # Reindex the final table and write to disk
    df = reindex_table(df)
    write_table(df, name, write_file_path)

    return df


# TODO - similar names to build_pset_tables, should I change it to be clearer?
def build_primary_tables(read_file_path, write_file_path):
    """
    Build all the primary tables, i.e., tables that require no joins,
    and return them in a dictionary.

    @param read_file_path: [`string`] The file path to the PSet tables
    @param write_file_path: [`string`] The file path to the final tables
    @return: [`dict(string: DataFrame)`] A dictionary of all the primary
                                        tables, with names as keys
    """
    # Load, concatenate, and write primary tables to disk
    tissue_df = load_concat_write('tissue', read_file_path, write_file_path)
    drug_df = load_concat_write('drug', read_file_path, write_file_path)
    gene_df = load_concat_write('gene', read_file_path, write_file_path)
    dataset_df = load_concat_write('dataset', read_file_path, write_file_path)

    return {'tissue': tissue_df, 'drug': drug_df, 'gene': gene_df, 'dataset': dataset_df}


def build_secondary_tables(join_dfs, read_file_path, write_file_path):
    """
    Build all secondary tables, i.e., all tables that have foreign keys corresponding
    to primary keys of primary tables. The function reads PSet tables from 
    read_file_path, concatenates and joins them with tables from primary_dfs, and 
    writes them to write_file_path.

    @param join_dfs: [`dict(string: DataFrame)`] A dictionary of all the primary
                                                    tables, with names as keys
    @param read_file_path: [`string`] The file path to the PSet tables
    @param write_file_path: [`string`] The file path to the final tables
    @return: [`dict(string: DataFrame)`] The updated dictionary of tables to
                                        join with (includes cell table)
    """
    # Build cell table and add to join_dfs dictionary
    cell_df = load_join_write('cell', ['tissue'], join_dfs, read_file_path, write_file_path)
    join_dfs['cell'] = cell_df

    # Build annotation tables
    load_join_write('drug_annotation', ['drug'], join_dfs,
                    read_file_path, write_file_path)
    load_join_write('gene_annotation', ['gene'], join_dfs, read_file_path, write_file_path)
    # symbol col causing issues with dtype; i think because lots of vals missing

    # Build all other secondary tables
    load_join_write('dataset_cell', ['dataset', 'cell'], join_dfs,
                        read_file_path, write_file_path)
    load_join_write('mol_cell', ['cell', 'dataset'], join_dfs,
                    read_file_path, write_file_path)
    # mol_cells has Kallisto. not sure why. from CTRPv2
    load_join_write('gene_drug', ['gene', 'drug', 'dataset', 'tissue'], join_dfs,
                        read_file_path, write_file_path)

    return join_dfs


def build_experiment_tables(join_dfs, read_file_path, write_file_path):
    """
    ??????

    @param join_dfs: [`dict(string: DataFrame)`]
    @param read_file_path: [`string`] The file path to the PSet tables
    @param write_file_path: [`string`] The file path to the final tables
    @return: [`None`]
    """
    # Load all experiments from PSets
    experiment_df = concat_tables(load_pset_tables('experiment', read_file_path))

    # Join experiment table and with primary tables
    for fk in ['cell', 'drug', 'dataset', 'tissue']:
        join_df = join_dfs[fk]
        experiment_df = safe_merge(experiment_df, join_df, f'{fk}_id')

    # Reindex experiment table
    experiment_df = reindex_table(experiment_df)

    # Load and concatenate dose response and profile tables, join with experiment table
    join_dict = {'experiment': experiment_df}
    load_join_write('dose_response', ['experiment'], join_dict, read_file_path, write_file_path)
    load_join_write('profile', ['experiment'], join_dict, read_file_path, write_file_path)

    # Drop experiment name from table after joins and write to disk
    experiment_df = experiment_df.drop(columns=['name'])
    write_table(experiment_df, 'experiment', write_file_path)
    # TODO: test this!!!


# TODO - better name
def build_final_tables(read_file_path, write_file_path):
    join_dfs = build_primary_tables(read_file_path, write_file_path)
    join_dfs = build_secondary_tables(join_dfs, read_file_path, write_file_path)
    build_experiment_tables(join_dfs, read_file_path, write_file_path)
    return join_dfs
