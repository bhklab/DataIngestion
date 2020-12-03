import glob
import os
import re
import numpy as np
import pandas as pd
import dask.dataframe as dd


read_file_path = os.path.join('data', 'procdata')
write_file_path = os.path.join('data', 'latest')
psets = ['CTRPv2', 'FIMM', 'gCSI', 'GDSC_v1', 'GDSC_v2', 'GRAY', 'UHNBreast']


def load_tables(name, file_path, psets=['CTRPv2', 'FIMM', 'gCSI', 'GDSC_v1', 'GDSC_v2', 'GRAY', 'UHNBreast']):
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
                dfs.append(dd.read_csv(f'{pset_path}/{pset}_{name}-*-*.csv'))
    
    return dfs


def concat_tables(df_list, chunksize=250000):
    """
    Concatenate all DataFrames in df_list. Also drops duplicate rows
    and resets index to start at 1 and increment each row.
    
    :param: df_list list[DataFrame] List of Dask DataFrames to be concatenated
    :return: DataFrame A single DataFrame containing all rows of the dfs in df_list
    """
    # Dask concat all dfs in df_list
    df = dd.concat(df_list, interleave_partitions=True)
        #TODO - when do i really need to compute it?
    
    # for now, until you rewrite the tables and fix the tree
    if 'Unnamed: 0' in df.columns:
        df = df.drop('Unnamed: 0', axis=1)
    
    # Drop duplicate rows (inplace parameter not supported in Dask)
    df = df.drop_duplicates()
    
    return df


def write_table(df, name, file_path):
    """
    Write table to file_path
    """
    # Reindex
    df = df.reset_index(drop=True) # This only works when it's one partition (TODO)
    df.index += 1
    df.index = df.index.rename('id', inplace=True)

    df_path = os.path.join(file_path, name)
    
    # Check that directory exists
    if not os.path.exists(df_path):
        os.mkdir(df_path)

    # TODO - empty directory first
    
    # Write to CSV
    dd.to_csv(df, os.path.join(df_path, f'{name}-*.csv'))


def load_concat_write(name, read_file_path, write_file_path):
    df_list = load_tables(name, read_file_path)
    df = concat_tables(df_list)
    write_table(df, name, write_file_path)
    # TODO - how to remove sthg from memory to clear up space once you write to disk
    # should clear automatically after function terminates


# I shouldn't do this because i'm making it run longer?
def load_join_table(name, file_path):
    df_files = glob.glob(os.path.join(file_path, name))
    if len(df_files) == 0:
        print('ERROR cannot find the join table')
        #raise ValueError
    
    return dd.read_csv(f'{os.path.join(file_path, name, name)}-*.csv')


def safe_merge(df1, df2, fk_name, right_on='name'):
    """
    Will always perform left join
    By default will join on 'name' in df2
    Will always delete join columns
    Will rename id column from df2 to fk_name
    Will check that left join worked well
    """
    # Make copy of relevant cols of df2 so you don't have to drop all the other ones after merge
    df2 = df2[['id', right_on]].copy()
    # Rename df2 FK column to avoid naming conflicts (since df1 will probably have 'name' col too)
    df2 = df2.rename(columns={right_on: f'{fk_name}_{right_on}'})
    right_on = f'{fk_name}_{right_on}'
    
    # Join DataFrames
    # NOTE: this automatically resets the index, so it starts at 0 and doesn't have the name 'id'
    join_df = dd.merge(df1, df2, left_on=fk_name, right_on=right_on, how='left')
    
    # Remove join columns
    join_df = join_df.drop(columns=[fk_name, right_on])
    # Rename df2 id col to fk_name
    join_df = join_df.rename(columns={'id': fk_name})
    
    # Check if any rows did not join properly with dataset
            # can't do this with variable? ---  join_df.query('dataset_id.isna()')
            #TODO - consider using validate param instead; consider throwing an error
    if join_df[join_df[fk_name].isna()].shape[0].compute() > 0:
        print(f'ERROR - some rows did not join to a {fk_name}!')
    
    return join_df


def load_join_write(name, fks, read_file_path, write_file_path):
    # Load and concatenate all 'name' tables for all PSets
    df_list = load_tables(name, read_file_path)
    df = concat_tables(df_list)

    # Load each join table and join with it
    for fk in fks:
        join_df = load_join_table(fk, write_file_path)
        df = safe_merge(df, join_df, f'{fk}_id')

    # Write to disk
    write_table(df, name, write_file_path)


# Primary tables can just be concatenated and written to disk
load_concat_write('tissue', read_file_path, write_file_path)
load_concat_write('drug', read_file_path, write_file_path)
load_concat_write('gene', read_file_path, write_file_path)
load_concat_write('dataset', read_file_path, write_file_path) # This seems extra

# Other tables need to be joined to get their foreign keys
load_join_write('cell', ['tissue'], read_file_path, write_file_path)
load_join_write('datasets_cells', ['dataset', 'cell'], read_file_path, write_file_path)
load_join_write('drug_annotations', ['drug'], read_file_path, write_file_path)
load_join_write('gene_annotations', ['gene'], read_file_path, write_file_path)
load_join_write('experiments', ['cell', 'drug', 'dataset', 'tissue'], read_file_path, write_file_path)
#load_join_write('dose_responses', ['experiment'], read_file_path, write_file_path)
#load_join_write('profiles', ['experiment'], read_file_path, write_file_path)
load_join_write('mol_cells', ['cell', 'dataset'], read_file_path, write_file_path)

# TODO - joins with experiment_id need special handling
# TODO - gene_drugs, dose_responses are multiple partitions; need to make sure these joins work
# and that the indices are updated appropriately (doesn't reset with each partition)
# TODO - need to load non-pset-specific tables and join those

# del or if it's inside a fxn then it's fine