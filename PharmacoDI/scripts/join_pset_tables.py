import glob
import os
import re
import numpy as np
import pandas as pd
from datatable import dt, fread, iread, join, by, rbind, cbind, f
import time

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


def load_join_write(name, data_dir, output_dir, foreign_keys=[], join_dfs=None):
    df = load_table(name, data_dir)
    if foreign_keys and not join_tables:
        print(f'ERROR: The {name} table has foreign keys {foreign_keys}'
                'but you have not passed any join_tables.')
        return None
    #rename_and_key(df, join_col, og_col='name')
    for fk in foreign_keys:
        df = join_tables(df, join_dfs[fk], fk+'_id')
    
    # Index datatable
    df = cbind(dt.Frame(id=np.arange(df.nrows) + 1), df)

    # Write to .csv
    df.to_csv(os.path.join(output_dir, f'{name}.csv'))
    return df

# from scripts.join_pset_tables import *
# tissue_df = load_join_write('tissue', data_dir, output_dir)
# join_dfs = { 'tissue': rename_and_key(tissue_df, 'tissue_id')}
# cell_df = load_join_write('cell', data_dir, output_dir, ['tissue'],join_dfs)

# deprecated
def load_pset_tables(name, file_path, psets=['CTRPv2', 'FIMM', 'gCSI', 'GDSC_v1', 'GDSC_v2', 'GRAY', 'UHNBreast', 'CCLE']):
    """
    Given the name of a table and a file_path to all PSet tables,
    read all the tables into DataFrames and return a list.

    :param: name string - Name of table
    :param: file_path string - File path to all processed data (all PSet tables)
    :return: list[DataFrame] - List of Dask DataFrames, one from each PSet
    """
    dfs = {pset: dt.fread(os.path.join(file_path, pset, name, f'{pset}_{name}.csv')) for name in tables for pset in psets}

    for pset in psets:
        # Check that this pset has this table
        pset_path = os.path.join(file_path, pset, name)
        if os.path.exists(pset_path):
            if len(os.listdir(pset_path)) == 1:
                # Written as a single CSV file with pandas read_csv
                dfs.append(dt.fread(f'{pset_path}/{pset}_{name}.csv', dtype={'symbol': 'object'}))
            else:
                # Written as many CSVs with Dask
                dfs.append(dt.fread(f'{pset_path}/{pset}_{name}-*.csv', dtype={'symbol': 'object'}))

    return dfs


# deprecated
def concat_tables(df_list, chunksize=500000):
    """
    Concatenate all DataFrames in df_list. Also drops duplicate rows
    and resets index to start at 1 and increment each row.

    @param df_list: [`list[DataFrame]`] The list of Dask DataFrames to be concatenated
    @return: [`DataFrame`] A single DataFrame containing all rows of the dfs in df_list
    """
    # Dask concat all dfs in df_list
    df = dd.concat(df_list)

    # Drop id column (old index) and any duplicate rows
    if 'id' in df.columns:
        df = df.drop('id', axis=1)
    df = df.drop_duplicates()

    # NOTE: originally concat was putting them all into the same partition
    # however now it seems to be preserving partitions so no longer necessary to repartition
    # Repartition into partitions that have chunksize rows
    #num_rows = df.shape[0].compute()  # This may take a while
    #df = df.repartition(npartitions=(num_rows // chunksize + 1))
    #df = df.repartition(npartitions=sum([df.npartitions for df in df_list]))
    #df = df.repartition(partition_size='15MB')

    return df


# deprecated
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


# deprecated
def write_table_old(df, name, file_path):
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


# deprecated
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


# deprecated
def prepare_join_table(join_df):
    """
    Transform the join_df so that its current index is stored in an 'id'
    column and its index is replaced by the 'name' column (for faster
    joins). Drop all other columns

    @param join_df: [`Dask DataFrame`] A table that has already been
                    written to disk, but which is kept in memory for
                    future joins
    @return: [`Dask DataFrame`] A copy of the table, but prepped for joins
    """
    join_df = join_df[['name']].copy()
    join_df['id'] = join_df.index
    join_df = join_df.set_index(['name'])
    return join_df


# deprecated
def safe_merge(df1, df2, fk_name, how='left'):
    """
    Given two Dask DataFrames, performs a left join to add a foreign key from
    df2 to df1. By default, the function will join on 'name' in df2 (ex. referring
    to the gene name, drug name, etc.) and on fk_name in df1.

    Once the join occurs, the forein will be renamed fk_name. safe_merge also prints a 
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
    # Reindex for better join speed (set_index is expensive in Dask)
    #df1 = df1.set_index([fk_name])

    # Join DataFrames
    join_df = dd.merge(df1, df2, left_on=fk_name, right_index=True, how=how)
    #join_df = join_df.drop(columns=fk_name)

    # Rename df2 id col to fk_name and cast as int
    join_df = join_df.rename(columns={'id': fk_name})
    #join_df = join_df[fk_name].astype(pd.Int64dfype())

    # Check if any rows did not join properly with dataset
    t1 = time()
    check_df = join_df.compute(scheduler='processes')
    t2 = time()
    # if join_df[join_df[fk_name].isna()].shape[0].compute(scheduler='processes') > 0:
    #     print(f'ERROR - some rows did not join to a {fk_name}!')

    return join_df


# deprecated
def load_join_write_old(name, fks, join_dfs, read_file_path, write_file_path):
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
    tissue_df = load_join_write('tissue', data_dir, output_dir)
    drug_df = load_join_write('drug', data_dir, output_dir)
    gene_df = load_join_write('gene', data_dir, output_dir)
    dataset_df = load_join_write('dataset', data_dir, output_dir)

    # Transform tables to be used for later joins
    tissue_df = prepare_join_table(tissue_df)
    drug_df = prepare_join_table(drug_df)
    gene_df = prepare_join_table(gene_df)
    dataset_df = prepare_join_table(dataset_df)
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
    join_dfs['cell'] = prepare_join_table(cell_df)

    # Build annotation tables
    load_join_write('drug_annotation', ['drug'], join_dfs,
                    read_file_path, write_file_path)
    load_join_write('gene_annotation', ['gene'], join_dfs, read_file_path, write_file_path)

    # Build all other secondary tables
    load_join_write('dataset_cell', ['dataset', 'cell'], join_dfs,
                        read_file_path, write_file_path)
    load_join_write('mol_cell', ['cell', 'dataset'], join_dfs,
                    read_file_path, write_file_path)
    # mol_cells has Kallisto. not sure why. from CTRPv2
    load_join_write('dataset_statistics', ['dataset'], join_dfs, read_file_path, write_file_path)
    load_join_write('gene_drug', ['gene', 'drug', 'dataset', 'tissue'], join_dfs,
                        read_file_path, write_file_path)

    return join_dfs


def build_experiment_tables(join_dfs, read_file_path, write_file_path):
    """
    Load and process experiment table, then use it to build the dose response
    and profile tables. Drop the 'name' column from the experiment table before
    writing to a CSV.

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


# TODO - better name
def build_final_tables(read_file_path, write_file_path):
    join_dfs = build_primary_tables(read_file_path, write_file_path)
    join_dfs = build_secondary_tables(join_dfs, read_file_path, write_file_path)
    build_experiment_tables(join_dfs, read_file_path, write_file_path)
    return join_dfs
