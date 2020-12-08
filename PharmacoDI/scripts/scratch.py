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


def load_tables(name, file_path, psets=['CTRPv2', 'FIMM', 'gCSI', 'GDSC_v1', 'GDSC_v2', 'GRAY', 'UHNBreast', 'CCLE']):
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
    df = df.drop('id', axis=1)

    # Drop duplicate rows
    df = df.drop_duplicates()

    # Repartition into partitions, each with chunksize rows
    num_rows = df.shape[0].compute()  # This may take a while
    df = df.repartition(npartitions=(num_rows // chunksize + 1))

    return df


def write_table(df, name, file_path, chunksize=500000):
    """
    Write table to file_path
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

    df_path = os.path.join(file_path, name)
    # Check that directory exists
    if not os.path.exists(df_path):
        os.mkdir(df_path)
    else:
        # Clear all old files
        for f in os.listdir(df_path):
            os.remove(os.path.join(df_path, f))

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
    # TEMP
    if fk_name == 'experiments_id':
        fk_name = 'experiment_id'
    # Make copy of relevant cols of df2 so you don't have to drop all the other ones after merge
    df2 = df2[['id', right_on]].copy()
    # Rename df2 FK column to avoid naming conflicts (since df1 will probably have 'name' col too)
    df2 = df2.rename(columns={right_on: f'{fk_name}_{right_on}'})
    right_on = f'{fk_name}_{right_on}'

    # Join DataFrames
    # NOTE: this automatically resets the index, so it starts at 0 and doesn't have the name 'id'
    join_df = dd.merge(df1, df2, left_on=fk_name,
                       right_on=right_on, how='left')

    # Remove join columns
    join_df = join_df.drop(columns=[fk_name, right_on])
    # Rename df2 id col to fk_name
    join_df = join_df.rename(columns={'id': fk_name})

    # Check if any rows did not join properly with dataset
    # can't do this with variable? ---  join_df.query('dataset_id.isna()')
    # TODO - consider using validate param instead; consider throwing an error
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


# TODO - similar names to build_pset_tables, should I change it to be clearer?
def build_primary_tables(read_file_path, write_file_path):
    # Primary tables can just be concatenated and written to disk
    load_concat_write('tissue', read_file_path, write_file_path)
    load_concat_write('drug', read_file_path, write_file_path)
    load_concat_write('gene', read_file_path, write_file_path)
    load_concat_write('dataset', read_file_path, write_file_path)


def build_secondary_tables(read_file_path, write_file_path):
    # Other tables need to be joined to get their foreign keys
    load_join_write('cell', ['tissue'], read_file_path, write_file_path)
    load_join_write('datasets_cells', [
                    'dataset', 'cell'], read_file_path, write_file_path)
    load_join_write('drug_annotations', [
                    'drug'], read_file_path, write_file_path)
    #load_join_write('gene_annotations', ['gene'], read_file_path, write_file_path)
    # symbol col causing issues with dtype; i think because lots of vals missing
    load_join_write('experiments', [
                    'cell', 'drug', 'dataset', 'tissue'], read_file_path, write_file_path)
    load_join_write('dose_responses', [
                    'experiments'], read_file_path, write_file_path)
    load_join_write('profiles', ['experiments'],
                    read_file_path, write_file_path)

    # Drop experiment name from table
    experiment_df = load_join_table('experiments', write_file_path)
    experiment_df = experiment_df.drop(columns=['name'])
    experiment_df.set_index('id', drop=True)
    write_table(experiment_df, 'experiments', write_file_path) 
    # for whatever reason this throws a file not found error even though it's the exact same write_fxn

    load_join_write('mol_cells', ['cell', 'dataset'],
                    read_file_path, write_file_path)
    load_join_write('gene_drugs', ['gene', 'drug', 'dataset', 'tissue'], read_file_path, write_file_path)

    # Load target table and join to gene table
    target_df = load_join_table('target', read_file_path) # if this works then i should rename the fxn
    gene_df = load_join_table('gene', write_file_path)
    target_df = safe_merge(target_df, gene_df, 'gene')
    del gene_df # Remove from memory once join is complete (TODO - more of this)

    # Load drug_targets table and join to target table
    drug_targets_df = load_join_table('drug_targets', read_file_path)
    drug_df = load_join_table('drug', write_file_path)
    drug_targets_df = safe_merge(drug_targets_df, drug_df, 'drug')
    drug_targets_df = safe_merge(drug_targets_df, target_df, 'target')
    del drug_df

    write_table(target_df, 'target', write_file_path)
    write_table(drug_targets_df, 'drug_targets', write_file_path)


# mol_cells has Kallisto. not sure why. from CTRPv2

# TODO - need to load non-pset-specific tables and join those
# del or if it's inside a fxn then it's fine
# TODO - figure out singular and plural; just choose one and then fix EDR
# TODO - refactor so you're not loading join tables so many times
# TODO - load cellosaurus
