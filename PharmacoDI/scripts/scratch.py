import glob
import os
import re
import numpy as np
import pandas as pd
import dask.dataframe as dd

#TODO - ask about file structure (procdata/pset/tissue) or (procdata/tissue/pset)

read_file_path = os.path.join('data', 'procdata')
write_file_path = os.path.join('data') # TODO - ????

# TODO - write function to load a table
# from all psets; return in a list (?)
def load_tables(name, file_path):
    """
    Given the name of a table and a file_path to all PSet tables,
    read all the tables into DataFrames and return a list.

    :param: name string - Name of table
    :param: file_path string - File path to all processed data (all PSet tables)
    :return: list[DataFrame] - List of Dask DataFrames, one from each PSet
    """
    dfs = []
    # Get all psets with this table
    psets = glob.glob(os.path.join(file_path, name))
    # Read csv files for each pset
    for pset in psets:
        dfs.append(dd.read_csv(f'{os.path.join(file_path, name, pset)}/{pset}_{name}-*-*.csv')) #TODO - modify

    return dfs


def concat_tables(df_list):
    """
    Concatenate all DataFrames in df_list. Also drops duplicate rows
    and resets index to start at 1 and increment each row.

    :param: df_list list[DataFrame] List of Dask DataFrames to be concatenated
    :return: DataFrame A single DataFrame containing all rows of the dfs in df_list
    """
    # Dask concat all dfs in df_list
    df = dd.concat(df_list, interleave_partitions=True).compute() #TODO - read more on this to make sure you're using it correctly
    # Drop duplicate rows (inplace parameter not supported in Dask)
    df = df.drop_duplicates()
    # Reindex
    df.reset_index(drop=True)
    df.index += 1 #TODO - see if this works

    return df


# ? 
def write_df(df, file_path):
    pass


# TODO - do I need this function
def load_concat_write(name, read_file_path, write_file_path):
    df_list = load_tables(name, read_file_path)
    df = concat_tables(df_list)
    write_df(df, write_file_path)
    # TODO - how to remove sthg from memory to clear up space once you write to disk


# start script

# Primary tables can just be concatenated and written to disk
load_concat_write('tissue', read_file_path, write_file_path)
load_concat_write('drug', read_file_path, write_file_path)
load_concat_write('gene', read_file_path, write_file_path)
load_concat_write('dataset', read_file_path, write_file_path) # This seems extra

# Now the fun part starts
def load_join_table(name, file_path):
    df_files = glob.glob(os.path.join(file_path, name))
    if len(df_files) == 0:
        print('ERROR cannot find the join table')
        #raise ValueError
    return dd.read_csv(f'{os.path.join(file_path, name, name)}-*-*.csv)
    

cell_df_list = load_tables('cell', read_file_path)
cell_df = concat_tables(cell_df_list)
tissue_df = load_join_table('tissue', write_file_path)
# TODO - join ! 