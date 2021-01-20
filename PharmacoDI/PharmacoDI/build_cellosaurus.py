import pandas as pd
import os
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import datatable as dt
from scripts.join_pset_tables import join_tables, index_and_write

if 'PharmacoDI' not in os.getcwd():
    os.chdir('PharmacoDI')

cellosaurus_path = 'data/metadata/Annotations/cellosaurus.txt'

# Using a default dict because it allows me to append duplicate indexes into a list
# Helper for build_cellosaurus_df


def build_defaultdict(tuple_list):
    def_dict = defaultdict(list)
    for tup in tuple_list:
        def_dict[tup[0]].append(tup[1])
    return def_dict


def build_cellosaurus_df(cellosaurus_path, output_dir, cell_df):
    """
    Build cellosaurus table.

    @param cellosaurus_path: [`string`] Full file path to cellosaurus file
    @param output_dir: [`string`] The directory to write the cellosaurus table
    @param cell_df: [`datatable.Frame`] The cell table; should be renamed, keyed,
                                        and shouldn't have 'tissue_id' column
    @return: [`datatable.Frame`] The cellosaurus table
    """

    with open(cellosaurus_path) as f:
        file = [line for line in f]

    file = file[55:]
    entries = ''.join(file).split('//\n')
    entry_list = [entry.split('\n') for entry in entries]
    entry_split_list = [[item.split('   ')
                         for item in entry] for entry in entry_list]
    entry_tuple_list = [[(item[0], item[1]) for item in entry if len(
        item) > 1] for entry in entry_split_list]

    pool = Pool(cpu_count() - 1)

    dict_list = pool.map(build_defaultdict, entry_tuple_list)
    dict_list = [dict(item) for item in dict_list]
    dict_list = [{key: '|||'.join(value)
                  for key, value in dct.items()} for dct in dict_list]

    cellosaurus_df = pd.DataFrame(dict_list)
    cellosaurus_df.dropna(axis=1, how='all', inplace=True)

    # Always close your pool or you will have a bunch of processes doing nothing
    pool.close()

    # Drop AG and DT columns (age of donor, date)
    cellosaurus_df.drop(columns=['AG', 'DT'], inplace=True)

    # Rename cols and add cell_id column
    rename_dict = {col: col.lower() for col in cellosaurus_df.columns}
    cellosaurus_df.rename(columns=rename_dict, inplace=True)
    cellosaurus_df.rename(
        columns={'id': 'identifier', 'ac': 'accession'}, inplace=True)
    cellosaurus_df['cell_id'] = cellosaurus_df['identifier']

    # Convert to datatable and join with cell_df
    df = join_tables(dt.Frame(cellosaurus_df), cell_df, 'cell_id')
    df = df[dt.f.cell_id >= 1, :]
    df = df[:, ['id', 'cell_id', 'identifier', 'accession', 'as', 'sy',
                'dr', 'rx', 'ww', 'cc', 'st', 'di', 'ox', 'hi', 'oi', 'sx', 'ca']]
    df = index_and_write(df, 'cellosaurus', output_dir)
    return df
