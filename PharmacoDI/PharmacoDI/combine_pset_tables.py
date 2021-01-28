import glob
import os
import re
import numpy as np
import pandas as pd
from datatable import dt, fread, iread, join, by, rbind, cbind, f


def combine_all_pset_tables(data_dir, output_dir):
    join_dfs = combine_primary_tables(data_dir, output_dir)
    join_dfs = combine_secondary_tables(data_dir, output_dir, join_dfs)
    join_dfs = combine_experiment_tables(data_dir, output_dir, join_dfs)
    return join_dfs


def combine_primary_tables(data_dir, output_dir):
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


def combine_secondary_tables(data_dir, output_dir, join_dfs):
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
    cell_df = load_join_write(
        'cell', data_dir, output_dir, ['tissue'], join_dfs)
    join_dfs['cell'] = rename_and_key(cell_df, 'cell_id')

    # Build drug annotation table
    load_join_write('drug_annotation', data_dir,
                    output_dir, ['drug'], join_dfs, add_index=False)
    # Build gene annotation table
    gene_annot_df = load_table('gene_annotation', data_dir)
    # Remove any rows with no actual annotations (no symbol)
    gene_annot_df = gene_annot_df[dt.f.symbol > "", :]
    # Join the other way so that genes that got cut out are included back in
    gene_annot_df.key = 'gene_id'
    gene_annot_df = join_tables(join_dfs['gene'], gene_annot_df, 'gene_id')
    write_table(gene_annot_df, 'gene_annotation', output_dir, add_index=False)

    # Build all other secondary tables
    load_join_write('dataset_cell', data_dir, output_dir,
                    ['dataset', 'cell'], join_dfs)
    load_join_write('mol_cell', data_dir, output_dir,
                    ['cell', 'dataset'], join_dfs)
    # mol_cells has Kallisto. not sure why. from CTRPv2 (TODO)
    load_join_write('dataset_statistics', data_dir,
                    output_dir, ['dataset'], join_dfs)
    load_join_write('gene_drug', data_dir, output_dir, [
                    'gene', 'drug', 'dataset', 'tissue'], join_dfs)

    return join_dfs


def combine_experiment_tables(data_dir, output_dir, join_dfs):
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
    experiment_df = load_join_write('experiment', data_dir, output_dir, [
                                    'cell', 'drug', 'dataset', 'tissue'], join_dfs)
    # Don't write the 'name' column
    experiment_df[:, ['id', 'cell_id', 'drug_id', 'dataset_id', 'tissue_id']].to_csv(
        os.path.join(output_dir, 'experiment.csv'))

    # Rename columns and key experiment table based on experiment name and dataset id
    experiment_df.names = {'name': 'experiment_id'}
    experiment_df = experiment_df[:, ['id', 'experiment_id', 'dataset_id']]
    experiment_df.key = ('dataset_id', 'experiment_id')
    join_dfs['experiment'] = experiment_df

    # Nearly the same code as in load_join_write but has special case handling
    for df_name in ['dose_response', 'profile']:
        df = load_table(df_name, data_dir)
        for fk in ['dataset', 'experiment']:
            df = join_tables(df, join_dfs[fk], fk+'_id')
        del df[:, 'dataset_id']
        write_table(df, df_name, output_dir,
                    add_index=(df_name == 'dose_response'))

    return join_dfs


def load_join_write(name, data_dir, output_dir, foreign_keys=[], join_dfs=None, add_index=True):
    df = load_table(name, data_dir)
    if foreign_keys and join_dfs is None:
        raise TypeError(f'The {name} table has foreign keys {foreign_keys} '
                        'but you have not passed any join tables.')

    for fk in foreign_keys:
        df = join_tables(df, join_dfs[fk], fk+'_id')

    df = write_table(df, name, output_dir, add_index)
    return df


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
    files = [file_name for file_name in files if re.search(
        data_dir + r'/(\w+)/\1_' + name + '.csv$', file_name)]
    # Read and concatenate tables
    df = rbind(*iread(files, sep=','))
    # Replace any empty strings with None/NA
    df.replace("", None)
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
    return df  # Not necessary? df passed by reference


def join_tables(df1, df2, join_col):
    """
    Join df2 and df1 based on join_col.

    @param df1: [`datatable.Frame`] The datatable with the foreign key
    @param df2: [`datatable.Frame`] The join table (ex. tissue datatable)
    @param join_col: [`string`] The name of the columns on which the tables
                            will be joined (ex. 'tissue_id')
    @return
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


def write_table(df, name, output_dir, add_index=True):
    """
    Add a primary key to df ('id' column) and write it to output_dir
    as a .csv file.

    @param df: [`datatable.Frame`] A PharmacoDB table
    @param name: [`string`] The name of the table
    @param output_dir: [`string`] The directory to write the table to
    @return: [`datatable.Frame`] The indexed PharmacoDB table
    """
    if add_index:
        # Index datatable
        df = cbind(dt.Frame(id=np.arange(df.nrows) + 1), df)
    # Write to .csv
    df.to_csv(os.path.join(output_dir, f'{name}.csv'), quoting="none")
    return df
