import glob
import os
import re
import numpy as np
import pandas as pd
import swifter  # Library to parallelize apply statements automagically

pset_name = 'GDSC_v1'
file_path = os.path.join('..', 'data', 'rawdata')
slot_names=['curation', 'drug', 'molecularProfiles', 'sensitivity', 'annotation', 'cell']

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
        raise ValueError(f'No PSet directory named {pset_name} could be found in {file_path}')
    
    slot_patterns = [f'.+@{slot}.+[^rds]$' for slot in slot_names]
    pset_files = pd.Series(os.listdir(pset_dir))
    pset_files_df = pset_files.str.split('$', expand=True)
    pset_files_df['file_paths'] = [os.path.join(pset_dir, file_name) for file_name in pset_files]
    pset_files_df.columns = ['slot', *[f'subitems{i}' for i in range(1, pset_files_df.shape[1] - 1)], 'file_paths']
    pset_files_df['data'] = pset_files_df['file_paths'].swifter.apply(read_pset_file)

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
        raise ValueError(f'Unsupported file type passed to this function from: {file_path}')