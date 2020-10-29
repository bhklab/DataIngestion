import glob
import os
import re
import numpy as np
import pandas as pd
import dask.dataframe as dd
from ..PharmacoDI import read_pset as rp

pset_name = 'GDSC_v1'
file_path = os.path.join('data', 'rawdata')
   
pset_files_df = rp.read_pset(pset_name, file_path)
pset_dict = rp.pset_df_to_nested_dict(pset_files_df)

pset_dict['molecularProfiles'].keys()
# ['Kallisto_0.46.1.isoforms', 'mutation_exome', 'cnv', 'rna', 'fusion', 
# 'Kallisto_0.46.1.rnaseq', 'Kallisto_0.46.1.isoforms.counts', 'mutation', 
# 'Kallisto_0.46.1.rnaseq.counts']

