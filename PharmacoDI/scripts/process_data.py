import os
from PharmacoDI.read_pset import pset_df_to_nested_dict, read_pset
from scripts.build_pset_tables import build_pset_tables

pset_names = ['CCLE', 'GDSC_v1', 'GDSC_v2', 'CTRPv2', 'FIMM', 'gCSI', 'GRAY', 'UHNBreast']

pset_file_path = os.path.join('data', 'rawdata')
procdata_file_path = os.path.join('data', 'procdata')

for pset_name in pset_names:
    print(f'Building tables for {pset_name}...')
    pset_dict = pset_df_to_nested_dict(read_pset(pset_name, pset_file_path))
    build_pset_tables(pset_dict, pset_name, procdata_file_path)


# TODO - process all tables that are not pset specific