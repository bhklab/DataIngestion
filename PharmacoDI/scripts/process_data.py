import os
from PharmacoDI.read_pset import pset_df_to_nested_dict, read_pset
from build_pset_tables import build_pset_tables
from join_pset_tables import build_final_tables

pset_names = ['GDSC_v1', 'GDSC_v2', 'CTRPv2', 'FIMM', 'gCSI', 'GRAY', 'CCLE', 'UHNBreast']

pset_file_path = os.path.join('data', 'rawdata')
procdata_file_path = os.path.join('data', 'procdata')

for pset_name in pset_names:
    print(f'Building tables for {pset_name}...')
    pset_dict = pset_df_to_nested_dict(read_pset(pset_name, pset_file_path))
    build_pset_tables(pset_dict, pset_name, procdata_file_path)

#build_final_tables()

# TODO: clear memory after this