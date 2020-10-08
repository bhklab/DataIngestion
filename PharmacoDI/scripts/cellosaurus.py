import pandas as pd
import os
from multiprocessing import Pool, cpu_count
from collections import defaultdict

if 'PharmacoDI' not in os.getcwd():
    os.chdir('PharmacoDI')
path = 'data/metadata/Annotations/cellosaurus.txt'

with open(path) as f:
    file = [line for line in f]

file = file[55:]
entries = ''.join(file).split('//\n')
entry_list = [entry.split('\n') for entry in entries]
entry_split_list = [[item.split('   ') for item in entry] for entry in entry_list]
entry_tuple_list = [[(item[0], item[1]) for item in entry if len(item) > 1 ] for entry in entry_split_list]

# Using a default dict because it allows me to append duplicate indexes into a list
def build_defaultdict(tuple_list):
    def_dict = defaultdict(list)
    for tup in tuple_list:
        def_dict[tup[0]].append(tup[1])
    return def_dict


pool = Pool(cpu_count() - 1)

dict_list = pool.map(build_defaultdict, entry_tuple_list)
dict_list = [dict(item) for item in dict_list]
dict_list = [{key: '|||'.join(value) for key, value in dct.items()} for dct in dict_list]

cellosaurus_df = pd.DataFrame(dict_list)
cellosaurus_df.dropna(axis=1, how='all', inplace=True)

# Always close your pool or you will have a bunch of processes doing nothing
pool.close()