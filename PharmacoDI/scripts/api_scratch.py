from chembl_webresource_client.new_client import new_client
import multiprocessing as mp
import numpy as np
import pandas as pd
import swifter
from collections import defaultdict
import requests
# Get a list of available tables to query
available_resources = [resource for resource in dir(new_client) if not
    resource.startswith('_')]
# Initiate connections to the tables we want to query
activity = new_client.activity
target = new_client.target
molecule = new_client.molecule
# Read in some sample genes
drug_info = pd.read_csv('data/rawdata/GDSC_v1_PSet/GDSC_v1@drug.csv.gz')
drug_ids = drug_info.inchikey.dropna()
# Setup a pool
pool = mp.Pool(mp.cpu_count())
# Retrieve all targets available in Chembl
target_result = target.filter(organism__in=[
    'Homo sapiens', 'Rattus norvegicus'])
results = list(target_result)
df = pd.DataFrame(results)
object_columns = df.dtypes[df.dtypes == 'object'].index.values
for column in object_columns:
    df = df.explode(column) # this only really explodes on cross_references
cross_references = pd.DataFrame(df.cross_references).dropna()

# what is the point of this function???
def group_dict_by_key(dictionary):
    new_dict = defaultdict(list)
    for key, value in dictionary:
        new_dict[key] = value
    return new_dict
    
# TODO:: Parse the object columns into their own DataFrames
        # @Chris idk what you mean by this --Evie
# need to make csv similar to the drugbank and save one



# This will just get a list of Uniprot IDs for all targets
# Exclude all targets with no cross-references (their target_components are NA anyways)
uniprots = pd.Series([component['accession'] for (id, component) in 
            df.query('cross_references.notna()')['target_components'].items()]).drop_duplicates()
# drop_duplicates reduces from 4949 rows to 2831!

# Map UniProt IDs to EnsemblGene IDs
# The followign code is copied directly from the Drugbank script
# So a helper needs to be made

# Split UniProt IDs into groups of 1000 entries
queries = [uniprots.iloc[i:i+1000]
           for i in np.arange(0, len(uniprots), 1000)]


def mapUniprotToEnsembl(idList):
    # Function to make GET request and process response
    
    # Make API call
    params = {
        'from': 'ID',
        'to': 'ENSEMBL_ID',
        'format': 'tab',
        'query': " ".join(idList)
    }
    r = requests.get('https://www.uniprot.org/uploadlists/', params=params)
    
    # Check that request was successful
    if r.status_code != 200:
        return f'ERROR: {r.status_code}'
    
    # Split text into rows (one row per ID) and build df
    # Exclude first row (header)
    gene_id_df = pd.DataFrame(r.text.split("\n")[1:])
    
    # Split into two columns, UniprotId and EnsemblId
    gene_id_df = gene_id_df[0].str.split(pat="\t", expand=True)
    
    # Drop empty rows and rename columns
    gene_id_df.dropna(inplace=True)
    gene_id_df.rename(columns={0:'UniProtId', 1: 'EnsemblGeneId'}, inplace=True)

    return gene_id_df


# Start pool to make API calls in parallel
pool = mp.Pool(mp.cpu_count() - 1)
uniprot_ensembls_df_list = pool.map(mapUniprotToEnsembl, queries)
pool.close()

# Combine all dfs into one and reset index
uniprot_ensembls_df = pd.concat(uniprot_ensembls_df_list)
uniprot_ensembls_df.reset_index(drop=True, inplace=True)
uniprot_ensembls_df.index = uniprot_ensembls_df.index + 1




## TODO:: Query the API to map Chembl ID to an identifier in drug_info
## - ideally it would map to IUPAC name, since we have that for every drug
## - if not we can map to smiles and inchikey then merge the non-redundant
## results

# FYI: these are the resources we want to query for targets Drugbank, Chembl, CCTRPV2, and clue.io (CMAP).
# available_resources contains a list of tables you can query
# This is the package for querying their API:
    #  https://github.com/chembl/chembl_webresource_client
# You can pip install it. Don't conda install it, the version on conda is broken.
# This is an example of the data you want to get:
    # https://www.ebi.ac.uk/chembl/g/#search_results/targets/query=erbB2
