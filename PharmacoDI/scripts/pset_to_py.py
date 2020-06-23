import glob
import os
from rpy2.robjects import r, pandas2ri
import pickle

if 'scripts' not in os.getcwd():
    os.chdir('scripts')

pset_files = glob.glob('../*/*rds')
pset_file = pset_files[0]


# TODO: Design PSet class for Python
def convert_pset_to_py(pset):
    """
    Read in .rds file of a PSet and parse each object slot into a nested dictionary replicating the PSet
    structure.

    :param pset_file: [string] Path to a PSet .rds file
    :return: [dict] Nested Python dictionaries holding each item of the R PSet at a key within
        it the respective dictionary.
    """
    molecular_profiles = rlist_to_dict(pset.slots['molecularProfiles'])

    py_pset = {
        'annotation': rlist_to_dict(pset.slots['annotation']),
        'molecularProfiles': {key: r_summarizedexperiment_to_dict(val) for key, val in molecular_profiles.items()},
        'sensitivity': rlist_to_dict(pset.slots['sensitivity']),
        'perturbation': rlist_to_dict(pset.slots['perturbation']),
        'drug': pset.slots['drug'],
        'datasetType': pset.slots['datasetType'],
        'cell': pset.slots['cell'],
        'curation': rlist_to_dict(pset.slots['curation'])
    }
    return py_pset


## Helper functions

def rlist_to_dict(robj):
    """
    Accepts an rpy2 R list object and converts it to a Python dictionary.

    :param robj: [R 'list'] Named R list object.
    :return: [dict] R list as a dict with all list elements coerced to Python objects.
    """
    # Handle empty list case
    if not list(robj):
        return list(robj)
    else:
        dictionary = dict(zip(robj.names, list(robj)))
        return dictionary


def try_catch(expr, error, environment):
    """
    One line function replicating the tryCatch function from R. Returns the result of `expr` if it successfully
    executes, otherwise returns `error`

    :param expr: [string] quoted version of the Python code you wish to execute
    :param error: [any] what to return if evaluation of `expr` fails
    :param environment: [dict] dictionary of variables needed to evaluate `expr`. Passed to `globals` argument of `eval`
    :return: `expr` evaluated with `environment` or `error` if the evaluation returns an error
    """
    try:
        return (eval(expr, environment))
    except:
        return (error)


def r_s4_to_dict(object, recursive=False):
    """
    Convert an R S4 object to a Python dictionary, with slot names as keys. Note this will NOT convert any additional
    S4 objects contained in the object and these will need to be converted manually.

    :param object: [R 'S4'] object to be coerced to a Python dictionary
    :return: [dict] with keys named for each slot. Converts R base types (vector, matrix, list, data.frame) to the
        appropriate Python object, but nested S4 classes must be managed manually.
    """
    # Convert S4 object to dictionary
    dct = {name: try_catch('object.slots[name]', None, {'name': name, 'object': object}) for name in
           list(object.slotnames())}

    # Convert any R list
    dct = {key: rlist_to_dict(val) if 'ListVector' in str(type(val)) else val for key, val in dct.items()}

    # Recursively convert nested S4 objects to dictionaries
    if recursive:
        dct = {key: r_s4_to_dict(val, recursive=True) if 'S4' in str(type(val)) else val for key, val in dct.items()}

    return dct


def r_summarizedexperiment_to_dict(object):
    """
    Convert an R SummarizedExperiment object into a nested Python dictionary

    :param object: A R SummarizedExperiment object or a dictionary made from one using the r_s4_to_dict function
    :return: a nested dictionary with keys corresponding to the items in each slot and values coerced to Python objects
    """
    # Convert to dictionary if necessary
    if 'dict' not in str(type(object)):
        se = r_s4_to_dict(object)
    else:
        se = object

    # Load R functions into Python
    as_data_frame = r['as.data.frame']

    # Deal with metadata
    metadata = se['metadata']
    metadata = {key: r_s4_to_dict(val, recursive=True) if 'S4' in str(type(val)) else val for key, val in
                metadata.items()}

    # Extract slot data and convert to Python
    py_se = {
        'colData': as_data_frame(se['colData']),
        ## TODO:: Determine if there is every useful information in the other slots of a SimpleList?
        'assays': r_s4_to_dict(se['assays'].slots['data'])['listData'],
        'NAMES': se['NAMES'],
        'elementMetadata': as_data_frame(se['elementMetadata']),
        'metadata': metadata
    }
    return py_se


def recursive_class(dct):
    return {key: recursive_class(val) if isinstance(val, dict) else type(val) for key, val in dct.items()}


## Script
# Allow compatible R objects to be returned as Python objects
# - This avoids repeatedly calling pandas2ri.ri2py() for conversions
# - Supported types include vectors, matrices, lists and data.frames
pandas2ri.activate()

readRDS = r["readRDS"]
pset = readRDS(pset_file)

pset_py = convert_pset_to_py(pset)
pset_py

# mprof = pset_py['molecularProfiles'].copy()
#
# with open('gCSI_mprof.pkl', 'wb') as f: pickle.dump(mprof, f, -1)
#
# # Open file for writing binary
# with open('gCSI_2017.pkl', 'wb') as outfile:
#     pickle.dump(pset_py, outfile, -1)  # -1 is shorthand for pickle.HIGHEST_PROTOCOL
