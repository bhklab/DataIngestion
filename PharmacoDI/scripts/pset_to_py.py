import glob
from collections import namedtuple
import pandas as pd
import numpy as np
from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr

pset_files = glob.glob('*/*rds')
pset_file = pset_files[1]

# TODO: Design PSet class for Python
def convert_pset_to_py(pset_file):
    """
    Read in .rds file of a PSet and parse each object slot into a nested dictionary replicating the PSet
    structure.

    :param pset_file: [string] Path to a PSet .rds file
    :return: [dict] Nested Python dictionaries holding each item of the R PSet at a key within
        it the respective dictionary.
    """
    # Allow compatible R objects to be returned as Python objects
    # - This avoids repeatedly calling pandas2ri.ri2py() for conversions
    # - Supported types include vectors, matrixes, lists and data.frames
    pandas2ri.activate()

    pgx = importr("PharmacoGx")
    readRDS = r["readRDS"]
    pset = readRDS(pset_file)

    py_pset = {
        'annotation': robj_to_dict(pset.slots['annotation']),
        'molecularProfiles': convert_molecular_prof_to_dict(pset),
        'sensitivity': robj_to_dict(pset.slots['sensitivity']),
        'perturbation': robj_to_dict(pset.slots['perturbation']),
        'drug': pset.slots['drug'],
        'datasetType': pset.slots['datasetType'],
        'cell': pset.slots['cell'],
        'curation': robj_to_dict(pset.slots['curation'])
    }

def robj_to_dict(robj):
    """
    Accepts an rpy2 R list object and converts it to a Python dictionary. All objects
    in the list are also converted to Python objects. Note the conversion may lead to

    :param robj: [R 'list'] Named R list object.
    :return: [dict] R list as a dict with all list elements coerced to Python objects.
    """
    dictionary = dict(zip(robj.names, list(robj)))
    return dictionary

def convert_molecular_prof_to_dict(pset):
    """
    Accepts a rpy2 robject containing an R PSet and coverts the molecularProfiles slot
    to a nested dictionary, storing all attributes and metadata from each
    SummarizedExperiment in its own dictionary.

    :param pset [rpy2.robject] containing an R PSet.
    :return: [dict] Pset@molecularProfiles slot as a Python dictionary, with each SummarizedExperiment
        as a dictionary, also named by slot.
    """


def convert_summarized_experiment_to_dict(sum_exp):
    """
    Convert an R SummarizedExperiment object to a Python dictionary

    :param sum_exp:
    :return:
    """
    py_sum_exp = dict(zip(sum_exp.slotnames(), list(sum_exp)))
    return py_sum_exp

def convert_dframe_to_df(dframe):
    """
    Take in an R DFrame object and convert to a Python dictionary

    :param dframe: [R 'DFrame'] R DFrame S4 class to be coerced to Python dictionary.
    :return: [DataFrame] containing the relevant information from the R DFrame.
    """
    df = pd.DataFrame(dframe.slots['listData']).T
    df.columns = list(dframe.names)

