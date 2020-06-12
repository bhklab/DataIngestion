import glob
import pandas as pd
from rpy2 import robjects
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
    # - Supported types include vectors, matrices, lists and data.frames
    pandas2ri.activate()

    pgx = importr("PharmacoGx")
    readRDS = r["readRDS"]
    pset = readRDS(pset_file)

    py_pset = {
        'annotation': rlist_to_dict(pset.slots['annotation']),
        'molecularProfiles': convert_molecular_prof_to_dict(pset),
        'sensitivity': rlist_to_dict(pset.slots['sensitivity']),
        'perturbation': rlist_to_dict(pset.slots['perturbation']),
        'drug': pset.slots['drug'],
        'datasetType': pset.slots['datasetType'],
        'cell': pset.slots['cell'],
        'curation': rlist_to_dict(pset.slots['curation'])
    }
    return py_pset


def rlist_to_dict(robj):
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
    mol_prof = pset.slots['molecularProfiles']
    mprof = dict(zip(mol_prof.names.tolist(), list(mol_prof)))


def r_summarized_experiment_to_dict(sum_exp):
    """
    Convert an R SummarizedExperiment object to a Python dictionary

    :param sum_exp: [R 'SummarizedExperiment'] to convert to a Python dictionary.
    :return: [dict] with keys named after the SummarizedExperiment slots
    """
    # @elementMetadata (rowData)
    dframe = sum_exp.slots['elementMetaData']
    rowdata_df = r_dframe_to_df(dframe)

    # @metadata
    metadata_dict = convert_se_metadata_to_dict(sum_exp.slots['metadata'])

    # @colData

    # @assays

    # @NAMES


def r_dframe_to_df(dframe):
    """
    Take in an R DFrame object and convert to a Pandas DataFrame

    :param dframe: [R 'DFrame'] R DFrame S4 class to be coerced to Python dictionary.
    :return: [DataFrame] containing the relevant information from the R DFrame object.
    """
    ##TODO: Determine if dframe.slots['metadata'] ever has anything useful it in?
    df = pd.DataFrame(dframe.slots['listData']).T
    df.columns = list(dframe.names)
    df.set_index('rownames')
    return df


def convert_se_metadata_to_dict(metadata):
    """
    Convert an R SummarizedExperiment@metadata slot into a dictionary, coercing subclasses to the appropriate type

    :param metadata: [R 'list] Containing metadata for the SummarizedExperiment, this may include other Bioconductor
        S4 classes.
    :return: [dict] A nested dictionary with keys corresponding to the names of slots, items or columns from the
        respective R object.
    """
    metadata = rlist_to_dict(metadata)

    # $experimentData
    expr_data = r_expr_data_to_df(metadata)

    # $annotation
    annotation = metadata['annotation'][0]

    # $protocolData
    prot_data = r_annotated_df_to_df(metadata)


def r_expr_data_to_df(expr_data):
    """
    Convert an R ExperimentData S4 object to a Pandas DataFrame.

    :param expr_data: [R 'ExperimentData'] to coerce to a DataFrame.
    :return: [DataFrame] containing t
    """


def r_annotated_df_to_df(annot_df):
    """
    Convert an R AnnotatedDataFrame to a Pandas DataFrame.

    :param annot_df: [R 'AnnotatedDataFrame'] to coerce to a DataFrame.
    :return: [DataFrame] containing the data from the AnnotatedDataFrame.
    """


def r_simple_assay_to_ndarray(assay):
    """
    Convert an R SimpleAssay S4 object to a Numpy ndarray.

    :param assay: [R 'SimpleAssay'] to coerce to an ndarray
    :return: [ndarray] with the same values and dimensions as the SimpleAssay
    """


def r_s4_to_dict(object):
    """
    Convert an R S4 object to a Python dictionary, with slot names as keys. Note this will NOT convert any additional
    S4 objects contained in the object and these will need to be converted manually.

    :param object: [R 'S4'] object to be coerced to a Python dictionary
    :return: [dict] with keys named for each slot. Converts R base types (vector, matrix, list, data.frame) to the
        appropriate Python object, but nested S4 classes must be managed manually.
    """
    ListVector = robjects.vectors.ListVector

    dct = {name: object.slots[name] for name in list(object.slotnames())}
    #dct = {key: (rlist_to_dict(value) if isinstance(value, ListVector) else value) for key, value in dct.items()}

    return dct


def r_s4_to_py(object):
    """
    Convert an R S4 object to a dictionary, with slot names as keys. If nested S4 objects are found,
    they are converted recursively to dictionaries until no more R S4 objects are contained in the
    dictionary tree.

    :param object: [R 'S4] object to coerce to Python objects.
    :return: [dict] of Python objects or additional dictionaries if there are nested S4 classes.
    """
    dct = r_s4_to_dict(object)

    classes = ['S4' in str(type(value)) for value in dct.values()]
    print(classes)

    if True in classes:
        return {key: (r_s4_to_py(value) if 'S4' in str(type(value)) else value) for key, value in dct.items()}
    else:
        return dct
