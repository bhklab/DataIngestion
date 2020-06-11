import pandas as pd
import numpy as np
import wget
import glob
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr

canonical_psets = ["CCLE_2015", "FIMM_2016", "GRAY_2017", "CTRPv2_2015", "gCSI_2017", "GDSC_2020(v1-8.2)",
                    "GDSC_2020(v2-8.2)", "GDSC_2020(v1-8.2)", "UHNBreast_2019"]

save_dir = "rawdata"

def download_psets(names, save_dir, api_url= "https://www.orcestra.ca/api/psets/available"):
    """
    Download the specified canonical PSets from the specified api_url

    :param names: [list] Names of PSets to download. Must match the 'names' in the api call JSON.
    :param save_dir: [string] Path to save the PSets in
    :param api_url: [string] URL where available PSets can be retrieved. Defaults to current Orcestra API.
    :return: [None] Downloads PSet .rds files into save_dir using wget.
    """
    pset_df = pd.read_json(api_url)
    pset_df = pset_df[pset_df.name.isin(names)]
    url_dict = pset_df.set_index('name').to_dict()['downloadLink']
    for name, url in url_dict.items():
        print("Downloading", name, "from", url, sep=" ")
        wget.download(url, save_dir + "/" + name + '.rds')

    return None


download_psets(canonical_psets, save_dir)

pset_files = glob.glob('*/*rds')

def pickle_pset(pset_file, pset_name, save_dir):
    """
    Read in .rds files of the PSet and parse each object slot into a nested dictionary replicating the PSet structure

    :param pset_file: [string] Path to a PSet .rds file
    :param save_dir: [string] Path to save the pickled PSet dictionary.
    :param pset_name: [string] Name of the PSet being pickled. Will be used to name the saved .pkl file.
    :return: [None] Save
    """
    r = ro.r
    pgx = importr("PharmacoGx")
    dt = importr("data.table")
    readRDS = r["readRDS"]
    pset = readRDS(pset_file)


def convert_pset_sensitivity_to_dict(pset):
    pgx = importr("PharmacoGx")
    sensitivity_info = pgx.sensitivityInfo(pset)
    sensitivty_raw = pgx.sensitivityRaw(pset)
    sensitivity_profiles = pgx.sensitivtyProfiles(pset)

def convert_pset_drug_to_df():
    pgx = importr("PharmacoGx")

def convert_pset_drug():
    pgx = importr("PharmacoGx")


def convert_pset_drug():
    pgx = importr("PharmacoGx")


def convert_pset_drug():
    pgx = importr("PharmacoGx")






