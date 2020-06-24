import pandas as pd
import wget

canonical_psets = ["CCLE_2015", "FIMM_2016", "GRAY_2017", "CTRPv2_2015", "gCSI_2017", "GDSC_2020(v1-8.2)",
                    "GDSC_2020(v2-8.2)", "GDSC_2020(v1-8.2)", "UHNBreast_2019"]

save_dir = "rawdata"

def download_psets(names, save_dir, api_url= "https://www.orcestra.ca/api/psets/available"):
    """
    Download the specified PSets from the specified api_url

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

def download_canonical_psets(save_dir, api_url= "https://www.orcestra.ca/api/psets/canonical"):
    """
    Download the specified canonical PSets from the specified api_url

    :param save_dir: [string] Path to save the PSets in
    :param api_url: [string] URL where available PSets can be retrieved. Defaults to current Orcestra API.
    :return: [None] Downloads PSet .rds files into save_dir using wget.
    """
    pset_df = pd.read_json(api_url)
    url_dict = pset_df.set_index('name').to_dict()['downloadLink']
    for name, url in url_dict.items():
        print("Downloading", name, "from", url, sep=" ")
        wget.download(url, save_dir + "/" + name + '.rds')

    return None

download_canonical_psets(save_dir)
