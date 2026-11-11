import os
import sys
import subprocess
import copernicusmarine
import socket
from datetime import datetime
import json

dataset_id="cmems_mod_arc_phy_anfc_6km_detided_P1M-m"

output_dir = "./data/"
files_to_download_filepath = "./files_in_api.txt"
credentials_filepath = "/home/a23236/scratch/EDITO/credentials.json"

credentials = None

if "anhur" in socket.gethostname():
    output_dir = "/data/osea/scratch/oyvjak/EDITO/topaz_data"
    files_to_download_filepath = "/data/osea/scratch/oyvjak/EDITO/topaz_data/files_in_api.txt"
    credentials_filepath = "/home/a23236/scratch/EDITO/credentials.json"

def get_all_nc_filepaths(root):
    nc_filepaths = set()
    for dir_item in os.listdir(root):
        dir_item_path = os.path.join(root, dir_item)
        if os.path.isfile(dir_item_path) and dir_item_path.endswith(".nc"):
            nc_filepaths.add(dir_item_path)
        if os.path.isdir(dir_item_path):
            nc_filepaths.update(get_all_nc_filepaths(dir_item_path))
    return nc_filepaths

def get_credentials():
    global credentials
    if credentials is None:
        with open(credentials_filepath) as credentials_file:
            credentials = json.load(credentials_file)

    return credentials



def update_filepaths():
    username = get_credentials()["username"]
    password = get_credentials()["password"]
    if os.path.exists(files_to_download_filepath):
        subprocess.run(f"rm {files_to_download_filepath}", shell=True)

    query_metadata = copernicusmarine.get(dataset_id = dataset_id, dry_run=True, username=username, password=password,
                                          create_file_list=files_to_download_filepath)
    print(query_metadata)

def download_fresh_files():

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    files_in_output_dir = get_all_nc_filepaths(output_dir)
    files_in_output_dir = [os.path.basename(x) for x in files_in_output_dir]
    with open(files_to_download_filepath) as files_to_download_file:
        files_to_download = files_to_download_file.readlines()
        files_to_download = [x.strip() for x in files_to_download if len(x) > 5]
        files_to_download = [x for x in files_to_download if os.path.basename(x) not in files_in_output_dir]

    # files_to_download = list(files_to_download)[0:3]
    with open(files_to_download_filepath, "w")as files_to_download_file:
        files_to_download_file.write("\n".join(files_to_download)+"\n")


    username = get_credentials()["username"]
    password = get_credentials()["password"]
    query_metadata = copernicusmarine.get(dataset_id=dataset_id, dry_run=False , username=username, password=password,
                                          file_list=files_to_download_filepath, output_directory=output_dir)



def main():
    update_filepaths()
    download_fresh_files()


if __name__ == '__main__':
    main()