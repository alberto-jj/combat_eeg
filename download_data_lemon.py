# Author: Denis Engemann (https://github.com/meeg-ml-benchmarks/brain-age-benchmark-paper/blob/main/download_data_lemon.py)
# Modified by: Alberto Jaramillo-Jimenez

import os
import pathlib
import urllib.request
import pandas as pd


DEBUG = False
url_lemon = ('https://ftp.gwdg.de/pub/misc/MPI-Leipzig_Mind-Brain-Body-LEMON'
             '/EEG_MPILMBB_LEMON/EEG_Raw_BIDS_ID/')

# Download and unzip behavioural data (https://fcp-indi.s3.amazonaws.com/data/Projects/INDI/MPI-LEMON/Compressed_tar/Behavioural_Data_MPILMBB_LEMON.tar.gz)
lemon_info = pd.read_csv("D:/eeg_datasets/raw/lemon/Behavioural_Data_MPILMBB_LEMON/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv")

#As the original names in LEMON have changed, the LEMON website suggests to use these names (https://fcp-indi.s3.amazonaws.com/data/Projects/INDI/MPI-LEMON/name_match.csv)
name_match = pd.read_csv("D:/eeg_datasets/raw/lemon/name_match.csv")

# Merge the dataframes to add the Initial_ID as participant_id to lemon_info
lemon_info = lemon_info.merge(name_match, left_on='Unnamed: 0', right_on='INDI_ID', how='left')

# Rename the Initial_ID column to participant_id
lemon_info['participant_id'] = lemon_info['Initial_ID']

# Drop the columns used for matching to clean up
lemon_info = lemon_info.drop(columns=['Initial_ID', 'INDI_ID'])

data_path = pathlib.Path("D:/eeg_datasets/raw/lemon")

if not data_path.exists():
    os.makedirs(data_path)

subjects = sorted(lemon_info['participant_id'])
if DEBUG:
    subjects = subjects[:1]

extensions = ["eeg", "vhdr", "vmrk"]
good_subjects = list()

for sub in subjects:
    for ext in extensions:
        sub_url = f"{sub}/RSEEG/{sub}.{ext}"
        url = f"{url_lemon}/{sub_url}"
        out_path = data_path / sub / "RSEEG"
        if not out_path.exists():
            os.makedirs(out_path)
        out_name = out_path / f"{sub}.{ext}"
        try:
            urllib.request.urlretrieve(url, out_name)
            good_subjects.append(sub)
        except Exception as err:
            print(err)

good_subs_df = pd.DataFrame(dict(subject=list(set(good_subjects))))
good_subs_df.to_csv('D:/eeg_datasets/raw/lemon/good_subjects.csv')
