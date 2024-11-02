# Code for: "ComBat models for harmonization of resting-state EEG features in multisite studies"


Benchmarking Combining Batches (ComBat) models on spectral parameterized features of resting-state EEG signals.


- The original sources/links to raw rsEEG data are available in [our manuscript](https://doi.org/10.1016/j.clinph.2024.09.019).
  - To download LEMON dataset, we recommend the following steps:
  - First, you will need to download and unzip the LEMON behavioral and demographics file (META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv) available at the [LEMON website](https://fcp-indi.s3.amazonaws.com/data/Projects/INDI/MPI-LEMON/Compressed_tar/Behavioural_Data_MPILMBB_LEMON.tar.gz).
  - Then, run download_data_lemon.py from [Engemann D, et al. NIMG 2022](https://www.sciencedirect.com/science/article/pii/S105381192200636X?via%3Dihub), available [here](https://github.com/meeg-ml-benchmarks/brain-age-benchmark-paper).
  - Finally, convert LEMON to BIDS.

    
- Automated preprocessing was performed using [sova-harmony](https://github.com/GRUNECO/eeg_harmonization) in Python.


- Spectral parameterization was made with Fitting Oscillations and One-over Frequency [(FOOOF - now specparams)](https://fooof-tools.github.io/fooof/).

- Benchmarked ComBat models comprised:
  - [neuroCombat (standard ComBat)](https://github.com/Jfortin1/neuroCombat).
  - [neuroHarmonize (variant with nonlinear adjustment of covariates)](https://github.com/rpomponio/neuroHarmonize).
  - [OPNested-GMM (variant based on Gaussian Mixture Models to fit bimodal feature distributions)](https://github.com/hannah-horng/opnested-combat).
  - [HarmonizR (variant based on resampling to handle missing feature values)](https://github.com/SimonSchlumbohm/HarmonizR).
 
  
- The full code to replicate figures and results is implemented in Python (*combat_harmonization_eeg.ipynb*) and R (*harmonization_figures.R*)


- The *tsne* folder has the .csv files required to generate Figures 3 - 5 (look *harmonization_figures.R*).




***Qualitative visualization of batch effects***
![image](https://github.com/user-attachments/assets/6a1667b3-5838-483e-9719-5c2fe8384533)


***Statistical testing across batches***

<img src="https://github.com/user-attachments/assets/3827fceb-9cc7-48d6-9124-807b10e7836b" width="50%" />

***Downstream analysis (age-related spectral changes)***
![image](https://github.com/user-attachments/assets/ac20195e-b045-4b1b-8b0c-d08d7e266ab6)


Cite as: Jaramillo-Jimenez, A., Tovar-Rios, D. A., Mantilla-Ramos, Y.-J., Ochoa-Gomez, J.-F., Bonanni, L., & Brønnick, K. (2024). ComBat models for harmonization of resting-state EEG features in multisite studies. Clinical Neurophysiology. 2024. https://doi.org/10.1016/j.clinph.2024.09.019

