# This repository contains PLF code

# ![Protein Locational Fingerprinter](assets/images/PLF.png) 
The code is utilised as a backend f MPLF code to analyse the structural datasets

This branch allows users to run MPLF locally. We have expanded the analysis by integrating the MOUSE, HORSE and RABBIT proteomes besides the HUMAN.

Furthermore the analysis can now be run on the HPC clusters.

## Quick start
1. Prepeare a file that lists **Protein** name (optional if source protein not determined), **Peptide** sequence (remove any special characters from these), **Sample** of protein belionging and **spectra** (can be multiple columns as per: **spectra_1**,**spectra_2**, etc. -- these will be added up): as per [this file](https://github.com/maxozo/MPLF/blob/mplf_package/Sample_Data/sample_inputs_small/Sample_Data_For_Analysis.csv).
2. Prepeare a tsv file that lists the experimental design - as per [this file](https://github.com/maxozo/MPLF/blob/mplf_package/Sample_Data/sample_inputs_small/Experiment_feed.tsv). If paired make sure that the rows list matching pairs, otherwise any order is ok.
3. Run the MPLF pipeline:
    `python MPLF.py --experimental_design Experiment_feed.tsv --peptides Sample_Data_For_Analysis.csv --spiecies HUMAN`