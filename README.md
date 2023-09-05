# This repository contains PLF code

# ![Protein Locational Fingerprinter](assets/images/PLF.png) 
**Peptide location fingerprinting (PLF)** is a technique capable of identifying modified proteins and potential causal mechanisms in complex biological samples. In standard proteomics, proteins are trypsinised which generate peptides whose sequence identities and relative abundances are measured by LC-MS/MS. During this process most proteins are only partially digested, due to differing solubilities, stabilities and enzyme susceptibilities related to their higher order structures. By mapping and quantifying LC-MS/MS-detected peptides within specific regions, PLF enables the detection of statistical differences in the regional digestibility along the protein structure due to ageing and disease mechanisms.

## Quick start
1. Prepeare a file that lists **Protein** name (optional if source protein not determined), **Peptide** sequence (remove any special characters from these), **Sample** of protein belionging and **spectra** (can be multiple columns as per: **spectra_1**,**spectra_2**, etc. -- these will be added up): as per [this file](https://github.com/maxozo/MPLF/blob/mplf_package/Sample_Data/sample_inputs_small/Sample_Data_For_Analysis.csv).
2. Prepeare a tsv file that lists the experimental design - as per [this file](https://github.com/maxozo/MPLF/blob/mplf_package/Sample_Data/sample_inputs_small/Experiment_feed.tsv). If paired make sure that the rows list matching pairs, otherwise any order is ok.
3. Run the MPLF pipeline:
    `python PLF.py --experimental_design Experiment_feed.tsv --peptides Sample_Data_For_Analysis.csv --spiecies HUMAN`