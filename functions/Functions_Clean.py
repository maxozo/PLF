import os
import re
import statistics
from .General_Functions import *
import pandas as pd
from .anova import Two_Way_mixed_Anova
import itertools
import numpy as np
from itertools import combinations



def peptide_start_within_domain(pep_start, domain_start, domain_end):
    return domain_start <= pep_start < domain_end       # checks if the pep start pos is within domain limits, including domain start pos, excluding domain end pos. Will return T if within limits

def peptide_end_within_domain(pep_end, domain_start, domain_end):
    return domain_start < pep_end <= domain_end     # checks if the pep end pos is within domain limits, excluding domain start pos, including domain end pos. Will return T if within limits

##### these, plus some changes to the func below, ensure that:
        # Peptides starting exactly at the domain start are included.
        # Peptides ending exactly at the domain start are excluded.             (these cases were not considered in og PLF)
        # Peptides ending exactly at the domain end are included.
        # Peptides starting exactly at the domain end are excluded.
        # Peptides which both start and end exactly at domain limits are included (i.e. peptides which are the exact span of the domain)
    # and all rules from before:
        # Peptides completely outside of the domain (no overlap) are excluded.
        # Peptides partially in the domain are included
        # Peptides completely within the domain are included
        # Peptides spanning across both sides of the domain are included. (this was attempted but not functional in og PLF)



def analysed_domain_coverage(Domain_name, Domain_start, Domain_finish, Protein_entries_experiment, sequence):

    susceptibility_of_domain = 0.0
    total_amino_acids = set()           # making it a set to start off with, saves converting twice at the end
    peptides_found = []

    for row in Protein_entries_experiment.itertuples(index=False):      # itertuples is faster than iterrows
        peptide_sequence = row.Peptide
        peptide_abundance = row.spectra
        peptide_length = len(peptide_sequence)

        sequence_positions = [m.start() for m in re.finditer(re.escape(peptide_sequence), sequence)]     # adding the re.escape() ensures that any special regex characters in the peptide sequence are escaped and treated as characters (recommended for AA sequence data)

        for start_pos in sequence_positions:        # save a step
            end_pos = start_pos + peptide_length

            start_within = peptide_start_within_domain(start_pos, Domain_start, Domain_finish)
            end_within = peptide_end_within_domain(end_pos, Domain_start, Domain_finish)

            if not (start_within or end_within):       # shorter way of saying if both F
                if start_pos <= Domain_start and end_pos >= Domain_finish:                # if peptide starts before or at domain start, and ends after or at domain end (better than above trying to match strings, and avoids slight mismatches in pep sequence hindering matching)
                    susceptibility_of_domain += peptide_abundance
                    total_amino_acids.update(range(Domain_start, Domain_finish))          # bc its a set, just using .update() and adding right away
                    peptides_found.append(peptide_sequence)
                continue

            if start_within and not end_within:  # Peptide starts in domain, ends outside
                susceptibility_of_domain += peptide_abundance
                total_amino_acids.update(range(start_pos, Domain_finish))
                peptides_found.append(peptide_sequence)

            elif not start_within and end_within:  # Peptide starts before domain, ends in domain
                susceptibility_of_domain += peptide_abundance
                total_amino_acids.update(range(Domain_start, end_pos))
                peptides_found.append(peptide_sequence)

            elif start_within and end_within:  # Entirely within domain
                susceptibility_of_domain += peptide_abundance
                total_amino_acids.update(range(start_pos, end_pos))
                peptides_found.append(peptide_sequence)

    amino_acids_covered = len(total_amino_acids)
    domain_length = Domain_finish - Domain_start
    percentage_covered = round((amino_acids_covered / domain_length) * 100, 2)
    peptides_found_str = ','.join(peptides_found)

    return susceptibility_of_domain, Domain_name, percentage_covered, peptides_found_str



def Fasta_Analysis_Arbitrarily_Domains(sequence, step_size):

    step_size = int(step_size)
    length = len(sequence)

    domain_info_list = [        # list comprehension
        {
            'Name': f'Domain_p{start}_to_{min(start + step_size, length)}',     # dynamic domain name - get start pos (from the for loop, below), then calculate end pos by adding step size to the start (e.g. 100 + 50AA step = 150),
                                                                                         # but wrapping min(.... ,length) ensures that if the calculated end pos is >length, it will put end pos as the length value
            'start': start,
            'finish': min(start + step_size, length)    # end pos calculated as above
        }
        for start in range(0, length, step_size)        # for every step_size increment along the full length of the protein, starting at 0
    ]

    df = pd.DataFrame(domain_info_list)
    df['Type'] = f'{step_size} AA STEP'
    return df



def produce_all_the_domains_for_protein(domains, fasta, domain_types):
    if not domains.empty:
        domains['Name'] = domains['Name'].str.replace("'", "", regex=False)

    domain_ranges = [
        float(s.removesuffix(" AA STEP")) 
        for s in domain_types or []         # (added empty list here to loop through if domain_types is empty - avoids errors)
        if "AA STEP" in s
    ]

    if domain_ranges:       # if there are numbers in this ^ list:
        new_domains = []
        for range_elem in domain_ranges:
            dom = Fasta_Analysis_Arbitrarily_Domains(sequence=fasta, step_size=range_elem)  
            dom['Name'] = f"{int(range_elem)}Step_" + dom['Name'].astype(str)
            new_domains.append(dom)     # add to list (more efficient than adding to pd df)
        if new_domains:     # if something was added
            domains = pd.concat([domains, *new_domains], ignore_index=True)     # add them to the biological domains df (input) (using * on new_domains list makes it into individual dfs to concat)

    if domain_types:
        domains = domains[domains['Type'].isin(domain_types)].copy()

    duplicated = domains['Name'].duplicated(keep=False)
    if duplicated.any():
        domains.loc[duplicated, 'Name'] += (
            "_" + domains.loc[duplicated, 'start'].astype(int).astype(str)
            + "_" + domains.loc[duplicated, 'finish'].astype(int).astype(str)
        )

    counts = domains['Type'].value_counts()
    unsuitable = counts[counts == 1].index
    domains = domains[~domains['Type'].isin(unsuitable)].copy()
    return domains
    

def MPLF_Domain_Quantifications(Protein=None, Domain_Types=None, Protein_Entries=None,Reference_Proteome=None,Reference_Domains=None):
    
    ###############
    ## This function counts the spectra within each of the defined domains for each of the counts
    ###############
    Fasta = Reference_Proteome[Reference_Proteome.Uniprot_ID==Protein]["FASTA"]
    if len(Fasta)==0:
        Fasta = Reference_Proteome[Reference_Proteome.Uniprot_AC.str.contains(f"^{Protein}|\\n{Protein}",regex=True)]["FASTA"]
    Fasta = Fasta.values[0]

    domains = Reference_Domains[Reference_Domains.Uniprot_ID==Protein]
    if len(domains)==0:
        domains = Reference_Domains[Reference_Domains.Uniprot_AC==Protein]
    # if len(domains)==0:
    #     domains = Reference_Domains[Reference_Domains.Uniprot_AC.str.contains(f"^{Protein}|\\n{Protein}",regex=True)]["FASTA"]
    # domains = Fasta.values[0]
    
    # domains=Reference_Domains[Reference_Domains.Uniprot_ID == Protein]
    # Prot1=Reference_Proteome[Reference_Proteome["Uniprot_AC"].str.contains(f"^{Protein}|\\n{Protein}",regex=True)]
    
    All_Protein_Domains_suitable_for_stats = produce_all_the_domains_for_protein(domains,Fasta,Domain_Types)
    
    experiment_names = Protein_Entries['Sample'].unique()
    Experiment_dict = {}
    Experiment_Coverages = pd.DataFrame()
    for experiment in experiment_names:
        Protein_entries_experiment = Protein_Entries[Protein_Entries['Sample'] == experiment]
        Exclusive_spectrum_count = Protein_entries_experiment['spectra'].sum()
        for index, row in All_Protein_Domains_suitable_for_stats.iterrows():
            Domain_name = row.Name
            Domain_start = row.start
            Domain_finish = row.finish
            Domain_type = row.Type

            susceptibility_of_domain, Domain_name, percentage_covered, peptides_found = analysed_domain_coverage(
                Domain_name, int(Domain_start), int(Domain_finish), Protein_entries_experiment, Fasta)

            Experiment_Coverages = Experiment_Coverages.append(
                {'Domain_Name': Domain_name,  #
                    'Domain_Start': Domain_start,  #
                    'Domain_Finish': Domain_finish,  #
                    'Domain Type': Domain_type,  #
                    'NumberOfSpectra': susceptibility_of_domain,  #
                    'Percent_Covered': percentage_covered,  #
                    'Exclusive_spectrum_count': Exclusive_spectrum_count,  #
                    'experiment_name': experiment,  #
                    'GeneAC': Protein,  #
                    'peptides_found': peptides_found  #
                    }, ignore_index=True)
            
        Experiment_dict[experiment] = Experiment_Coverages
    return Experiment_Coverages


def per_domain_quantification_matrix(key, value2, Results):
    
    ##############
    ## This function emits 2 matrises - 
    # 2) Matrix with the per sample counts 
    # 3) Peptides identified for each of the domains in the experimental group
    ##############
    
    experimental_group_counts = Results[
        ['Domain_Name', 'Domain_Start', 'Domain_Finish', 'Domain Type']].drop_duplicates()
    experimental_group_counts = experimental_group_counts.set_index(experimental_group_counts.Domain_Name)
    experimental_group_counts[
        'Domain_Length'] = experimental_group_counts.Domain_Finish - experimental_group_counts.Domain_Start
    col = 'Peptides ' + key
    per_domain_peptides = pd.DataFrame([''] * experimental_group_counts.__len__(), columns=[col]).set_index(
        experimental_group_counts.index)
    i = -1
    per_domain_peptides['Peptides ' + key] = ''
    # per_domain_peptides['Peptides c1' + key] = 0
   
    for sample in value2:
        i = i + 1
        ##Get only this sample entries
        experimental_group_counts[sample] = 0.0
        # experimental_group_counts_Peptide_Counts [ sample ] = 0.0
        records = Results[Results.experiment_name == sample]
        records = records.set_index(records.Domain_Name)
        if not records.empty:
            Domain_Spectral_Count = records.NumberOfSpectra
            # total = 0
            for id1 in records.peptides_found.index:
                # print(id1)
                l1 = records.peptides_found[id1].split(',')
                l1 = [element for element in l1 if element != '']
                
                # total+=len(l1)
                # print(total)
                # print(l1)
                s1 =set(records.peptides_found[id1].split(','))
                s1.discard('')
                # per_domain_peptides.loc[id1,'Peptides c1' + key] = per_domain_peptides.loc[id1,'Peptides c1' + key] +len(l1)
                if len(s1)>0:
                    per_domain_peptides.loc[id1,'Peptides ' + key] = per_domain_peptides.loc[id1,'Peptides ' + key]+' '+' '.join(s1)
            experimental_group_counts[sample] = Domain_Spectral_Count
    return experimental_group_counts, per_domain_peptides


def normalise(Spectra_Matrix, Spectral_total_counts):
    # This function scales the counts based on the normalisation factors.
    Spectra_Matrix_output = Spectra_Matrix.copy()
    Spectra_Values = Spectra_Matrix.iloc[:, 5:]
    for column in Spectra_Values:        # print(column)
        if column in Spectral_total_counts.index:

            Factor = float(Spectral_total_counts.loc[column,"Factor"])
            value_for_Analysis = round(Spectra_Matrix[column] * Factor)
            Spectra_Matrix_output[column] = value_for_Analysis
        else:
            print("There is no sample present in the: "+column)
            continue
    return Spectra_Matrix_output


def determine_differences(Norm_Stats_Spectra, Domain_lengths):
    df = Norm_Stats_Spectra
    uniques = [*df.keys()]
    # todo have to change the sorting function. curently we pick the column names randomly which may cause issues.
    differences = pd.DataFrame()
    averages = {}
    for comb in combinations(uniques, 2):
        Average_Spectra1 = round(Norm_Stats_Spectra[comb[0]].iloc[:, 5:].mean(axis=1))

        Average_Spectra2 = round(Norm_Stats_Spectra[comb[1]].iloc[:, 5:].mean(axis=1))

        averages[comb[0]] = Average_Spectra1
        averages[comb[1]] = Average_Spectra2
        Difference = (Average_Spectra1 - Average_Spectra2) / Domain_lengths
        differences['Diff ' + comb[1] + ' vs ' + comb[0]] = Difference
    return differences, averages


def drop_one_of_the_overlapping_domains(domains, Results):       
    Combinations = list(itertools.combinations(domains.Domain_Name, 2))
    to_drop = []
    for combo in Combinations:
        y = range(int(domains[domains.Domain_Name == combo[0]].Domain_Start),
                  int(domains[domains.Domain_Name == combo[0]].Domain_Finish))
        x = range(int(domains[domains.Domain_Name == combo[1]].Domain_Start),
                  int(domains[domains.Domain_Name == combo[1]].Domain_Finish))
        xs = set(x)
        overlap = xs.intersection(y)
        if len(overlap) > 0:
            ##Here we pick only one of the 2 domains
            if len(x) > len(y):
                to_drop.append(combo[0])
            else:
                to_drop.append(combo[1])
        else:
            continue
    # Here we drop any overlapping domains
    domains = domains[~domains.Domain_Name.isin(to_drop)]
    Results = Results[Results.Domain_Name.isin(domains.Domain_Name)]
    return Results


def keys_to_Pandas(experiment_feed):
    Data = []
    for key, value2 in experiment_feed.items():
        try:
            Data = pd.concat([Data, value2], axis=1)
        except:
            Data = value2
    Data = Data.loc[:, ~Data.columns.duplicated()]
    return Data

def Add_values_to_missing_experiments(Results,experiment_feed):
    
    ##################
    # There may be samples determined in the experimental design but not quantified in the PLF
    # For these we are adding 0 to the normalisation
    ###################
    value=[]
    Spectral_total_counts = Results[["experiment_name", "Exclusive_spectrum_count"]].drop_duplicates()

    for key,value2 in experiment_feed.items():
        value.extend(value2)
    for value1 in value:
        if Spectral_total_counts[(value1 == Spectral_total_counts.experiment_name)].experiment_name.count()>0:
            if (float(Spectral_total_counts[(value1 == Spectral_total_counts.experiment_name)].Exclusive_spectrum_count)==0.0):
                Spectral_total_counts[(value1 == Spectral_total_counts.experiment_name)].Exclusive_spectrum_count=0.0000001
            else:
                continue
        else:
            Spectral_total_counts=Spectral_total_counts.append({"experiment_name": value1, "Exclusive_spectrum_count": 0.0000001}, ignore_index=True)
    # Here have to check whether all the samples are in the file, if not we add 0 as a norm factor.
    Spectral_total_counts = Spectral_total_counts.set_index("experiment_name")
    Median_Norm_Factor = statistics.median(
        [Spectral_total_counts["Exclusive_spectrum_count"].max(),
        Spectral_total_counts["Exclusive_spectrum_count"].min()])
    
    Spectral_total_counts["Median_Norm"] = Median_Norm_Factor
    Spectral_total_counts["Factor"] = Spectral_total_counts["Median_Norm"] / Spectral_total_counts[
        "Exclusive_spectrum_count"]
    return Spectral_total_counts

def MPLF_Statistical_Analyisis(experiment_feed=None, Results=None, Protein=None,paired=True,cuttoff_p_val=0.05):
    Combined_Data_Frame_With_Statistics = pd.DataFrame()
    Unique_Domains = Results['Domain Type'].unique()
    Datafile = Results

    
    Spectral_total_counts = Add_values_to_missing_experiments(Results,experiment_feed)
    p_val_name = None
    # Loop through each of the domain types
    for Data_Type in Unique_Domains:
        Results = Datafile[Datafile['Domain Type'] == Data_Type]
        domains = Results[["Domain_Name", "Domain_Finish", "Domain_Start"]]
        domains = domains.drop_duplicates()
        # Here Check for the overlap of the domains and if present then just pick one that covers the largest area
        Results = drop_one_of_the_overlapping_domains(domains, Results)
        if not Results.empty:
            experiments_all_Spectra = {}

            experiments_all_Peptides = {}
            experiments_all_Spectra_For_Stats = {}

            # Get the spectral total counts for different samples.
            All_initial_peptide_quantifications =pd.DataFrame()
            for key, value2 in experiment_feed.items():

                Spectra_Matrix, Peptides = per_domain_quantification_matrix(key, value2,
                                                                Results)
                Norm_Stats_Spectra = normalise(Spectra_Matrix,Spectral_total_counts)
                experiments_all_Spectra[key] = Norm_Stats_Spectra
                experiments_all_Spectra_For_Stats[key] = Norm_Stats_Spectra.drop(
                    columns=["Domain_Start", "Domain_Finish", "Domain Type", "Domain_Length"])
                Original = Spectra_Matrix.iloc[:,5:].add_prefix("Original_")
                All_initial_peptide_quantifications=pd.concat([All_initial_peptide_quantifications,Original], axis=1, sort=False)
                experiments_all_Peptides[key] = Peptides

            Domain_lengths = Norm_Stats_Spectra.Domain_Length
            Differences, Averages = determine_differences(experiments_all_Spectra, Domain_lengths)
            
            Peptides_PD = keys_to_Pandas(experiments_all_Peptides)
            All_Spectra = keys_to_Pandas(experiments_all_Spectra).drop_duplicates()
            Averages = pd.DataFrame(Averages)
            
            p_values_adjusted = Two_Way_mixed_Anova(experiments_all_Spectra_For_Stats,paired=paired)
            p_values_adjusted = p_values_adjusted.set_index('Domain_Name')
            p_val_name = p_values_adjusted.columns.values[0]
            Combined_Data_Frame_for_this_experiment = pd.concat([All_Spectra, Differences, p_values_adjusted, Averages, Peptides_PD,All_initial_peptide_quantifications], axis=1)
            Combined_Data_Frame_With_Statistics = Combined_Data_Frame_With_Statistics.append(Combined_Data_Frame_for_this_experiment)
        else:
            continue
    
    if Combined_Data_Frame_With_Statistics.__len__()>0:
        col = Combined_Data_Frame_With_Statistics.pop(p_val_name)
        Combined_Data_Frame_With_Statistics.insert(3, col.name, col)
        if (Combined_Data_Frame_With_Statistics[p_val_name] <= cuttoff_p_val).any():
            Combined_Data_Frame_With_Statistics['GeneAC'] = Protein
            Combined_Data_Frame_With_Statistics = Combined_Data_Frame_With_Statistics.loc[:, ~Combined_Data_Frame_With_Statistics.columns.duplicated()]
        else:
            return None, None
    return Combined_Data_Frame_With_Statistics, Spectral_total_counts


def Master_Run_Score_Calculations(Structural_Analysis_Results, Protein):
    Scores = {}
    Results_Step = Structural_Analysis_Results[Structural_Analysis_Results['Domain Type'].str.contains('STEP')]
    each_p_val = [x for x in Structural_Analysis_Results.keys().values if x.startswith('p:')]
    for p_elem in each_p_val:
        keys_h = p_elem.replace("p: ", "").split(" vs ")
        Diff_Exp_name = [x for x in Structural_Analysis_Results.keys().values if
                         ((keys_h[0]) in x and (keys_h[1]) in x and 'Diff ' in x)][0]
        Diff = Results_Step[Diff_Exp_name][Results_Step[p_elem] < 0.05]
        Trend = Diff.mean()
        Score = abs(Diff).mean()
        # if math.isnan(Score):
        #     Score=0
        # if math.isnan(Trend):
        #     Trend=0
        Scores['Score: ' + keys_h[0] + ' vs ' + keys_h[1]] = Score
        Scores['Trend: ' + keys_h[0] + ' vs ' + keys_h[1]] = Trend

    return pd.DataFrame(Scores, index=[Protein])
