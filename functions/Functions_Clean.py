import os
import re
import statistics
from functions.General_Functions import *
import pandas as pd
from functions.anova import Two_Way_mixed_Anova
import itertools
import numpy as np
from itertools import combinations

def Within_limits(point_to_analyse, z, q):
    if point_to_analyse > int(z) and point_to_analyse < int(q):
        return 1
    else:
        return 0


def analysed_domain_coverage(Domain_name, Domain_start, Domain_finish, Protein_entries_experiment, sequence):
    ##########
    ## This function determines the coverage of all the domain by peptides.
    ##########

    susceptibility_of_domain = 0
    total_amino_acids = []
    peptides_found = []
    domain_sequence=sequence[Domain_start:Domain_finish]
    for number3 in Protein_entries_experiment.iterrows():
        peptide_sequence = number3[1]['Peptide']
        peptide_abundance = number3[1]['spectra']

        length_of_peptide = len(peptide_sequence)
        Sequence_Positions = [m.start() for m in re.finditer(peptide_sequence, sequence)]
        for start_position in Sequence_Positions:
            # loops through each of the start positions, detected in the sequence
            value1 = start_position
            value2 = value1 + length_of_peptide
            
            # instead of doing this we just produce ranges of indexes
            domain_start_within_limits = Within_limits(value1, Domain_start, Domain_finish)
            domain_end_within_limits = Within_limits(value2, Domain_start, Domain_finish)
            
            if domain_end_within_limits == False and domain_start_within_limits == False:
                # This peptide doesnt fall within the domain ranges.
                pass
            else:
                if domain_start_within_limits == True and domain_end_within_limits == False:
                    # peptide start is in the domain, while the end expands beyond
                    #                  PEPTIDE       
                    #           DDDDDDDDDD
                    susceptibility_of_domain = susceptibility_of_domain + peptide_abundance
                    peptide_coverage = Domain_finish - value1
                    amino_acids_covered = list(range(value1, int(Domain_finish)))
                    total_amino_acids = total_amino_acids + amino_acids_covered
                    peptides_found.append(peptide_sequence)
                elif domain_start_within_limits == False and domain_end_within_limits == True:
                    # peptide start is prior to the domain start and the peptide finishes within domain limits
                    #       PEPTIDE       
                    #           DDDDDDDDDD
                    susceptibility_of_domain = susceptibility_of_domain + peptide_abundance
                    peptide_coverage = (value2 - Domain_start)
                    amino_acids_covered = list(range(int(Domain_start), value2))
                    peptides_found.append(peptide_sequence)
                    if type(amino_acids_covered) != int:
                        total_amino_acids = total_amino_acids + amino_acids_covered
                elif domain_start_within_limits == True and domain_end_within_limits == True:
                    # both start and finish is within limits.
                    #          PEPTIDE       
                    # DDDDDDDDDDDDDDDDDDDDDDDDDD
                    susceptibility_of_domain = susceptibility_of_domain + peptide_abundance
                    peptide_coverage = length_of_peptide
                    amino_acids_covered = list(range(value1, value2))
                    total_amino_acids = total_amino_acids + amino_acids_covered
                    peptides_found.append(peptide_sequence)
                elif domain_sequence in peptide_sequence:
                    # here the entire peptide has covered the domain - this was lacking in v1
                    #          PEPTIDE       
                    #           DDDDD
                    susceptibility_of_domain = susceptibility_of_domain + peptide_abundance
                    amino_acids_covered = list(range(Domain_start, Domain_finish))
                    total_amino_acids = total_amino_acids + amino_acids_covered
                    peptides_found.append(peptide_sequence)

    myset_of_covered_amino_acids = list(set(total_amino_acids))
    amino_acids_covered = len(myset_of_covered_amino_acids)
    percentage_covered = float(amino_acids_covered) / (float(Domain_finish) - float(Domain_start)) * 100
    percentage_covered = round(percentage_covered, 2)
    peptides_found = ','.join(peptides_found)
    return susceptibility_of_domain, Domain_name, percentage_covered, peptides_found


def Fasta_Analysis_Arbitarely_Domains(sequence=None, step_size=None):
    step_size = float(step_size)
    length_of_sequence = float(len(sequence))
    number_of_for_loops = int(length_of_sequence / step_size)

    df_with_doamin_info = pd.DataFrame()
    arbitary_domain_start = 0
    arbitary_domain_end = int(step_size)

    for process_number in range(0, number_of_for_loops + 1):

        if process_number == 0:
            arbitary_domain_start = arbitary_domain_start
            arbitary_domain_end = arbitary_domain_end
        else:
            arbitary_domain_start = arbitary_domain_start + step_size
            arbitary_domain_end = arbitary_domain_end + step_size
        if arbitary_domain_end > length_of_sequence:
            arbitary_domain_end = length_of_sequence
        Domain_start = arbitary_domain_start
        Domain_finish = arbitary_domain_end
        domain_name = 'Domain_p' + str(Domain_start) + '_to_' + str(Domain_finish)
        df_with_doamin_info = df_with_doamin_info.append(
            {'Name': domain_name, 'start': Domain_start, 'finish': Domain_finish}, ignore_index=True)
    df_with_doamin_info['Type'] = str(step_size) + ' AA STEP'
    return df_with_doamin_info



def produce_all_the_domains_for_protein(domains,Fasta,Domain_Types):
    
    # This function produces all the domains needed for analysis.
    if not domains.empty:
        domains.Name=domains.Name.str.replace("'", "")
    domain_ranges = [float(s.replace(" AA STEP", "")) for s in Domain_Types if "AA STEP" in s]
    for range_elem in domain_ranges:
        dom= Fasta_Analysis_Arbitarely_Domains(sequence=Fasta,step_size=range_elem)
        dom.Name = str(int(range_elem)) + "Step_" + dom.Name
        domains = pd.concat([domains, dom])
    # Here filter out the data
    if Domain_Types != None:
        domains = domains[domains.Type.isin(Domain_Types)]
    # Deal with repeated domain names
    domains.loc[domains.Name.duplicated(), "Name"] = domains[domains.Name.duplicated()].Name + "_" + domains[domains.Name.duplicated()].start.astype(int).astype(str) + "_" + domains[domains.Name.duplicated()].finish.astype(int).astype(str)
    
    # Remove any domain types that has only 1 entry as stats can not be performed on these.
    Domain_type_frequencies = domains['Type'].value_counts()
    Domain_types_not_suitable_for_stats = list(Domain_type_frequencies[Domain_type_frequencies==1].keys())
    domains = domains[~domains['Type'].isin(Domain_types_not_suitable_for_stats)]
    
    return domains
    

def MPLF_Domain_Quantifications(Protein=None, Domain_Types=None, Protein_Entries=None,Reference_Proteome=None,Reference_Domains=None):
    
    ###############
    ## This function counts the spectra within each of the defined domains for each of the counts
    ###############

    Fasta = Reference_Proteome[Reference_Proteome.Uniprot_ID==Protein]["FASTA"][0]
    domains=Reference_Domains[Reference_Domains.Uniprot_ID == Protein]
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
    for sample in value2:
        i = i + 1
        ##Get only this sample entries
        experimental_group_counts[sample] = 0.0
        # experimental_group_counts_Peptide_Counts [ sample ] = 0.0
        records = Results[Results.experiment_name == sample]
        records = records.set_index(records.Domain_Name)
        if not records.empty:
            Domain_Spectral_Count = records.NumberOfSpectra
            per_domain_peptides['Peptides ' + key] = records.peptides_found.str.replace(",", " ")
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
