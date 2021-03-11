import os
PASSWORD = os.getenv('MYSQL_PASSWORD')
USER = os.getenv('MYSQL_USER')
HOST = os.getenv('MYSQL_HOST')
import re
import statistics
from General_Functions import *
import pandas as pd

from anova import Two_Way_mixed_Anova
import itertools



def retrieve_reviewed(Protein_list):
    ##In this function we should use MySQL
    mydb = mysql.connector.connect(
        host=HOST,
        user=USER,
        passwd=PASSWORD,
        database="Uniprot",
        auth_plugin='mysql_native_password'
    )
    mycursor = mydb.cursor()
    sql = "SELECT `Uniprot_ID` FROM `Uniprot_Analysed`"
    # Reviewed_proteins = (pd.read_csv('Structural/Structural_Algorythm/data/Reviewed', index_col=False))
    mycursor.execute(sql)
    Reviewed_proteins = pd.DataFrame(mycursor.fetchall())
    if not Reviewed_proteins.empty:
        field_names = [i[0] for i in mycursor.description]
        Reviewed_proteins.columns = field_names
    mydb.disconnect()
    mycursor.close()
    Protein = ""
    for protein1 in Protein_list:
        if protein1 in Reviewed_proteins['Uniprot_ID'].values.tolist():
            return protein1
    if Protein == "":
        return Protein_list[0]


def retrieve_FASTA(Gene_Name_global):
    # TODO -- if we plan to expand to different spiecies here we need to process the other spiecies uniprot database.
    mydb = mysql.connector.connect(
        host=HOST,
        user=USER,
        passwd=PASSWORD,
        database="Uniprot",
        auth_plugin='mysql_native_password'
    )
    mycursor = mydb.cursor()
    try:
        mycursor.execute("SELECT FASTA FROM `Uniprot_Analysed` WHERE `Uniprot_ID` LIKE '" + Gene_Name_global + "' ")
        FASTA = pd.DataFrame(mycursor.fetchall())
        if FASTA.empty:
            mycursor.execute("SELECT FASTA FROM `Trembl_Entries` WHERE `Uniprot_ID` LIKE '" + Gene_Name_global + "' ")
            FASTA = pd.DataFrame(mycursor.fetchall())
        FASTA = FASTA.iloc[0][0]
    except:
        FASTA = None

    mydb.disconnect()
    mycursor.close()
    return FASTA


def Within_limits(point_to_analyse, z, q):
    if point_to_analyse > int(z) and point_to_analyse < int(q):
        return 1
    else:
        return 0


def analysed_domain_coverage(Domain_name, Domain_start, Domain_finish, Protein_entries_experiment, sequence):
    susceptibility_of_domain = 0
    total_amino_acids = []
    peptides_found = []

    for number3 in Protein_entries_experiment.iterrows():
        peptide_sequence = number3[1]['Peptide']
        peptide_abundance = number3[1]['spectra']

        length_of_peptide = len(peptide_sequence)
        Sequence_Positions = [m.start() for m in re.finditer(peptide_sequence, sequence)]
        for start_position in Sequence_Positions:
            # loops through each of the start positions, detewcted in the sequence
            value1 = start_position
            value2 = value1 + length_of_peptide
            susceptibility_of_domain_start = Within_limits(value1, Domain_start, Domain_finish)
            susceptibility_of_domain1_end = Within_limits(value2, Domain_start, Domain_finish)
            if susceptibility_of_domain1_end == 0 and susceptibility_of_domain_start == 0:
                pass
            else:
                if susceptibility_of_domain_start == 1 and susceptibility_of_domain1_end == 0:
                    susceptibility_of_domain = susceptibility_of_domain + peptide_abundance
                    peptide_coverage = Domain_finish - value1
                    amino_acids_covered = list(range(value1, int(Domain_finish)))
                    total_amino_acids = total_amino_acids + amino_acids_covered
                    peptides_found.append(peptide_sequence)
                elif susceptibility_of_domain_start == 0 and susceptibility_of_domain1_end == 1:

                    susceptibility_of_domain = susceptibility_of_domain + peptide_abundance
                    peptide_coverage = (value2 - Domain_start)
                    amino_acids_covered = list(range(int(Domain_start), value2))
                    peptides_found.append(peptide_sequence)
                    if type(amino_acids_covered) != int:
                        total_amino_acids = total_amino_acids + amino_acids_covered
                elif susceptibility_of_domain_start == 1 and susceptibility_of_domain1_end == 1:
                    susceptibility_of_domain = susceptibility_of_domain + peptide_abundance
                    peptide_coverage = length_of_peptide
                    amino_acids_covered = list(range(value1, value2))
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


def Get_Domains_SQL(AC):
    DOMAINS = pd.DataFrame()
    mydb = mysql.connector.connect(
        host=HOST,
        user=USER,
        passwd=PASSWORD,
        database="Uniprot",
        auth_plugin='mysql_native_password'
    )

    mycursor = mydb.cursor()
    try:

        mycursor.execute(
            f"SELECT `Name`,`finish`,`start`,`Type` FROM `Domains_Uniprot_all` WHERE `Uniprot_ID` LIKE '{AC}'")
        df_with_doamin_info = pd.DataFrame(mycursor.fetchall())
        if not df_with_doamin_info.empty:
            field_names = [i[0] for i in mycursor.description]
            df_with_doamin_info.columns = field_names

        if not df_with_doamin_info.empty:
            field_names = [i[0] for i in mycursor.description]
            df_with_doamin_info.columns = field_names
        else:
            df_with_doamin_info = pd.DataFrame()
        DOMAINS = df_with_doamin_info
    except:
        DOMAINS = pd.DataFrame()

    mydb.disconnect()
    mycursor.close()
    return DOMAINS


def Get_all_Experiment_Domains(Gene_Name_global, sample, Analysis_Type, Results):

    Domain_lengths = []
    mydb = mysql.connector.connect(
        host="localhost",
        user="root",
        passwd="Weafrae1",
        database="Experiments_MS"
    )
    mycursor = mydb.cursor()
    sql_select_Query = "SELECT `Domain_Name`, `NumberOfSpectra`, `Exclusive_spectrum_count`,`Domain_Start`,`Domain_Finish`,`Domain Type` FROM `Counting` WHERE `GeneAC` = '" + Gene_Name_global + "'  AND `experiment_name` = '" + sample + "' AND `Domain Type` = '" + Analysis_Type + "';";
    mycursor.execute(sql_select_Query)
    Data = pd.DataFrame(mycursor.fetchall())
    mydb.disconnect()
    mycursor.close()


def Master_Run_Counting_Algorythm_Clean(Protein=None, Domain_Types=None, Protein_peptides=None,Reference_Proteome=None,Reference_Domains=None):

    Fasta = Reference_Proteome[Reference_Proteome.Uniprot_ID==Protein]["FASTA"][0]
    domains=Reference_Domains[Reference_Domains.Uniprot_ID == Protein]

    
    if domains.empty:
        pass
    else:
        domains.Name=domains.Name.str.replace("'", "")
    domain_ranges = [float(s.replace(" AA STEP", "")) for s in Domain_Types if "AA STEP" in s]
    for range_elem in domain_ranges:
        dom= Fasta_Analysis_Arbitarely_Domains(sequence=Fasta,step_size=range_elem)
        dom.Name = str(int(range_elem)) + "Step_" + dom.Name
        domains = pd.concat([domains, dom])

    # Here filter out the data
    if Domain_Types != None:
        domains = domains[domains.Type.isin(Domain_Types)]

    domains.loc[domains.Name.duplicated(), "Name"] = domains[domains.Name.duplicated()].Name + "_" + domains[
        domains.Name.duplicated()].start.astype(int).astype(str) + "_" + domains[
                                                         domains.Name.duplicated()].finish.astype(int).astype(str)

    Protein_Entries = Protein_peptides[
        Protein_peptides['Protein'].str.contains(Protein, na=False)]
    experiment_names = Protein_Entries['Sample'].unique()
    Experiment_dict = {}
    Experiment_Coverages = pd.DataFrame()
    for experiment in experiment_names:
        Protein_entries_experiment = Protein_Entries[Protein_Entries['Sample'] == experiment]
        Exclusive_spectrum_count = Protein_entries_experiment['spectra'].sum()
        for index, row in domains.iterrows():
            Domain_name = row.Name
            Domain_start = row.start
            Domain_finish = row.finish
            Domain_type = row.Type
            try:
                susceptibility_of_domain, Domain_name, percentage_covered, peptides_found = analysed_domain_coverage(
                    Domain_name, Domain_start, Domain_finish, Protein_entries_experiment, Fasta)

                Experiment_Coverages = Experiment_Coverages.append(
                    {'Domain_Name': Domain_name,  #
                     'Domain_Start': Domain_start,  #
                     'Domain_Finish': Domain_finish,  #
                     'Domain Type': Domain_type,  #
                     'NumberOfSpectra': susceptibility_of_domain,  #
                     'Percent_Covered': percentage_covered,  #
                     'Exclusive_spectrum_count': Exclusive_spectrum_count,  #
                     'experiment_name': experiment,  #
                     # 'Experiment_Setup': domains , #????
                     'GeneAC': Protein,  #
                     'peptides_found': peptides_found  #
                     }, ignore_index=True)
            except:
                continue
        Experiment_dict[experiment] = Experiment_Coverages
    return Experiment_Coverages, Fasta


def analyse_sample(key, value2, Results):
    Spectral_total_counts = pd.DataFrame(index=value2)
    experiment_setup_pandas = Results[
        ['Domain_Name', 'Domain_Start', 'Domain_Finish', 'Domain Type']].drop_duplicates()
    experiment_setup_pandas = experiment_setup_pandas.set_index(experiment_setup_pandas.Domain_Name)
    experiment_setup_pandas[
        'Domain_Length'] = experiment_setup_pandas.Domain_Finish - experiment_setup_pandas.Domain_Start
    Spectral_total_counts[key] = [0] * value2.__len__()
    # experiment_setup_pandas_Peptide_Counts = experiment_setup_pandas.copy ( )

    col = 'Peptides ' + key
    per_domain_peptides = pd.DataFrame([''] * experiment_setup_pandas.__len__(), columns=[col]).set_index(
        experiment_setup_pandas.index)
    i = -1

    for sample in value2:
        i = i + 1
        ##Get only this sample entries
        experiment_setup_pandas[sample] = 0.0
        # experiment_setup_pandas_Peptide_Counts [ sample ] = 0.0
        records = Results[Results.experiment_name == sample]
        records = records.set_index(records.Domain_Name)
        if not records.empty:
            Domain_Spectral_Count = records.NumberOfSpectra
            Total_spectral_count = records.Exclusive_spectrum_count.reset_index().Exclusive_spectrum_count[0]
            Spectral_total_counts[key][sample] = Total_spectral_count

            for i, value in enumerate(records.peptides_found):
                Domain = records.Domain_Name.iloc[i]
                if value == '':
                    number_of_peptides = 0
                else:
                    # number_of_peptides = value.split ( "," ).__len__ ( )
                    value = value.replace(",", " ")
                per_domain_peptides['Peptides ' + key][Domain] = per_domain_peptides['Peptides ' + key][
                                                                     Domain] + ' ' + value
                # experiment_setup_pandas_Peptide_Counts [ sample , Domain ] = number_of_peptides
            experiment_setup_pandas[sample] = Domain_Spectral_Count
    return Spectral_total_counts, experiment_setup_pandas, per_domain_peptides


def normalise(Spectra_Matrix, Spectral_total_counts):
    Spectra_Matrix_output = Spectra_Matrix.copy()

    Spectra_Values = Spectra_Matrix.iloc[:, 5:]
    for column in Spectra_Values:        # print(column)
        if column in Spectral_total_counts.index:
            # print("yes")
            Factor = float(Spectral_total_counts.loc[column,"Factor"])
            value_for_Analysis = round(Spectra_Matrix[column] * Factor)
            Spectra_Matrix_output[column] = value_for_Analysis
        else:
            # print("no")
            print("There is no sample present in the: "+column)


    return Spectra_Matrix_output


def determine_differences(Norm_Stats_Spectra, Domain_lengths):
    from itertools import combinations
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


def drop_repeats(domains, Results):
    import collections
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
            if len(x) < len(y):
                to_drop.append(combo[0])
            else:
                to_drop.append(combo[1])
        else:
            continue

    ctr = collections.Counter(to_drop)
    for elem in ctr.items():
        if elem[1] > 1:
            domains = domains[domains.Domain_Name != elem[0]]
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


def Master_Run_Structural_Analysis(experiment_feed=None, Results=None, Protein=None,paired=True):
    Data2 = pd.DataFrame()
    Unique_Domains = Results['Domain Type'].unique()
    Datafile = Results
    DataVal = []
    Spectral_total_counts = Results[["experiment_name", "Exclusive_spectrum_count"]].drop_duplicates()
    value=[]

    # adding 0 to the normalisation
    for key,value2 in experiment_feed.items():
        value.extend(value2)
    for value1 in value:
        # print(value1)
        if Spectral_total_counts[(value1 == Spectral_total_counts.experiment_name)].experiment_name.count()>0:
            if (float(Spectral_total_counts[(value1 == Spectral_total_counts.experiment_name)].Exclusive_spectrum_count)==0.0):
                Spectral_total_counts[(value1 == Spectral_total_counts.experiment_name)].Exclusive_spectrum_count=0.0000001
            else:
                continue
        
        else:
            Spectral_total_counts=Spectral_total_counts.append({"experiment_name": value1, "Exclusive_spectrum_count": 0.0000001}, ignore_index=True)
    # here have to check whether all the samples are in the file, if not we add 0 as a norm factor.

    Spectral_total_counts = Spectral_total_counts.set_index("experiment_name")

    Median_Norm_Factor = statistics.median(
        [Spectral_total_counts["Exclusive_spectrum_count"].max(),
         Spectral_total_counts["Exclusive_spectrum_count"].min()])
         
    Spectral_total_counts["Median_Norm"] = Median_Norm_Factor
    Spectral_total_counts["Factor"] = Spectral_total_counts["Median_Norm"] / Spectral_total_counts[
        "Exclusive_spectrum_count"]
    # loop through each of the domain types
    for Data_Type in Unique_Domains:
        # f = open("log_file.log", "a")
        # f.write(f"Analysing: {Data_Type}\n")
        # f.close()
        # try:
        Results = Datafile[Datafile['Domain Type'] == Data_Type]
        domains = Results[["Domain_Name", "Domain_Finish", "Domain_Start"]]
        domains = domains.drop_duplicates()

        # Here first check the number of domains and then
        if domains.Domain_Name.__len__() > 1:
            # Here Check for the overlap of the DataSets.
            Results = drop_repeats(domains, Results)
            if not Results.empty:
                experiments_all_Spectra = {}
                experiments_all_counts = {}
                experiments_all_Peptides = {}
                experiments_all_Spectra_For_Stats = {}

                # get the spectral total counts for different samples.
                Dataframes =pd.DataFrame()
                for key, value2 in experiment_feed.items():
                    # had to fix here, we have now used the toatal sample normalisation
                    _, Spectra_Matrix, Peptides = analyse_sample(key, value2,
                                                                 Results)

                    Norm_Stats_Spectra = normalise(Spectra_Matrix,Spectral_total_counts)
                    experiments_all_Spectra[key] = Norm_Stats_Spectra
                    experiments_all_Spectra_For_Stats[key] = Norm_Stats_Spectra.drop(
                        columns=["Domain_Start", "Domain_Finish", "Domain Type", "Domain_Length"])
                    Original = Spectra_Matrix.iloc[:,5:].add_prefix("Original_")
                    Dataframes=pd.concat([Dataframes,Original], axis=1, sort=False)
                    experiments_all_Peptides[key] = Peptides

                Domain_lengths = Norm_Stats_Spectra.Domain_Length
                Differences, Averages = determine_differences(experiments_all_Spectra, Domain_lengths)
                Peptides_PD = keys_to_Pandas(experiments_all_Peptides)
                All_Spectra = keys_to_Pandas(experiments_all_Spectra).drop_duplicates()
                # All_Spectra.to_csv(
                #     f"/run/user/1000/gvfs/smb-share:server=10.2.82.9,share=bmhrss$/snapped/replicated/Sherratt_Lab/Matiss Ozols/Structural_Debug/{Protein}.csv")
                Averages = pd.DataFrame(Averages)
                p_values_adjusted = Two_Way_mixed_Anova(experiments_all_Spectra_For_Stats,paired=paired)
                p_values = p_values_adjusted.set_index(p_values_adjusted.Domain_Name)
                p_values.fillna(value=1,inplace=True)
                Data_Values = p_values_adjusted.iloc[:, 1:]
                DataVal.extend(Data_Values.values)
                # if((Data_Values <0.05).any()[0]):
                Data = pd.concat([All_Spectra, Differences, p_values, Averages, Peptides_PD,Dataframes], axis=1)
                Data2 = Data2.append(Data)

                # else:
                #     print("No p values with significance!")
            else:
                # print(Data_Type)
                print("result empty")
        else:
            print("Only 1 domain avaliable. Cant do stats on 1 Domain")

        # except:
        #     print(Data_Type)
        #     print("failed with domains!")
    if DataVal.__len__()>0:
        if (pd.DataFrame(DataVal) < 0.05).any()[0]:
            Data2['GeneAC'] = Protein
            # Data2 = pd.concat([Data2, Dataframes], axis=1, sort=False)
            Data3 = Data2.loc[:, ~Data2.columns.duplicated()]
            return Data3, Spectral_total_counts
        else:
            print("No p values with significance here!")
            return [], []
    else:
        print("No p values with significance here!")
        return [],[]
    # here have to check whether there are any p values that are significant.


def Master_Run_Score_Calculations(Structural_Analysis_Results, Protein):
    # import math

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
