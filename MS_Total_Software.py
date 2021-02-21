##This is the new Python 3 version of the software to analyse the MS Files
import pandas as pd

import sys
from sqlalchemy import create_engine

from Structural.Structural_Algorythm.Functions_Clean import retrieve_reviewed, Master_Run_Counting_Algorythm_Clean, \
    Master_Run_Structural_Analysis, Master_Run_Score_Calculations
import json
import mysql.connector
from MSP.settings import PASSWORD, USER, HOST, DB

n = 15
k = 0
record_every_nr_steps =1
connection_db = mysql.connector.connect(host=HOST,
                                        database=DB,
                                        user=USER,
                                        password=PASSWORD,
                                        auth_plugin='mysql_native_password',
                                        )
cursor_query = connection_db.cursor()


def record_data(Structural_Json,Coverage_Json,Owner_ID,id):
    Significant_Protein_Count = Structural_Json.keys().__len__()
    experiment_coverages = json.dumps(Coverage_Json)
    structural_analysis_results = json.dumps(Structural_Json)
    sql_query = f"UPDATE `Structural_userdata` SET Structural_Results=JSON_MERGE_PATCH(`Structural_Results`, '{structural_analysis_results}')," \
                f" Experimental_Coverages=JSON_MERGE_PATCH(`Experimental_Coverages`,'{experiment_coverages}')," \
                f" Progress='Analysing', Significant_Protein_Count=`Significant_Protein_Count`+{Significant_Protein_Count} " \
                f"WHERE (`id` like {id} and `owner_id` LIKE {Owner_ID})"
    cursor_query.execute(sql_query)
    connection_db.commit()
    # clear the memory
    Structural_Json = {}
    Coverage_Json = {}
    print("successfuly recorded significant domain")
    return Structural_Json, Coverage_Json

def run_full_analysis(Proteins_to_Analyse, Domain_types, Protein_peptides, experiment_feed, Owner_ID, id, paired=False):
    # If we do decide to remove the protein entry then we have to look up each peptide in the library and find all the peptides for the protein thatr are provided.

    k_val=0
    # paired=True
    Structural_Json = {}
    Coverage_Json = {}
    Experimental_coverages_all = {}

    for i, entry in enumerate(Proteins_to_Analyse):
        # if i>50:
        #     break
        try:
            print(f"Here: {i}: {entry}")
            # if "F213A_HUMAN" in entry:
            Protein_list = entry.split(',')
            Protein = retrieve_reviewed(Protein_list)

            # Domain_types=["50.0 AA STEP"] #experiment coverages is a pure peptide counts without normalisation
            Experiment_Coverages, Fasta = Master_Run_Counting_Algorythm_Clean(Gene_Name_global=Protein,
                                                                              Domain_Types=Domain_types,
                                                                              Protein_peptides=Protein_peptides,Protein_list=Protein_list)

            Structural_Analysis_Results, Norm_Factors = Master_Run_Structural_Analysis(experiment_feed=experiment_feed,
                                                                                       Results=Experiment_Coverages,
                                                                                       Protein=Protein, paired=paired)

            if Structural_Analysis_Results.__len__()>0:
                Structural_Analysis_Results.drop("GeneAC", axis=1, inplace=True)
                Experiment_Coverages.drop("GeneAC", axis=1, inplace=True)
                k_val = k_val + 1
                Experiment_Coverages = Experiment_Coverages[
                    Experiment_Coverages["Domain Type"].isin(Structural_Analysis_Results["Domain Type"].unique())]
                Experiment_Coverages.set_index("Domain_Name", drop=False, inplace=True)
                Coverage_Json[Protein] = Experiment_Coverages.to_dict()
                Structural_Json[Protein] = {"Data": Structural_Analysis_Results.to_dict('index'),
                                            "Norm_Factors": Norm_Factors.to_dict()}

                if k_val % record_every_nr_steps == 0:
                    Structural_Json, Coverage_Json=record_data(Structural_Json, Coverage_Json, Owner_ID,id)
            else:
                print("Next: no p values to record")

        except:

            # Structural_Json = {}
            # Coverage_Json = {}
            print(sys.exc_info()[0])
            continue
        # else:
        #     continue
    if (Structural_Json.keys().__len__()>0):
        Structural_Json, Coverage_Json=record_data(Structural_Json, Coverage_Json, Owner_ID,id)

    connection_db.disconnect()
    cursor_query.close()
    return "success"

# if __name__ == '__main__':
#
#     file4='Skin_Epidermis_' #no need
#     Data_Table = 'Experiment_setup2' #no need
#     file_name='Epidermis_Structural3.csv' #no need #'Peptide Report for Elastase_HDF Fbn MF_SSR_UVB5.csv'#'Peptide Report for Dermis_Full_Skin_MS.csv'
#     database='Experiments_MS'  #no need
#     table_prefix = file_name.split('.')[0] #no need
#     table_prefix=table_prefix.replace(" ","");table_prefix=table_prefix.replace("_","")+" " #no need
#     Protein_peptides = pd.read_csv('/home/matiss/Documents/Projects/Structural_Analysis/Working_Data_File.csv', index_col=False, sep=',', ) #no need
#     Experiments = Protein_peptides['MS/MS sample name'].unique()
#
#     #Epidermis
#     experiment_feed = {'Buttock':['20180601_SherrattM_MO_01.raw (Full_Skin_1)','20180601_SherrattM_MO_03.raw (Full_Skin_3)','20180601_SherrattM_MO_05.raw (Full_Skin_5)','20180601_SherrattM_MO_07.raw (Full_Skin_7)','20180601_SherrattM_MO_09.raw (Full_Skin_9)','20180601_SherrattM_MO_11.raw (Full_Skin_11)','20180601_SherrattM_MO_13.raw (Full_Skin_13)'],
#                        'Forearm':['20180601_SherrattM_MO_02.raw (Full_Skin_2)','20180601_SherrattM_MO_04.raw (Full_Skin_4)','20180601_SherrattM_MO_06.raw (Full_Skin_6)','20180601_SherrattM_MO_08.raw (Full_Skin_8)','20180601_SherrattM_MO_10.raw (Full_Skin_10)','20180601_SherrattM_MO_12.raw (Full_Skin_12)','20180601_SherrattM_MO_14.raw (Full_Skin_14)']} #epidermis
#
#
#     Proteins_to_Analyse = pd.Series(Protein_peptides['Protein accession numbers'].unique()).dropna()
#     Reviewed_proteins = (pd.read_csv('Uploaded_MS_files/Reviewed', index_col=False))
#
#     Domain_types = ['REGIONS','TOPO_DOM','DOMAINS','REPEAT','50.0 AA STEP']
#     Domain_types = ['REGIONS']
#     engine = create_engine("mysql+pymysql://{user}:{pw}@localhost/{db}"
#                            .format(user="root",
#                                    pw="Weafrae1",
#                                    db="Experiments_MS"))
#
#     Proteins_to_Analyse=["APOE_HUMAN"]
#     run_full_analysis(Proteins_to_Analyse)
