import pandas as pd
import sys

from Functions_Clean import retrieve_reviewed, Master_Run_Counting_Algorythm_Clean, \
    Master_Run_Structural_Analysis, Master_Run_Score_Calculations
import json
Structural_Json = {}
Coverage_Json = {}
Experimental_coverages_all = {}
Reference_Proteome=None
Reference_Domains=None
# Protein_peptides2=pd.DataFrame(columns=['Sample', 'Peptide', 'Protein', 'spectra'])
# Protein_peptides2=pd.DataFrame(columns=['Sample', 'Peptide', 'Protein', 'spectra'])
Protein_peptides2=[]
# Reference_Proteome=[]
# Reference_Domains=[]

def collect_result(result):
    Protein1=result[2]
    Coverage_Json1=result[0]
    Structural_Json1=result[1]
    Coverage_Json[Protein1] = Coverage_Json1
    Structural_Json[Protein1] = Structural_Json1

def record_data(Structural_Json,Coverage_Json,Owner_ID,id):
    import mysql.connector
    from secret import HOST, PORT, PASSWORD, DB, USER
    connection = mysql.connector.connect(host=HOST,
                                        database=DB,
                                        user=USER,port=PORT,
                                        password=PASSWORD,
                                        auth_plugin='mysql_native_password')
    cursor_query = connection.cursor()
    Significant_Protein_Count = Structural_Json.keys().__len__()
    experiment_coverages = json.dumps(Coverage_Json)
    structural_analysis_results = json.dumps(Structural_Json)
    # sql_query = f"UPDATE `Structural_userdata` SET Structural_Results=JSON_MERGE_PATCH(`Structural_Results`, '{structural_analysis_results}')," \
    #             f" Experimental_Coverages=JSON_MERGE_PATCH(`Experimental_Coverages`,'{experiment_coverages}')," \
    #             f" Progress='Analysing', Significant_Protein_Count=`Significant_Protein_Count`+{Significant_Protein_Count} " \
    #             f"WHERE (`id` like {id} and `owner_id` LIKE {Owner_ID})"

    sql_query = f"UPDATE `Structural_userdata` SET Structural_Results='{structural_analysis_results}'," \
            f" Progress='Finished', Significant_Protein_Count=`Significant_Protein_Count`+{Significant_Protein_Count} " \
            f"WHERE (`id` like {id} and `owner_id` LIKE {Owner_ID})"

    cursor_query.execute(sql_query)
    connection.commit()
    connection.disconnect()
    cursor_query.close()
    # clear the memory
    # Structural_Json = {}
    # Coverage_Json = {}
    print("successfuly recorded significant domain")
    return Structural_Json, Coverage_Json

def run_protein(Protein,Reference_Proteome,Reference_Domains,Domain_types,Protein_peptides,experiment_feed,paired):
    print(f"Here: : {Protein}")

    Experiment_Coverages, Fasta = Master_Run_Counting_Algorythm_Clean(Protein=Protein,
                                                                    Domain_Types=Domain_types,
                                                                    Protein_peptides=Protein_peptides,Reference_Proteome=Reference_Proteome,
                                                                    Reference_Domains=Reference_Domains)

    Structural_Analysis_Results, Norm_Factors = Master_Run_Structural_Analysis(experiment_feed=experiment_feed,
                                                                            Results=Experiment_Coverages,
                                                                            Protein=Protein, paired=paired)


    if Structural_Analysis_Results.__len__()>0:
        Structural_Analysis_Results.drop("GeneAC", axis=1, inplace=True)
        Experiment_Coverages.drop("GeneAC", axis=1, inplace=True)
        # k_val = k_val + 1
        Experiment_Coverages = Experiment_Coverages[
            Experiment_Coverages["Domain Type"].isin(Structural_Analysis_Results["Domain Type"].unique())]
        Experiment_Coverages.set_index("Domain_Name", drop=False, inplace=True)


        Coverage = Experiment_Coverages.to_dict()
        Structural_Json= {"Data": Structural_Analysis_Results.to_dict('index'),
                                    "Norm_Factors": Norm_Factors.to_dict()}
        return (Coverage,Structural_Json,Protein)


        # if k_val % record_every_nr_steps == 0:
        #     Structural_Json, Coverage_Json=record_data(Structural_Json, Coverage_Json, Owner_ID,id)
    
    else:
        print("Next: no p values to record")
        next()
        
def run_full_analysis( Domain_types, Protein_peptides, experiment_feed, Owner_ID, id, paired=False, Spiecies="HUMAN"):
    # If we do decide to remove the protein entry then we have to look up each peptide in the library and find all the peptides for the protein thatr are provided.
    #
    
    Reference_Proteome = pd.read_csv(f"./outputs/Uniprot_{Spiecies}.tsv",sep="\t",index_col=0)
    Reference_Domains = pd.read_csv(f"./outputs/Domains_Uniprot_{Spiecies}.tsv",sep="\t",index_col=0)

    # Protein_peptides=match_peptide_to_protein(Protein_peptides,Reference_Proteome)
    # Protein_peptides.to_csv("Protein_peptides_Mouse_Aorta.tsv",sep="\t")
    Protein_peptides=pd.read_csv("Protein_peptides_Mouse_Aorta.tsv",sep="\t",index_col=0)
    k_val=0
    # paired=True
    import multiprocessing as mp
    # pool.close()
    pool = mp.Pool(mp.cpu_count()-1)
     # pool = mp.Pool(1)
    
    Protein_peptides=Protein_peptides.dropna(subset=['spectra']) # Drop the NaN vales on spectra. ie - peptides are not detected in that sample
    for i, Protein in enumerate(Protein_peptides.Protein.str.split(";").explode().unique()):
        try:
            pool.apply_async(run_protein, args=([Protein,Reference_Proteome,Reference_Domains,Domain_types,Protein_peptides,experiment_feed,paired]),callback=collect_result) #paralel runs - uses all the cores available
            # run_protein(Protein,Reference_Proteome,Reference_Domains,Domain_types,Protein_peptides,experiment_feed,paired) #individual cores - uses 1 core hence 1 protein at a time is analysed.
        except:
            print(sys.exc_info()[0])
            continue
        if i>10:
            break
    pool.close()
    pool.join() 
    # print(Coverage_Json)


    with open("./bin/Coverage_Json.json", 'w') as json_file:
        json.dump(Coverage_Json, json_file)
    with open("./bin/Structural_Json.json", 'w') as json_file:
        json.dump(Structural_Json, json_file)

    # with open("./bin/Structural_Json.json", 'r') as myfile:
    #     Structural_Json2=myfile.read()
    # Structural_Json2=json.loads(Structural_Json2)
    # here we wait for the process to finish and then HPC should send the data back.

    record_data(Structural_Json, Coverage_Json, Owner_ID,id)
    # # this is to record last entries that are not processed
    # if (Structural_Json.keys().__len__()>0):
    #     Structural_Json, Coverage_Json=record_data(Structural_Json, Coverage_Json, Owner_ID,id)

    # connection_db.disconnect()
    # cursor_query.close()
    return "success"
    
def retrieve_all_proteins(peptide,Reference_Proteome,Protein_peptides,i):
    # print(f"{peptide} {i}")
    Proteins_containing_peptide = Reference_Proteome[Reference_Proteome.FASTA.str.contains(peptide)]["Uniprot_ID"]
    All_Proteins = ";".join(Proteins_containing_peptide)
    
    # All_Prots = Protein_peptides[Protein_peptides.Peptide == peptide]
    # All_Prots["Protein"] = All_Proteins

    return All_Proteins

def append_results(result):
    for i,row in result.iterrows():
        Protein_peptides2.append(row)

def match_peptide_to_protein(Protein_peptides,Reference_Proteome):
    # limit to the protein of interest.
    # Reference_Proteome = Reference_Proteome[Reference_Proteome.Uniprot_Type=="Uniprot"]

    Protein_peptides["Protein"]=""
    # sometimes data outputs contain a peptide sequence in brackets - the folowing will remove this
    Protein_peptides.Peptide=Protein_peptides.Peptide.str.replace(".*\]\.",'')
    Protein_peptides.Peptide=Protein_peptides.Peptide.str.replace("\.\[.*","")

    # import multiprocessing as mp
    # pool = mp.Pool(mp.cpu_count()-1)

    count=0
    # from time import process_time
    # t1_start = process_time()

    for peptide in Protein_peptides.Peptide.unique():
        count+=1
        # len(Protein_peptides.Peptide.unique())
        All_Proteins = retrieve_all_proteins(peptide,Reference_Proteome,Protein_peptides,count)
        Protein_peptides.loc[Protein_peptides.Peptide == peptide,"Protein"]=All_Proteins

        # pool.apply_async(retrieve_all_proteins, args=([peptide,Reference_Proteome,Protein_peptides,count]),callback=append_results) #paralel runs - uses all the cores available
        # if count>100:
        #     break
    
    
    # pool.close()
    # pool.join() 

    # Protein_peptides=pd.DataFrame(Protein_peptides2)

    # t1_stop = process_time() 
    # print("Elapsed time:", t1_stop, t1_start)  
   
    # print("Elapsed time during the whole program in seconds:", 
    #                                         t1_stop-t1_start)  

    return Protein_peptides

def retrieve_mysql_data():

    # import requests
    # response = requests.get("http://www.manchesterproteome.manchester.ac.uk/run_api/MSP_api/?page=1&search=")
    # print(response.json())
    # d=response.json()
    # with open("test_output.json", 'w') as json_file:
    #     json.dump(d, json_file)


    import mysql.connector
    from secret import HOST, PORT, PASSWORD, DB, USER
    connection = mysql.connector.connect(host=HOST,
                                        database=DB,
                                        user=USER,port=PORT,
                                        password=PASSWORD,
                                        auth_plugin='mysql_native_password')
    cursor_query = connection.cursor()

    # here could select the jobs that are qued - do this every 6h and if a new job is qued then process on the HPC cloud
    sql="SELECT id FROM `Structural_userdata` WHERE Progress LIKE 'Que'"
    cursor = connection.cursor()
    cursor.execute(sql)
    Data_ids = pd.DataFrame(cursor.fetchall())
    if not Data_ids.empty:
        field_names = [i[0] for i in cursor.description]
        Data_ids.columns = field_names

    sql=f"SELECT * FROM `Structural_userdata` WHERE `id` LIKE '{Data_ids.id[1]}'"
    cursor = connection.cursor()
    cursor.execute(sql)
    Data = pd.DataFrame(cursor.fetchall())
    if not Data.empty:
        field_names = [i[0] for i in cursor.description]
        Data.columns = field_names
    Domain_types=json.loads(Data.Domain_Types[0])
    experiment_feed = json.loads(Data.experimental_design[0])
    Owner_ID = Data.owner_id[0]
    id = Data.id[0]
    Data_Val = Data.data[0]
    Data_Val = json.loads(Data_Val)
    Protein_peptides = pd.DataFrame(Data_Val)
    Paired= Data.Paired[0]
    Spiecies="MOUSE"
    run_full_analysis(Domain_types, Protein_peptides, experiment_feed,Owner_ID,id,paired=Paired, Spiecies=Spiecies)

def retrieve_save_and_process():
    import mysql.connector
    from secret import HOST, PORT, PASSWORD, DB, USER
    connection = mysql.connector.connect(host=HOST,
                                        database=DB,
                                        user=USER,port=PORT,
                                        password=PASSWORD,
                                        auth_plugin='mysql_native_password')
    cursor_query = connection.cursor()

    # here could select the jobs that are qued - do this every 6h and if a new job is qued then process on the HPC cloud
    sql="SELECT id FROM `Structural_userdata` WHERE Progress LIKE 'Que'"
    cursor = connection.cursor()
    cursor.execute(sql)
    Data_ids = pd.DataFrame(cursor.fetchall())
    if not Data_ids.empty:
        field_names = [i[0] for i in cursor.description]
        Data_ids.columns = field_names
    print(Data_ids)
    Data_ids.to_csv("tmp_ids.csv")


def retrieve_mysql_data_test():
    import requests
    response = requests.get('http://www.manchesterproteome.manchester.ac.uk/run_api/MSP_api/?page=1&search=')
    d = response.json()
    with open("test_output.json", 'w') as json_file:
        json.dump(d, json_file)

    print("Done")
        



if __name__ == '__main__':
    # retrieve_mysql_data_test()
    retrieve_save_and_process()
    # import json
    # paired=1
    # id='1'
    # Owner_ID='1'
    # with open("./sample_input/experiment_feed.json", 'r') as myfile:
    #     experiment_feed=myfile.read()
    # experiment_feed=json.loads(experiment_feed)
    # Protein_peptides=pd.read_csv("./sample_input/Protein_peptides.tsv",sep="\t",index_col=0)
    # Domain_types=pd.read_csv("./sample_input/Domain_types.tsv",sep="\t",index_col=0,names=["Dom"],header=1)
    # Domain_types=Domain_types.Dom.values.tolist()
    # Protein_peptides=pd.read_csv("tmp_working_file.csv",sep="\t",index_col=0)
    # run_full_analysis(Domain_types, Protein_peptides, experiment_feed,Owner_ID,id,paired=paired)
