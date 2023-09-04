import pandas as pd
import sys
import mysql.connector
from Functions_Clean import MPLF_Domain_Quantifications, MPLF_Statistical_Analyisis
import json
Structural_Json = {}
Coverage_Json = {}
Experimental_coverages_all = {}
Reference_Proteome=None
Reference_Domains=None
# Protein_peptides2=pd.DataFrame(columns=['Sample', 'Peptide', 'Protein', 'spectra'])
# Protein_peptides2=pd.DataFrame(columns=['Sample', 'Peptide', 'Protein', 'spectra'])
Protein_peptides2={}
# Reference_Proteome=[]
# Reference_Domains=[]

def collect_result(result):
    print("collected")
    Protein1=result[2]
    Coverage_Json1=result[0]
    Structural_Json1=result[1]
    Coverage_Json[Protein1] = Coverage_Json1
    Structural_Json[Protein1] = Structural_Json1

def record_data(Structural_Json,Owner_ID,id,Domain_types):
    import mysql.connector
    from secret import HOST, PORT, PASSWORD, DB, USER
    connection = mysql.connector.connect(host=HOST,
                                        database=DB,
                                        user=USER,port=PORT,
                                        password=PASSWORD,
                                        auth_plugin='mysql_native_password')
    cursor_query = connection.cursor()
    Significant_Protein_Count = Structural_Json.keys().__len__()
    # experiment_coverages = json.dumps(Coverage_Json)
    structural_analysis_results = json.dumps(Structural_Json)
    Domain_types1=json.dumps(Domain_types)
    # sql_query = f"UPDATE `Structural_userdata` SET Structural_Results=JSON_MERGE_PATCH(`Structural_Results`, '{structural_analysis_results}')," \
    #             f" Experimental_Coverages=JSON_MERGE_PATCH(`Experimental_Coverages`,'{experiment_coverages}')," \
    #             f" Progress='Analysing', Significant_Protein_Count=`Significant_Protein_Count`+{Significant_Protein_Count} " \
    #             f"WHERE (`id` like {id} and `owner_id` LIKE {Owner_ID})"

    sql_query = f"UPDATE `Structural_userdata` SET Structural_Results='{structural_analysis_results}'," \
            f" Progress='Finished',Domain_types='{Domain_types1}', Significant_Protein_Count='{Significant_Protein_Count}' " \
            f"WHERE (`id` like {id} and `owner_id` LIKE {Owner_ID})"
    # sql_query = f"UPDATE `Structural_userdata` SET Structural_Results='{structural_analysis_results}'," \
    #         f" Domain_types='{Domain_types1}', Significant_Protein_Count='{Significant_Protein_Count}' " \
    #         f"WHERE (`id` like {id} and `owner_id` LIKE {Owner_ID})"

    cursor_query.execute(sql_query)
    connection.commit()
    connection.disconnect()
    cursor_query.close()
    # clear the memory
    # Structural_Json = {}
    # Coverage_Json = {}
    print("successfuly recorded significant domain")
    return Structural_Json, Coverage_Json

def MPLF(Protein,Reference_Proteome,Reference_Domains,Domain_types,Protein_peptides,experiment_feed,paired):
    ###########
    ## First all the coverages are calculated for each of the experiments within domains
    ## Then the statistical analysis is performed to determine significant hits
    ###########
    print(f"Analysing protein: : {Protein}")
    Experiment_Coverages = MPLF_Domain_Quantifications(Protein=Protein,
                                                                    Domain_Types=Domain_types,
                                                                    Protein_peptides=Protein_peptides,Reference_Proteome=Reference_Proteome,
                                                                    Reference_Domains=Reference_Domains)
        
    Structural_Analysis_Results, Norm_Factors = MPLF_Statistical_Analyisis(experiment_feed=experiment_feed,
                                                                            Results=Experiment_Coverages,
                                                                            Protein=Protein, paired=paired,cuttoff_p_val=0.05)

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

        
def run_full_analysis( Domain_types, Protein_peptides, experiment_feed, Owner_ID=1, id=1, cpus=1,paired=False, Spiecies="HUMAN"):
    # If we do decide to remove the protein entry then we have to look up each peptide in the library and find all the peptides for the protein thatr are provided.
    import os
    Reference_Proteome = pd.read_csv(f"outputs/Uniprot_{Spiecies}.tsv",sep="\t",index_col=0)
    Reference_Domains = pd.read_csv(f"outputs/Domains_Uniprot_{Spiecies}.tsv",sep="\t",index_col=0)
    if not 'Protein' in list(Protein_peptides.columns):
        Protein_peptides=match_peptide_to_protein(Protein_peptides,Reference_Proteome,cpus=cpus)
    elif  Protein_peptides.Protein.unique()[0]=='undefined':
        Protein_peptides=match_peptide_to_protein(Protein_peptides,Reference_Proteome,cpus=cpus)

    # Protein_peptides=replace_with_ids(Protein_peptides,Reference_Proteome)
    Protein_peptides=Protein_peptides.dropna(subset=['spectra']) # Drop the NaN vales on spectra. ie - peptides are not detected in that sample
    Protein_peptides.Protein = Protein_peptides.Protein.str.replace(',',';')
    All_proteins={}
    for i, Protein in enumerate(Protein_peptides.Protein.str.split(";").explode().unique()):
        # Here gather all the unique gene names - all the revirewed entries + each unique non uniprot entry
        Prot1=Reference_Proteome[Reference_Proteome["Uniprot_ID"]==Protein]
        try:
            Gene=Prot1["Uniprot_Gene"][0].split(" {")[0]
        except:
            Gene=Protein
        try:
            Type=Prot1["Uniprot_Type"][0]
        except:
            print(f'failed with {Protein} and hence continuing')
            continue
            # next
        try:
            All_proteins[Gene][Type].append(Protein)
        except:
            try:
                All_proteins[Gene][Type]=[]
                All_proteins[Gene][Type].append(Protein)
            except:
                All_proteins[Gene]={}
                All_proteins[Gene][Type]=[]
                All_proteins[Gene][Type].append(Protein)
    i=0
    
    
    import multiprocessing as mp
    pool = mp.Pool(cpus)
    print(f"CPUS: {cpus}")
    for key in All_proteins.keys():
        try:
            Protein= All_proteins[key]['Uniprot'][0]
        except:
            Protein= All_proteins[key]['Trembl'][0]
        
        if cpus>1:
            pool.apply_async(MPLF, args=([Protein,Reference_Proteome,Reference_Domains,Domain_types,Protein_peptides,experiment_feed,paired]),callback=collect_result) #paralel runs - uses all the cores available
        else:
            result = MPLF(Protein,Reference_Proteome,Reference_Domains,Domain_types,Protein_peptides,experiment_feed,paired)
            collect_result(result)
        i+=1

    pool.close()
    pool.join() 
    # print(Coverage_Json)


    # with open(f"bin/Coverage_Json_{Spiecies}_{Owner_ID}_{id}.json", 'w') as json_file:
    #     json.dump(Coverage_Json, json_file)
    with open(f"bin/Structural_Json_{Spiecies}_{Owner_ID}_{id}.json", 'w') as json_file:
        json.dump(Structural_Json, json_file)

    # here we wait for the process to finish and then HPC should send the data back.

    # Structural_Json2=Structural_Json
    # here have to add a visualisation module as per https://github.com/maxozo/ManchesterProteome/blob/e5fb1a1385b2bf11ddbc514d6ca3f0db6b2f272d/frontend/src/components/Structural/BarChart.js#L888-L890
    record_data(Structural_Json, Owner_ID,id,Domain_types)


    # # this is to record last entries that are not processed
    # if (Structural_Json.keys().__len__()>0):
    #     Structural_Json, Coverage_Json=record_data(Structural_Json, Coverage_Json, Owner_ID,id)

    # connection_db.disconnect()
    # cursor_query.close()
    return "success"
    
def retrieve_all_proteins(peptide,Reference_Proteome):
    Proteins_containing_peptide = Reference_Proteome[Reference_Proteome.FASTA.str.contains(peptide)]["Uniprot_ID"]
    All_Proteins = ";".join(Proteins_containing_peptide)
    return {'peptide':peptide,'All_Proteins':All_Proteins}

def append_results(result):
    Protein_peptides2[result['peptide']]=result['All_Proteins']

def replace_with_ids(Protein_peptides,Reference_Proteome):
    for Prot_string in Protein_peptides.Protein.unique():
        AC_String=[]
        try:
            if Prot_string is None:
                print("skip")
            else:
                for Prot1 in Prot_string.split(","):
                    try:
                        UID=Reference_Proteome[Reference_Proteome.Uniprot_AC.str.contains(Prot1)]["Uniprot_ID"][0]
                        if isinstance(UID, pd.Series):
                            for key in UID:
                                print(key)
                                AC_String.append(key)
                        else:
                            AC_String.append(UID)
                    except:
                        AC_String.append(Prot1)
                New_String = ";".join(AC_String)
                Protein_peptides.loc[Protein_peptides.Protein==Prot_string,"Protein"]=New_String
        except:
            print("not handled exceptation")
    return Protein_peptides

class Determine_Peptide_Protein_Identities:
    # This class takes in the peptide sequences
    # and determines the protein it may belong to
    def __init__(self,peptides,Reference_Proteome,cpus=1):
        self.peptides = peptides
        self.Reference_Proteome = Reference_Proteome
        self.protein_peptides={}
        self.cpus = cpus
        
    def append_results(self,result):
        # We append the results here
        self.protein_peptides.update(result)
        
    def retrieve_all_proteins(self,peptide,Reference_Proteome):
        # Here we find all the protein IDs that contain the peptide sequence
        Proteins_containing_peptide = Reference_Proteome[Reference_Proteome['FASTA'].str.contains(peptide)]["Uniprot_ID"]
        All_Proteins = ";".join(Proteins_containing_peptide)
        return {'peptide':peptide,'All_Proteins':All_Proteins}
    
    def find_matches(self,to_process):
        all_peptides={}
        for peptide in to_process:
            results = self.retrieve_all_proteins(peptide,self.Reference_Proteome)
            all_peptides[results['peptide']]=results['All_Proteins']
        return all_peptides
        
    def determine_all_proteins(self):
        import multiprocessing as mp  
        pool = mp.Pool(self.cpus)
        count= 0
        # Here we split the peptides in batches
        to_process=[]
        for peptide in self.peptides:
            count+=1
            if (count % 1000 == 0):
                if self.cpus>1:
                    pool.apply_async(self.find_matches,args=([to_process]),callback=self.append_results)
                else:
                    result = self.find_matches(to_process)
                    self.append_results(result)
                to_process=[]
            else:
                # Here we append to the list to process
                to_process.append(peptide)
        # Here we process the remaining ones that dint get triggered by exacution 
        if self.cpus>1:
            pool.apply_async(self.find_matches,args=([to_process]),callback=self.append_results)
        else:
            result = self.find_matches(to_process)
            self.append_results(result)
        pool.close()
        pool.join() 
        return self.protein_peptides
    
def match_peptide_to_protein(Protein_peptides,Reference_Proteome,cpus=1):
    # Here we determine all the protrein identities if not already provided in file.
    print("Assigning protein IDs based on search")
    Protein_peptides["Protein"]=""
    # sometimes data outputs contain a peptide sequence in brackets - the folowing will remove this
    Protein_peptides.Peptide=Protein_peptides.Peptide.str.replace(".*\]\.",'')
    Protein_peptides.Peptide=Protein_peptides.Peptide.str.replace("\.\[.*","")
    peptide_protein_identities = Determine_Peptide_Protein_Identities(Protein_peptides.Peptide.unique(),Reference_Proteome,cpus=cpus).determine_all_proteins()
    for key in peptide_protein_identities.keys():
        Protein_peptides.loc[Protein_peptides.Peptide == key,"Protein"]=peptide_protein_identities[key]
    return Protein_peptides

def retrieve_mysql_data(id_to_process,cpus=1):

    # import requests
    # response = requests.get("http://www.manchesterproteome.manchester.ac.uk/run_api/MSP_api/?page=1&search=")
    # print(response.json())
    # d=response.json()
    # with open("test_output.json", 'w') as json_file:
    #     json.dump(d, json_file)


    
    from secret import HOST, PORT, PASSWORD, DB, USER
    connection = mysql.connector.connect(host=HOST,
                                        database=DB,
                                        user=USER,port=PORT,
                                        password=PASSWORD,
                                        auth_plugin='mysql_native_password')
    cursor_query = connection.cursor()

    # # here could select the jobs that are qued - do this every 6h and if a new job is qued then process on the HPC cloud
    # sql="SELECT id,name FROM `Structural_userdata` WHERE Progress LIKE 'Que'"
    # cursor = connection.cursor()
    # cursor.execute(sql)
    # Data_ids = pd.DataFrame(cursor.fetchall())
    # if not Data_ids.empty:
    #     field_names = [i[0] for i in cursor.description]
    #     Data_ids.columns = field_names
    # id_to_process="138"
    sql=f"SELECT * FROM `Structural_userdata` WHERE `id` LIKE '{id_to_process}'"
    cursor = connection.cursor()
    cursor.execute(sql)
    Data = pd.DataFrame(cursor.fetchall())
    if not Data.empty:
        field_names = [i[0] for i in cursor.description]
        Data.columns = field_names
    print(f"id = {id_to_process}")

    Owner_ID = Data.owner_id[0]
    
    Domain_types=json.loads(Data.Domain_Types[0])
    # Domain_types=['DOMAINS', 'REGIONS', 'TOPO_DOM', 'TRANSMEM', 'REPEAT', '50.0 AA STEP']
    experiment_feed = json.loads(Data.experimental_design[0])
    
    id = Data.id[0]
    Data_Val = Data.data[0]
    Data_Val = json.loads(Data_Val)
    Protein_peptides = pd.DataFrame(Data_Val)
    Paired= Data.Paired[0]
    Spiecies=Data.Spiecies[0]
    print(f"Spiecies: {Spiecies}")
    # Domain_types = ['DOMAINS', 'REGIONS', 'TOPO_DOM', 'TRANSMEM', 'REPEAT', '50.0 AA STEP']
    # Protein_peptides=pd.read_csv('Sample_Data/Protein_peptides.tsv',sep='\t',index_col=[0])
    # Protein_peptides.to_csv('Sample_Data/Protein_peptides.tsv',sep='\t')
    
    # experiment_feed = pandas_to_experiment(pd.read_csv('Sample_Data/Experiment_feed.tsv',sep='\t',index_col=[0]))
    # experiment_feed_pd = pd.DataFrame(experiment_feed)
    # experiment_feed_pd.to_csv('Sample_Data/Experiment_feed.tsv',sep='\t')

    run_full_analysis(Domain_types, Protein_peptides, experiment_feed,Owner_ID=Owner_ID,cpus=cpus,id=id,paired=Paired, Spiecies=Spiecies)

    # '''
    # with open(f"bin/Structural_Json_{Spiecies}_{Owner_ID}_{id}.json", 'r') as myfile:
    #     Structural_Json=myfile.read()
    # Structural_Json=json.loads(Structural_Json)
    # record_data(Structural_Json, Owner_ID,id,Domain_types)
    # '''

def pandas_to_experiment(df):
    dict={}
    dict[df.iloc[:,0].name]=list(df.iloc[:,0])
    dict[df.iloc[:,1].name]=list(df.iloc[:,1])
    return dict

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
    
    Data_ids.to_csv("tmp_ids.csv")

def retrieve_mysql_data_test():
    import requests
    response = requests.get('http://www.manchesterproteome.manchester.ac.uk/run_api/MSP_api/?page=1&search=')
    d = response.json()
    with open("test_output.json", 'w') as json_file:
        json.dump(d, json_file)

    print("Done")
        
def update_user_id():
    from secret import HOST, PORT, PASSWORD, DB, USER
    connection = mysql.connector.connect(host=HOST,
                                        database=DB,
                                        user=USER,port=PORT,
                                        password=PASSWORD,
                                        auth_plugin='mysql_native_password')
    cursor_query = connection.cursor()
    sql="SELECT id,owner_id,name FROM `Structural_userdata` WHERE 1"
    
    cursor = connection.cursor()
    cursor.execute(sql)
    Data_ids = pd.DataFrame(cursor.fetchall())
    if not Data_ids.empty:
        field_names = [i[0] for i in cursor.description]
        Data_ids.columns = field_names
    # print(Data_ids)
    id=138
    sql=f"UPDATE `Structural_userdata` SET `owner_id`=11 WHERE id LIKE {id}"
    cursor = connection.cursor()
    cursor.execute(sql)
    
def local_run():
    import os
    Spiecies='HUMAN'
    Paired=1
    Domain_types = ['DOMAINS', 'REGIONS', 'TOPO_DOM', 'TRANSMEM', 'REPEAT', '50.0 AA STEP']
    cpus = os.cpu_count()
    cpus = 1
    try:
        Protein_peptides=pd.read_csv('Sample_Data/sample_inputs_small/Sample_Data_For_Analysis.csv',index_col=False)
    except:
        Protein_peptides=pd.read_csv('Sample_Data/sample_inputs_small/Sample_Data_For_Analysis.csv',sep='\t',index_col=False)
    #  If more then 1 spectra columns, then add them up.
    matching = [s for s in Protein_peptides.columns if "spectra" in s]
    if len(matching)>1:
        summed_spectra = Protein_peptides[matching].sum(axis=1)
        Protein_peptides.drop(matching,axis=1,inplace=True)
        Protein_peptides['spectra']=summed_spectra
        del summed_spectra
    
    experiment_feed = pandas_to_experiment(pd.read_csv('Sample_Data/sample_inputs_small/Experiment_feed.tsv',sep='\t',index_col=False))
    run_full_analysis(Domain_types, Protein_peptides, experiment_feed,cpus=cpus,paired=Paired, Spiecies=Spiecies)

if __name__ == '__main__':
    # '''bsub -o exercise5.output -n10 -R"select[mem>2500] rusage[mem=2500]" -M2500 python MS_Total_Software.py '''
    # export  LSB_DEFAULTGROUP=hgi
    # id_to_process = sys.argv[1]
    # cpus = sys.argv[2]
    # print(cpus)
    # cpus=int(cpus)
    # print(f"using {cpus} cpus")

    # # update_user_id()
    # # retrieve_mysql_data_test()
    # # retrieve_save_and_process()
    # # retrieve_save_and_process()
    # retrieve_mysql_data(id_to_process,cpus=cpus)
    
    local_run()


    # import json
    # paired=1
    # id='1'
    # Owner_ID='1'
    # with open("sample_input/experiment_feed.json", 'r') as myfile:
    #     experiment_feed=myfile.read()
    # experiment_feed=json.loads(experiment_feed)
    # Protein_peptides=pd.read_csv("sample_input/Protein_peptides.tsv",sep="\t",index_col=0)
    # Domain_types=pd.read_csv("sample_input/Domain_types.tsv",sep="\t",index_col=0,names=["Dom"],header=1)
    # Domain_types=Domain_types.Dom.values.tolist()
    # Protein_peptides=pd.read_csv("tmp_working_file.csv",sep="\t",index_col=0)
    # run_full_analysis(Domain_types, Protein_peptides, experiment_feed,Owner_ID,id,paired=paired)
