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

def collect_result(result):
    
    Protein1=result[2]
    Coverage_Json1=result[0]
    Structural_Json1=result[1]
    
    Coverage_Json[Protein1] = Coverage_Json1

    Structural_Json[Protein1] = Structural_Json1

def record_data(Structural_Json,Coverage_Json,Owner_ID,id):
    import mysql.connector
    from MSP.settings import PASSWORD, USER, HOST, DB
    connection_db = mysql.connector.connect(host=HOST,
                                        database='MSP',
                                        user=USER,
                                        password=PASSWORD,
                                        auth_plugin='mysql_native_password')

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

def run_protein(Protein,Reference_Proteome,Reference_Domains):
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
        return (0,0,0)
        
def run_full_analysis( Domain_types, Protein_peptides, experiment_feed, Owner_ID, id, paired=False):
    # If we do decide to remove the protein entry then we have to look up each peptide in the library and find all the peptides for the protein thatr are provided.
    import multiprocessing as mp
    # pool.close()
    pool = mp.Pool(mp.cpu_count())
    # pool = mp.Pool(1)
    Reference_Proteome = pd.read_csv("/home/matiss/work/expansion/MPLF/outputs/Uniprot_HUMAN.tsv",sep="\t",index_col=0)
    Reference_Domains = pd.read_csv("/home/matiss/work/expansion/MPLF/outputs/Domains_Uniprot_HUMAN.tsv",sep="\t",index_col=0)

    k_val=0
    # paired=True


    for i, Protein in enumerate(Protein_peptides.Protein.str.split(";").explode().unique()):


        try:
            pool.apply_async(run_protein, args=([Protein,Reference_Proteome,Reference_Domains]),callback=collect_result) #paralel runs - uses all the cores available
            # run_protein(Protein,Reference_Proteome,Reference_Domains) #individual cores - uses 1 core hence 1 protein at a time is analysed.
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

    with open("./bin/Structural_Json.json", 'r') as myfile:
        Structural_Json2=myfile.read()
    Structural_Json2=json.loads(Structural_Json2)
    # here we wait for the process to finish and then HPC should send the data back.


    # # this is to record last entries that are not processed
    # if (Structural_Json.keys().__len__()>0):
    #     Structural_Json, Coverage_Json=record_data(Structural_Json, Coverage_Json, Owner_ID,id)

    # connection_db.disconnect()
    # cursor_query.close()
    return "success"

def match_peptide_to_protein(Protein_peptides):
    
    # limit to the protein of interest.
    Reference_Proteome = Reference_Proteome[Reference_Proteome.Uniprot_Type=="Uniprot"]


    Protein_peptides["Protein"]=""
    count=0
    for peptide in Protein_peptides.Peptide.unique():
        count+=1
        # print(count)
        Proteins_containing_peptide = Reference_Proteome[Reference_Proteome.FASTA.str.contains(peptide)]["Uniprot_ID"]
        All_Proteins = ";".join(Proteins_containing_peptide)
        Protein_peptides.loc[Protein_peptides.Peptide == peptide,"Protein"]=All_Proteins

    return Protein_peptides


if __name__ == '__main__':

    import json
    paired=1
    id='1'
    Owner_ID='1'
    with open("/home/matiss/work/expansion/MPLF/sample_input/experiment_feed.json", 'r') as myfile:
        experiment_feed=myfile.read()
    experiment_feed=json.loads(experiment_feed)
    Protein_peptides=pd.read_csv("/home/matiss/work/expansion/MPLF/sample_input/Protein_peptides.tsv",sep="\t",index_col=0)
    Domain_types=pd.read_csv("/home/matiss/work/expansion/MPLF/sample_input/Domain_types.tsv",sep="\t",index_col=0,names=["Dom"],header=1)
    Domain_types=Domain_types.Dom.values.tolist()
    Protein_peptides=pd.read_csv("tmp_working_file.csv",sep="\t",index_col=0)
    run_full_analysis(Domain_types, Protein_peptides, experiment_feed,Owner_ID,id,paired=paired)
