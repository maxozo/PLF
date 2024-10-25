#!/usr/bin/env python
__date__ = '2023-09-13'
__version__ = '0.0.1'

import pandas as pd
import json
import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

Experimental_coverages_all = {}
Reference_Proteome=None
Reference_Domains=None
Protein_peptides2={}

class PLF:
    
    ######################
    # This is a main PLF analysis processing class.
    ######################
    
    def __init__(self,Protein_peptides,experiment_feed,Spiecies='HUMAN',Domain_Types=['DOMAINS', 'REGIONS', 'TOPO_DOM', 'TRANSMEM', 'REPEAT', '50.0 AA STEP'],paired=False,cpus=1,p_threshold=0.05,protein_list=None):
        self.Spiecies = Spiecies
        self.Domain_types = Domain_Types
        self.Protein_peptides=Protein_peptides
        self.Coverage_Json={}
        self.experiment_feed=experiment_feed
        self.Structural_Json = {}
        self.cpus = cpus
        self.paired=paired
        self.Full_MPLF_Results = pd.DataFrame()
        self.p_threshold = p_threshold
        self.protein_list = protein_list
        

    def MPLF(self,Protein,Reference_Proteome,Reference_Domains,Protein_Entries,count,total_count):

        ###########
        ## First all the coverages are calculated for each of the experiments within domains
        ## Then the statistical analysis is performed to determine significant hits
        ###########

        from functions.Functions_Clean import MPLF_Domain_Quantifications, MPLF_Statistical_Analyisis

        Domain_types = self.Domain_types
        experiment_feed = self.experiment_feed
        paired = self.paired

        Data_entries = Reference_Proteome[Reference_Proteome.Uniprot_ID==Protein]
        if len(Data_entries)==0:
            Data_entries = Reference_Proteome[Reference_Proteome.Uniprot_AC.str.contains(f"^{Protein}|\\n{Protein}",regex=True)]
        Uniprot_ID = Data_entries.Uniprot_ID.values[0].replace('\n','')
        Uniprot_AC = Data_entries.Uniprot_AC.values[0].replace('\n','')

        print(f"Analysing protein: : {Protein} | {count}/{total_count}")
        # Each of the domain coverages per experiment is calculated in the folowing code
        Experiment_Coverages = MPLF_Domain_Quantifications(Protein=Protein,
                                                                        Domain_Types=Domain_types,
                                                                        Protein_Entries=Protein_Entries,Reference_Proteome=Reference_Proteome,
                                                                        Reference_Domains=Reference_Domains)
        # Once the quantification has been performed then we perform the statistical analysis based on experimental design provided.
        Structural_Analysis_Results, Norm_Factors = MPLF_Statistical_Analyisis(experiment_feed=experiment_feed,
                                                                                Results=Experiment_Coverages,
                                                                                Protein=Protein, paired=paired,cuttoff_p_val=self.p_threshold)
        Structural_Analysis_Results2 = Structural_Analysis_Results.copy()
        Structural_Analysis_Results2['Uniprot_ID']=Uniprot_ID
        Structural_Analysis_Results2['Uniprot_ACs']=Uniprot_AC
        Structural_Analysis_Results.drop("GeneAC", axis=1, inplace=True)
        Experiment_Coverages.drop("GeneAC", axis=1, inplace=True)
        Experiment_Coverages = Experiment_Coverages[
            Experiment_Coverages["Domain Type"].isin(Structural_Analysis_Results["Domain Type"].unique())]
        Experiment_Coverages.set_index("Domain_Name", drop=False, inplace=True)
        Coverage = Experiment_Coverages.to_dict()
        Structural_Json= {"Data": Structural_Analysis_Results.to_dict('index'),
                                    "Norm_Factors": Norm_Factors.to_dict()}
        
        return (Coverage,Structural_Json,Protein,Structural_Analysis_Results2)

    def append_protein_to_dictionary(self,Protein_peptides,Reference_Proteome):
        
        #######################
        # This function determines how many isoforms of the same gene is analysed. 
        # For the performence and downstream analyisis purpose the pipeline is run 
        # with only one of the isoforms. If uniprot reviewed entry is availavle then 
        # this will be utilised, otherwise a random splicing version will be analysed.
        #######################
        
        All_proteins={}
        # Protein_peptides[Protein_peptides['Protein'].str.contains(";", na=False)]
        for i, Protein in enumerate(Protein_peptides['Protein'].str.split(";").explode().unique()):
            # Here gather all the unique gene names - all the revirewed entries + each unique non uniprot entry
            if Protein!=Protein:
                continue
            if self.protein_list:
                if Protein not in self.protein_list:
                    continue
            Prot1=Reference_Proteome[Reference_Proteome["Uniprot_ID"]==Protein]
            if len(Prot1)==0:
                Prot1=Reference_Proteome[Reference_Proteome["Uniprot_AC"].str.contains(f"^{Protein}|\\n{Protein}",regex=True)]
                if len(Prot1)>1:
                    Prot1=Prot1[Prot1["Uniprot_Type"]=='Uniprot']
                    if len(Prot1)==0:
                        Prot1=Prot1.iloc[0]
            
            try:
                Gene=Prot1["Uniprot_Gene"].values[0].split(" {")[0]
                Gene=Gene.split(' ORFN')[0]
            except:
                # Gene name is not listed in the reference proteome, hence proceeding with the ID represented in the peptide file
                Gene=Protein
                
            try:
                Type=Prot1["Uniprot_Type"][0]
            except:
                # If the type is not available w asume that this is an unreviewed entity.
                Type='Trembl'    
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
        return All_proteins

        
    def collect_result(self,result):
        Protein1=result[2]
        Coverage_Json1=result[0]
        Structural_Json1=result[1]
        self.Coverage_Json[Protein1] = Coverage_Json1
        self.Structural_Json[Protein1] = Structural_Json1
        self.Full_MPLF_Results = pd.concat([self.Full_MPLF_Results,result[3]])
           
    def PLF_Analysis(self):
        i=0
        import multiprocessing as mp
        import time
        Reference_Proteome = pd.read_csv(f"{dir_path}/outputs/Uniprot_{self.Spiecies}.tsv",sep="\t",index_col=0)
        Reference_Domains = pd.read_csv(f"{dir_path}/outputs/Domains_Uniprot_{self.Spiecies}.tsv",sep="\t",index_col=0)

        print('Prepearing file for analysis >>>')
        t = time.time()
        Protein_isoform_grouping = self.append_protein_to_dictionary(self.Protein_peptides,Reference_Proteome)
        Nr_proteins = len(Protein_isoform_grouping.keys())
        print(f"Analysing {Nr_proteins} proteins >>>")
        pool = mp.Pool(self.cpus)
        proteins = set(Protein_isoform_grouping.keys())
        try:
            import numpy as np
            proteins.remove(np.nan)
        except:
            _='no nan'
        count=0
        total_count=len(proteins)
        for key in proteins:
            
            count+=1
            # if count>40:
            #     break
            try:
                Protein= Protein_isoform_grouping[key]['Uniprot'][0]
            except:
                Protein= Protein_isoform_grouping[key]['Trembl'][0]
            # if (Protein!='FBLN1_HUMAN'):
            #     continue              
            Protein_Entries = self.Protein_peptides[self.Protein_peptides['Protein'].str.contains(Protein, na=False)]
            
            try:
                Fasta = Reference_Proteome[Reference_Proteome.Uniprot_ID==Protein]["FASTA"]
                if len(Fasta)==0:
                    Fasta = Reference_Proteome[Reference_Proteome.Uniprot_AC.str.contains(f"^{Protein}|\\n{Protein}|{Protein}.*",regex=True)]["FASTA"]
                Fasta = Fasta.values[0]
            except:
                print(f"Protein {Protein}, doesnt have a fasta reference in Uniprot version utilised")
                print(f"Skipping {Protein} analysis.")
                continue
            if self.cpus>1:
                pool.apply_async(self.MPLF, args=([Protein,Reference_Proteome,Reference_Domains,Protein_Entries,count,total_count]),callback=self.collect_result) #paralel runs - uses all the cores available
            else:

                result = self.MPLF(Protein,Reference_Proteome,Reference_Domains,Protein_Entries,count,total_count)
                self.collect_result(result)
            i+=1
        pool.close()
        pool.join()
        elapsed = time.time() - t
        print(f'execution took {elapsed} seconds')
        return (self.Coverage_Json,self.Structural_Json,self.Full_MPLF_Results)      

   
def run_full_analysis( Domain_types, Protein_peptides, experiment_feed, cpus=1,paired=False, Spiecies="HUMAN",outname='Default_MPLF',p_threshold=0.05,protein_list=None):
    # If we do decide to remove the protein entry then we have to look up each peptide in the library and find all the peptides for the protein thatr are provided.

    if not 'Protein' in list(Protein_peptides.columns):
        Protein_peptides=match_peptide_to_protein(Protein_peptides,Reference_Proteome,cpus=cpus)
    elif  Protein_peptides.Protein.unique()[0]=='undefined':
        Protein_peptides=match_peptide_to_protein(Protein_peptides,Reference_Proteome,cpus=cpus)

    Protein_peptides=Protein_peptides.dropna(subset=['spectra']) # Drop the NaN vales on spectra. ie - peptides are not detected in that sample
    Protein_peptides.Protein = Protein_peptides.Protein.str.replace(',',';').str.replace(' ',';')
    print(f'For analysis we are using: {cpus} cpus')
    print(f'Spiecies has been set to: {Spiecies}')
    ##################
    # PLF processing of each of the proteins
    ##################
    Coverage_Json,Structural_Json,Full_Results_TSV = PLF(Protein_peptides,experiment_feed,Spiecies=Spiecies,paired=paired,Domain_Types=Domain_types,cpus=cpus,p_threshold=p_threshold,protein_list=protein_list).PLF_Analysis()
    Full_Results_TSV.to_csv(f'{outname}.tsv',sep='\t',index=False)

    mplf_export_data={}
    mplf_export_data['Structural_Results']=Structural_Json
    mplf_export_data['experiment_feed']=experiment_feed
    mplf_export_data['Domain_types']=Domain_types
    mplf_export_data['Significant_Protein_Count']=len(Structural_Json.keys())
    mplf_export_data['Paired']=paired
    mplf_export_data['Spiecies']=Spiecies
    with open(f'{outname}.mplf', 'w') as f:
        json.dump(mplf_export_data, f,allow_nan=False,)
    print('Analysis completed ... ') 
    print(f'Check your results at {outname}.tsv')
    return "success"
    
def retrieve_all_proteins(peptide,Reference_Proteome):
    Proteins_containing_peptide = Reference_Proteome[Reference_Proteome.FASTA.str.contains(peptide)]["Uniprot_ID"]
    All_Proteins = ";".join(Proteins_containing_peptide)
    return {'peptide':peptide,'All_Proteins':All_Proteins}


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


def pandas_to_experiment(df):
    dict={}
    dict[df.iloc[:,0].name]=list(df.iloc[:,0])
    dict[df.iloc[:,1].name]=list(df.iloc[:,1])
    for k1 in dict.keys():
        dict[k1] = [x for x in dict[k1] if str(x) != 'nan']
    return dict

def PLF_run(options):
    print(f'Analysing file {options.peptides} with experimental design file {options.experimental_design}.......')
    Spiecies=options.spiecies
    Paired= not options.paired.lower()=='false'
    p_threshold=float(options.p_threshold)
    Domain_types_pre = set(options.domain_types.split(','))
    Domain_types=[]
    for s in Domain_types_pre:
        if "AA" in s:
            s=s.replace('AA','.0 AA STEP')
        Domain_types.append(s)
    cpus = options.cpus
    outname=options.outname
    if cpus=='max':
        cpus=os.cpu_count()
    else:
        cpus=int(cpus)
    try:
        Protein_peptides=pd.read_csv(options.peptides,index_col=False)
    except:
        Protein_peptides=pd.read_csv(options.peptides,sep='\t',index_col=False)
    
    #  If more then 1 spectra columns, then add them up.
    matching = [s for s in Protein_peptides.columns if "spectra" in s]
    if len(matching)>1:
        summed_spectra = Protein_peptides[matching].sum(axis=1)
        Protein_peptides.drop(matching,axis=1,inplace=True)
        Protein_peptides['spectra']=summed_spectra
        del summed_spectra
    

    experiment_feed = pandas_to_experiment(pd.read_csv(options.experimental_design,sep='\t',index_col=False))
    run_full_analysis(Domain_types, Protein_peptides, experiment_feed,cpus=cpus,paired=Paired, Spiecies=Spiecies,outname=outname,p_threshold=p_threshold,protein_list=options.protein_list)

     
def test_run(options):
    print('Test Run.......')
    dir1= os.path.dirname(os.path.abspath(sys.argv[0]))
    Spiecies=options.spiecies
    Paired=not options.paired.lower()=='false'
    p_threshold=options.p_threshold
    Domain_types_pre = set(options.domain_types.split(','))
    Domain_types=[]
    for s in Domain_types_pre:
        if "AA" in s:
            s=s.replace('AA','.0 AA STEP')
        Domain_types.append(s)
    cpus = options.cpus
    outname=options.outname
    if cpus=='max':
        cpus=os.cpu_count()
    else:
        cpus=int(cpus)
    try:
        Protein_peptides=pd.read_csv(f'{dir1}/sample_data/sample_inputs_small/Sample_Data_For_Analysis.csv',index_col=False)
    except:
        Protein_peptides=pd.read_csv(f'{dir1}/sample_data/sample_inputs_small/Sample_Data_For_Analysis.csv',sep='\t',index_col=False)
    
    #  If more then 1 spectra columns, then add them up.
    matching = [s for s in Protein_peptides.columns if "spectra" in s]
    if len(matching)>1:
        summed_spectra = Protein_peptides[matching].sum(axis=1)
        Protein_peptides.drop(matching,axis=1,inplace=True)
        Protein_peptides['spectra']=summed_spectra
        del summed_spectra
    
    experiment_feed = pandas_to_experiment(pd.read_csv(f'{dir1}/sample_data/sample_inputs_small/Experiment_feed.tsv',sep='\t',index_col=False))
    run_full_analysis(Domain_types, Protein_peptides, experiment_feed,cpus=cpus,paired=Paired, Spiecies=Spiecies,outname=outname,p_threshold=p_threshold)

def cli():
    parser = argparse.ArgumentParser(
        description="""
            Combines all the Raw cellranger metadata to be passed in the h5ad
            """
    )  
    parser.add_argument(
        '--experimental_design',
        action='store',
        dest='experimental_design',
        required=False,
        default=None,
        help='Path')
    
    parser.add_argument(
        '--peptides',
        action='store',
        dest='peptides',
        required=False,
        default=None,
        help='Path')
    
    parser.add_argument(
        '--spiecies',
        action='store',
        dest='spiecies',
        required=False,
        default='HUMAN',
        help='spiecies')
    
    parser.add_argument(
        '--domain_types',
        action='store',
        dest='domain_types',
        required=False,
        default='DOMAINS,REGIONS,TOPO_DOM,TRANSMEM,REPEAT,50AA,100AA',
        help='domain_types')
    
    parser.add_argument(
        '--paired',
        action='store',
        dest='paired',
        required=False,
        default='False',
        help='paired')
    
    parser.add_argument(
        '--cpus',
        action='store',
        dest='cpus',
        required=False,
        default='max',
        help='cpus')
    
    parser.add_argument(
        '--test',
        action='store_true',
        dest='test',
        help='test run')

    parser.add_argument(
        '--outname',
        action='store',
        dest='outname',
        required=False,
        default='MPLF_Results',
        help='outname')
    
    parser.add_argument(
        '--p_threshold',
        action='store',
        dest='p_threshold',
        required=False,
        default=0.05,
        help='p_threshold')
    
    parser.add_argument(
        '--protein_list',
        action='store',
        dest='protein_list',
        required=False,
        default=None,
        help='protein_list')

    options = parser.parse_args()
    return options


if __name__ == '__main__':
    # PLF main processing 

    options=cli()
    
    print(options.test)
    if options.test:
        test_run(options)
    else:
        PLF_run(options)
