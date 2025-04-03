import pandas as pd
import json
import argparse
import os
import sys
import warnings
dir_path = f"{os.path.dirname(os.path.realpath(__file__))}/../"


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

        from .Functions_Clean import MPLF_Domain_Quantifications, MPLF_Statistical_Analyisis

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
        if Structural_Analysis_Results is None:
            return None
        Structural_Analysis_Results2 = Structural_Analysis_Results.copy()
        Structural_Analysis_Results2['Uniprot_ID']=Uniprot_ID
        Structural_Analysis_Results2['Uniprot_ACs']=Uniprot_AC
        try:
            Structural_Analysis_Results.drop("GeneAC", axis=1, inplace=True)
        except:
            print('no gene ac')
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
        if result is not None:
        
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