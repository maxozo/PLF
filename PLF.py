#!/usr/bin/env python
__date__ = '2025-04-03'
__version__ = '0.0.1'

import pandas as pd
import json
import argparse
import os
import sys
import warnings
from functions.PLF import PLF
from functions.Determine_Peptide_Protein_Identities import match_peptide_to_protein
warnings.filterwarnings("ignore")

import os 

Experimental_coverages_all = {}
Reference_Proteome=None
Reference_Domains=None
Protein_peptides2={}

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
