#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-07-28'
__version__ = '0.0.1'

from os import listdir
import os
import pandas as pd

if __name__ == "__main__":

    print('Lets convert the excel files to a csv file')

    dir='/Volumes/GoogleDrive-110351538372930204855/My Drive/HK_Analysis'
    Folders = listdir(dir)
    try:
        os.mkdir(f'{dir}/csv')
    except:
        print("dir exists") 


    for file in Folders:
        print(file)
        df =  pd.read_excel(f"{dir}/{file}")
        read_from = df[df.iloc[:,0]=='Experiment name'].index[0]
        df =  pd.read_excel(f"{dir}/{file}",skiprows=read_from+1)
        df.to_csv(f"{dir}/csv/{file.split('.')[0]} Mouse IVDs.tsv",sep='\t')
        
