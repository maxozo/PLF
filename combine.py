import pandas as pd

Data1=pd.read_csv("/home/matiss/work/MPLF/bin/Uniprot_Entries_HUMAN.csv",sep="\t")
Data1 = Data1.drop(['Min_Ft', 'Max_Ft','Main_AC'], axis=1)
Data2=pd.read_csv("/home/matiss/work/MPLF/bin/Trembl_Entries_HUMAN.csv",sep="\t")
print(_)