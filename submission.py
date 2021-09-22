import pandas as pd
import sys
import mysql.connector
from Functions_Clean import retrieve_reviewed, Master_Run_Counting_Algorythm_Clean, \
    Master_Run_Structural_Analysis, Master_Run_Score_Calculations
import json
from MS_Total_Software import retrieve_mysql_data

def start_tasks():
    from secret import HOST, PORT, PASSWORD, DB, USER
    connection = mysql.connector.connect(host=HOST,
                                        database=DB,
                                        user=USER,port=PORT,
                                        password=PASSWORD,
                                        auth_plugin='mysql_native_password')
    cursor_query = connection.cursor()
    cpus=4
    # here could select the jobs that are qued - do this every 6h and if a new job is qued then process on the HPC cloud
    sql="SELECT id,name FROM `Structural_userdata` WHERE name LIKE 'mouse aorta - 22 vs 72 w all'"

    cursor = connection.cursor()
    cursor.execute(sql)
    Data_ids = pd.DataFrame(cursor.fetchall())
    if not Data_ids.empty:
        field_names = [i[0] for i in cursor.description]
        Data_ids.columns = field_names
    print(Data_ids)
    import os

    count=0
    for id in Data_ids.id:
        os.system(f"qsub run_on_HPC.sh {id} {cpus}")
        print(id)
        # retrieve_mysql_data(id,cpus=8)
        # if count>0:
        #     break
        count+=1
        


if __name__ == '__main__':
    '''bsub -o exercise5.output -n10 -R"select[mem>2500] rusage[mem=2500]" -M2500 python MS_Total_Software.py '''
    # export  LSB_DEFAULTGROUP=hgi
    
    start_tasks()
    # retrieve_mysql_data_test()
    # retrieve_save_and_process()
    # retrieve_mysql_data()
