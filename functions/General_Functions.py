import matplotlib as mpl
import pandas as pd
import os
import urllib
import urllib.request  as urllib2
import mysql.connector

def find_str(s, char):
    # from http://stackoverflow.com/questions/21842885/python-find-a-substring-in-a-string-and-returning-the-index-of-the-substring
    index = 0

    if char in s:
        c = char[0]
        for ch in s:
            if ch == c:
                if s[index:index + len(char)] == char:
                    return index

            index += 1

    return -1

def grep_if_is_in(x, y):
    # print 'grep_if_is_in'
    import re
    count = 0
    array = list()

    for line in x:
        line = str(line)
        y = str(y)
        #y=re.compile(y)
        if bool(re.search(y,line)):
            # if y == line:
            array.append(count)

        count = count + 1
    return array

def connect_With_Database(mycursor):
    
    mydb = mysql.connector.connect(
        host="localhost",
        user="root",
        passwd="Weafrae1")

    mycursor = mydb.cursor(buffered=True)

    try:
        # create a new database
        mycursor.execute("CREATE DATABASE analysis;")
    except:
        print("exists")
    mycursor.execute("use analysis;")
    return mycursor

def mysql_record_experiment_mappings(mapping=None, database=None,table=None):

    '''This records the experiment name mappings'''

    
    mydb = mysql.connector.connect(
        host="localhost",
        user="root",
        passwd="Weafrae1")
    mycursor = mydb.cursor(buffered=True)

    try:
        # create a new database
        mycursor.execute("CREATE DATABASE `"+database+"`;")
    except:
        print("exists")

    mycursor.execute("use `"+database+"`;")

    try:
        mycursor.execute(
            " CREATE TABLE   "+table+" (id INT(11) NOT NULL AUTO_INCREMENT, Old_exp_name VARCHAR(200), New_exp_name VARCHAR(200),    PRIMARY KEY (id))")
        for ind in range(0, mapping.__len__(), 1):
            New_exp_name = mapping['exp_name'][ind]
            Old_exp_name = mapping['t'][ind]
            querry="INSERT INTO "+table+" (id, Old_exp_name, New_exp_name) VALUES (NULL,'" + Old_exp_name + "','" + New_exp_name + "')"
            mycursor.execute(querry)

    except:
        print("table exists")

    mydb.commit()

    mycursor.close()
    mydb.close()

def Get_Domains(data=None):
    import re
    Domain_bound_index1 = grep_if_is_in(data, "FT   DOMAIN   ")
    Repeat_bound_index = grep_if_is_in(data, "FT   REPEAT ")
    Repeat_bound_index2 = grep_if_is_in(data, "FT   TOPO_DOM ")
    Repeat_bound_index3 = grep_if_is_in(data, "FT   REGION ")
    Domain_bound_index_group = [Domain_bound_index1,Repeat_bound_index,Repeat_bound_index2,Repeat_bound_index3]

    #FT   REGION
    #Repeat_bound_index3 = grep_if_is_in(data, "FT   REGION ")
    Domain_bound_index=Domain_bound_index1+Repeat_bound_index+Repeat_bound_index2#+Repeat_bound_index3#+Repeat_bound_index3
    df_with_doamin_info = pd.DataFrame(columns=['Domain', 'Domain Start', 'Domain Finish'])
    for Domain_bound_index in Domain_bound_index_group:
        for element in Domain_bound_index:
            line = re.split(r'\s\s+', data[element])
            Domain_name = line[4].split(".")[0]
            Domain_name = Domain_name.split(";")[0]
            Domain_name=Domain_name.replace("\n","")
            Domain_name = Domain_name.replace("\r", "")
            Domain_type = line[1]
            Domain_start = line[2]
            Domain_finish = line[3]

            if (Domain_name in df_with_doamin_info.Domain.tolist()):
                Domain_name=Domain_name+"_p_"+str(Domain_start)+"_"+str(Domain_finish)


            df_with_doamin_info = df_with_doamin_info.append(
                {'Domain': Domain_name, 'Domain Start': Domain_start, 'Domain Finish': Domain_finish,'Domain Type':Domain_type},
                ignore_index=True)
    return df_with_doamin_info

def mysql_record_results_Experiment_analysis(experiment_name=None,susceptibilities_of_domains_for_UVA=None,database=None,protein=None,Exclusive_spectrum_count=None):

    '''This records all experiment entries for the domain analysis'''
    print(experiment_name)
    
    mydb = mysql.connector.connect(
        host="localhost",
        user="root",
        passwd="Weafrae1")
    mycursor = mydb.cursor(buffered=True)

    try:
        # create a new database
        mycursor.execute("CREATE DATABASE `"+database+"`;")
    except:
        print("database exists")

    mycursor.execute("use `"+database+"`;")

    #First select a correct mapping name

    querry = " CREATE TABLE Analysis (id INT(11) NOT NULL AUTO_INCREMENT, GeneAC VARCHAR(100),experiment_name VARCHAR(5000), Domain_Name VARCHAR(100), Domain_Start INT(11), Domain_Finish INT(11), NumberOfSpectra INT(11), Colours VARCHAR(30),Percent_Covered INT(11), Exclusive_spectrum_count INT(6), peptides_found VARCHAR(10000), PRIMARY KEY (id))"
    try:
        mycursor.execute(querry)
    except:
       pass


    #Loop through each dataset and put the data information in database
    for ind in range(0,susceptibilities_of_domains_for_UVA.__len__(),1):
        Domain_name=susceptibilities_of_domains_for_UVA['Domain Name'][ind]
        Domain_Start=susceptibilities_of_domains_for_UVA['Domain Start'][ind]
        Domain_Finish=susceptibilities_of_domains_for_UVA['Domain Finish'][ind]
        NumberOfSpectra=susceptibilities_of_domains_for_UVA['NumberOfSpectra'][ind]
        Colours=susceptibilities_of_domains_for_UVA['Colours'][ind]
        Colours=tuple(map(lambda x: isinstance(x, float) and round(x, 2) or x, Colours))
        Percent_Covered=susceptibilities_of_domains_for_UVA['percentage_covered'][ind]
        peptides_found=susceptibilities_of_domains_for_UVA['peptides_found'][ind]
        str_peptides_found =','.join(peptides_found)
        querry="INSERT INTO Analysis (`id`, `GeneAC`,`experiment_name`,`Domain_Name`, `Domain_Start`, `Domain_Finish`, `NumberOfSpectra`, `Colours`,`Percent_Covered`,`Exclusive_spectrum_count`,`peptides_found`) VALUES (NULL, '"+ str(protein)+"','"+ str(experiment_name)+"','"+ str(Domain_name)+"', '"+ str(Domain_Start)+"', '"+ str(Domain_Finish)+"', '"+ str(NumberOfSpectra)+"', '"+ str(Colours)+"', '"+ str(Percent_Covered)+"', '"+ str(Exclusive_spectrum_count)+"', '"+ str_peptides_found+"');"
        #print querry
        mycursor.execute(querry)

    mydb.commit()
    mycursor.close()
    mydb.close()

def record_Pandas_toMYSQL(pandasDF=pd.DataFrame,database="test",table="test"):
    
    mydb = mysql.connector.connect(
        host="localhost",
        user="root",
        passwd="Weafrae1")
    mycursor = mydb.cursor(buffered=True)
    try:
        # create a new database
        print("creating Database")
        mycursor.execute("CREATE DATABASE `"+database+"`;")
    except:
        print("database exists")

    mycursor.execute("use `" + database + "`;")
    querry = "CREATE TABLE `"+table+"` (id INT(11) NOT NULL AUTO_INCREMENT, "

    for column in pandasDF.columns.tolist():
        column2=column

        if (isinstance(pandasDF[column][0], float)):
            querry = querry+"`"+column2+"`"+" DOUBLE(10,4), "
        elif (isinstance(pandasDF[column][0], int)):
            querry = querry + "`"+column2+"`" + " INT(11), "

        elif (isinstance(pandasDF[column][0], str)):
            '''Sometimes this is alittle bit too long, therefore perhaps it would be good to check if the name is peptide to set the length to larger.'''

            if (column2 == 'peptides_found'):
                querry = querry + "`" + column2 + "`" + " TEXT(10000), "
            else:
                querry = querry + "`" + column2 + "`" + " TEXT(10000), "

        elif (isinstance(pandasDF[column][0], bool)):
            querry = querry + "`"+column2+"`" + " TEXT(1000), "
        elif (isinstance(pandasDF[column][0], tuple)):
            querry = querry + "`"+column2+"`" + " TEXT(10000), "
    querry=querry+" PRIMARY KEY (id))"
    try:
        mycursor.execute(querry)
    except:
        pass
    import math
    count_rec=0;
    for i_rows_recording in range(0,pandasDF.__len__()):
        row3 = pandasDF.iloc[count_rec, :]
        count_rec=count_rec+1
        querry_first_part = "INSERT INTO `"+table+"` (`id"
        querry_second_part="`) VALUES (NULL"
        for column in pandasDF.columns.tolist():
            value =row3[column]
            try:
                if (math.isnan(value)):
                    value='NULL'
            except:
                pass


            querry_first_part=querry_first_part+"` ,`"+column
            querry_second_part=querry_second_part+", '"+ str(value)+"'"
        #print querry_first_part
        querry = querry_first_part+querry_second_part+");"
        querry=querry.replace("'NULL'", "NULL")
        mycursor.execute(querry)

        mydb.commit()
    mycursor.close()
    mydb.close()

def ID_to_Gene_Name(ID):


    url = 'https://www.uniprot.org/uploadlists/'

    params = {
        'from': 'ACC+ID',
        'to': 'GENENAME',
        'format': 'tab',
        'query': ID
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = "matiss.ozols@manchester.ac.uk"  # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read(200000)
    gene_name = ((page.split("\n"))[1]).split("\t")[1]

    return gene_name

def record_Results_MS(susceptibilities_of_domains_for_UVA, MAX, MIN, Full_susceptibility=None, process=None, Result_Basket=None,experiment_name='',mapping=None):
    # this part records information that will be displayed on the html file in a Text file.

    file_object = open('/home/mbchpmo2/Write/'+experiment_name+'MS_FASTA.seq', 'w+') #database neme
    sequence = Result_Basket.sequence.iloc[0]
    file_object.write(str(Full_susceptibility['NumberOfSpectra'][0])) #'this probably should be coverage'
    file_object.write("\n")
    file_object.write(sequence)
    file_object.write("\n")
    file_object.write(str(MAX))
    file_object.write("\n")
    file_object.write(str(MIN))
    file_object.close()
    cols = ['Domain Name', 'Domain Start', 'Domain Finish', 'NumberOfSpectra','Colours']
    susceptibilities_of_domains_for_UVA=susceptibilities_of_domains_for_UVA[cols]
    susceptibilities_of_domains_for_UVA.to_csv("/home/mbchpmo2/Write/"+experiment_name+"MS_Output.csv")

    cols = ['Domain Name', 'Domain Start', 'Domain Finish', 'NumberOfSpectra']
    Downoladable = susceptibilities_of_domains_for_UVA[cols]
    Downoladable.to_csv("/home/mbchpmo2/Write/"+experiment_name+"MS_Domains.csv",index=False)
    os.chmod("/home/mbchpmo2/Write/"+experiment_name+"MS_FASTA.seq", 0o777)
    os.chmod("/home/mbchpmo2/Write/"+experiment_name+"MS_Output.csv", 0o777)
    os.chmod("/home/mbchpmo2/Write/"+experiment_name+"MS_Domains.csv", 0o777)

def record_Results(susceptibilities_of_domains_for_UVA, MAX, MIN, Full_susceptibility=None, process=None, Result_Basket=None):
    # this part records information that will be displayed on the html file.

    file_object = open('/home/mbchpmo2/Write/FASTA.seq', 'w+')
    sequence = Result_Basket.sequence.iloc[0]
    file_object.write(str(Full_susceptibility[process]))
    file_object.write("\n")
    file_object.write(sequence)
    file_object.write("\n")
    file_object.write(str(MAX))
    file_object.write("\n")
    file_object.write(str(MIN))
    file_object.close()
    cols = ['Domain Name', 'Domain Start', 'Domain Finish', 'Susceptibility','Colours']
    susceptibilities_of_domains_for_UVA=susceptibilities_of_domains_for_UVA[cols]
    susceptibilities_of_domains_for_UVA.to_csv("/home/mbchpmo2/Write/Output.csv")

    cols = ['Domain Name', 'Domain Start', 'Domain Finish', 'Susceptibility']
    Downoladable = susceptibilities_of_domains_for_UVA[cols]
    Downoladable.to_csv("/home/mbchpmo2/Write/Domains.csv",index=False)
    os.chmod("/home/mbchpmo2/Write/FASTA.seq", 0o777)
    os.chmod("/home/mbchpmo2/Write/Output.csv", 0o777)
    os.chmod("/home/mbchpmo2/Write/Domains.csv", 0o777)

def grep_regex(x,y):
    import re
    count = 0
    array = list()
    pattern = re.compile(y)
    for line in x:
        line = str(line)
        if pattern.match(line):
            array.append(count)
        count = count + 1

    return array

def reverse_colourmap(cmap, name='my_cmap_r'):
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1 - t[0], t[2], t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k, reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL)
    return my_cmap_r

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def num(s):
    try:
        return float(s)
    except ValueError:
        return int(s)

def grep(x, y):
    # print 'grep'
    count = 0
    array = list()

    for line in x:

        if y == line:
            array.append(count)

        count = count + 1
    return array

def Fasta_Analysis_Arbitarely_Domains(sequence=None, step_size=None):
    step_size=float(step_size)
    length_of_sequence = float(len(sequence))
    number_of_for_loops = length_of_sequence / step_size


    k = (number_of_for_loops).is_integer()

    if k == False:
        number_of_for_loops = int(number_of_for_loops) + 1
    elif k == True:
        number_of_for_loops = int(number_of_for_loops)

    df_with_doamin_info = pd.DataFrame()
    arbitary_domain_start = 0
    arbitary_domain_end = int(step_size)
    domain_names=[]
    domain_star=[]
    domain_finish=[]

    for process_number in range(0, number_of_for_loops):

        if process_number == 0:
            arbitary_domain_start = arbitary_domain_start
            arbitary_domain_end = arbitary_domain_end
        else:
            arbitary_domain_start = arbitary_domain_start + step_size
            arbitary_domain_end = arbitary_domain_end + step_size

        #domain_string = sequence[(arbitary_domain_start + 1):(arbitary_domain_end + 1)]
        Domain_start = arbitary_domain_start
        Domain_finish = arbitary_domain_end
        #Size_of_domain = len(domain_string)
        Domain_name = 'Domain_p'+str(Domain_start)+'_to_'+str(Domain_finish)

        if Domain_finish > length_of_sequence:
            Domain_finish = length_of_sequence+1
        if Domain_start < 1:
            Domain_start = 1

        df_with_doamin_info = df_with_doamin_info.append(
                {'Domain': Domain_name, 'Domain Start': Domain_start, 'Domain Finish': Domain_finish}, ignore_index=True)

        # this is for glycation:
    df_with_doamin_info['Domain Type'] = str(step_size)+' AA STEP'
    return df_with_doamin_info

def process_php_domain_inputs(domains=None,sequence=None):
    import sys
    import json
    df_with_doamin_info = pd.DataFrame()


    try:
        length_of_sequence = float(len(sequence))

        domains_json = json.loads(domains)
        number_loops=len(domains_json)

        for i in range(0, number_loops):

            data_array=(str(domains_json[i])).split(",")
            #print data_array
            Domain_name=data_array[2]
            Domain_name = Domain_name.replace("\n", "")
            Domain_name = Domain_name.replace("\r", "")
            #print Domain_name
            Domain_start=data_array[0]
            Domain_finish = data_array[1]

            if float(Domain_finish) > length_of_sequence:
                Domain_finish = length_of_sequence
                print ("the information entered contains domains that exceed Fasta sequence length entered")
                print(sys.exc_info())
            if float(Domain_start) < 1:
                print("the domains entered are negarive, please specify range of 0-100")
                Domain_start = 1
            df_with_doamin_info = df_with_doamin_info.append(
                {'Domain': Domain_name, 'Domain Start': Domain_start, 'Domain Finish': Domain_finish},
                ignore_index=True)

        #here have to loop through each of the json elements and create a gene entries basket.
    except:

        print ("")

    #print df_with_doamin_info
    return df_with_doamin_info



