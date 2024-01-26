import pandas as pd
import re
import sys
import json
    
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


def Get_Domains(data=None):
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

def grep_regex(x,y):
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



