'''The only thing here now left to do is to figure out how to adjust all the measures for the
to do this have to implement methods from the
SA Glantz and BK Slinker, Primer of Applied Regression and Analysis of Variance, McGraw-Hill, second edition, 2000.
https://accessbiomedicalscience.mhmedical.com/book.aspx?bookid=2117

http://www.real-statistics.com/anova-repeated-measures/two-within-subjects-factors/


'''

'''https://www.vlebooks.com/vleweb/Product/Index/840280?page=0
page 429, 522

https://www3.nd.edu/~rwilliam/stats1/x61.pdf


https://www.graphpad.com/guides/prism/7/statistics/index.htm?stat_anova_table_in_two_ways_rm_ano.htm
https://statistics.laerd.com/statistical-guides/repeated-measures-anova-statistical-guide-2.php

'''
import pandas as pd
from scipy import stats
import numpy as np
import math

#https://www.marsja.se/three-ways-to-carry-out-2-way-anova-with-python/
#http://jpktd.blogspot.com/2013/03/multiple-comparison-and-tukey-hsd-or_25.html

def data_set_concerter(df):
    import string
    
    names = list(string.ascii_lowercase)
    print("data_set_concerter")
    df_returned = pd.DataFrame(columns=['Individual','Treatment',"Domain_Name",'Yield'])
    print("data_set_concerter2")
    key = list(df.keys())
    print("data_set_concerter3")
    for i in range(0,df.__len__()):
        Working_Experimental_Condition = df[key[i]]
        Treatment = key[i]
        for ind in range(0,Working_Experimental_Condition.__len__()):
            Data = Working_Experimental_Condition.iloc[ind]
            Domain=Data['Domain_Name']
            Data = Data.drop(["Domain_Name"])
            for Data_line_idx in range(0,Data.__len__()):
                Yield = Data[Data_line_idx]
                Individual=names[Data_line_idx]
                df_returned=df_returned.append({'Individual':Individual,'Treatment':Treatment,"Domain_Name":Domain,'Yield':Yield},ignore_index=True)
    print("data_set_concerter_finish")
    return df_returned


def calculate_P_values(MS_Resid,DF_Resid, df):

    '''This has to be adjusted according to the 3 or more values to be compared with'''

    '''
    Now we loop through each of the domain pairs and calculate the p -values using ms_w
    SE = sqrt(ms_w*(1/Na)+(1/Nb))
    t = mean diff / SE
    then we use the degrees of freedom residual to calculate the p value
    Multiple comparisons tab here now

    '''
    list_p_values = pd.DataFrame(columns=["Domain_Name"]);

    for i in range (0,df[list(df.keys())[0]].__len__()): #for each domain
        Domain = df[list(df.keys())[0]].iloc[i]["Domain_Name"]
        list_p_values = list_p_values.append({'Domain_Name': Domain}, ignore_index=True)
        for i1 in range(0, list(df.keys()).__len__()): #for each experimental condition
            for i2 in range(i1+1,list(df.keys()).__len__()):
                if(i2<list(df.keys()).__len__()): #for each next experimental condition

                    colname = 'p: ' + list(df.keys())[i2]+' vs '+list(df.keys())[i1]
                    if(i==0):#first time going through; make a new columns for each experimental condition
                        list_p_values[colname] = ''
                    Df_Data = {list(df.keys())[i1]: df[list(df.keys())[i1]].iloc[i].drop("Domain_Name").tolist(),
                                            list(df.keys())[i2]:df[list(df.keys())[i2]].iloc[i].drop("Domain_Name").tolist()}
                    

                    # D1 = Df_Data[list(df.keys())[i1]]
                    # D2 = Df_Data[list(df.keys())[i2]]
                    M1=np.mean(Df_Data[list(df.keys())[i1]])
                    M2=np.mean(Df_Data[list(df.keys())[i2]])
                    M_diff=M1-M2
                    SE = math.sqrt(MS_Resid*(1.0/Df_Data[list(df.keys())[i1]].__len__()+1.0/Df_Data[list(df.keys())[i2]].__len__()))
                    t_value = M_diff/SE
                    Degrees_of_Freedom = DF_Resid
                    pval = stats.t.sf(np.abs(t_value), Degrees_of_Freedom)*2

                    list_p_values[colname][i] = round(pval,4)
                else:
                    continue
                    #list_p_values=list_p_values.append({'Domain':Domain,colname:pval},ignore_index=True)
                    #now test if this puts everything in the right places.
    return list_p_values

def Two_Way_mixed_Anova(df,paired=True):
    print("anova here 0")
    
    print("anova here 01")
    df2 = data_set_concerter(df)
    print("anova here 02")
    N = len(df2)
    df1=len(df2["Domain_Name"].unique())-1
    df_2=len(df2['Treatment'].unique())-1
    df_axb = df1*df_2
    print("anova here 1")
    df_w = N-(len(df2['Treatment'].unique())*(len(df2["Domain_Name"].unique())))
    grand_mean = df2['Yield'].mean() #1
    ssq_t = sum((df2.Yield - grand_mean)**2) #2

    ssq_a = sum([(df2[df2.Domain_Name ==l].Yield.mean()-grand_mean)**2 for l in df2.Domain_Name])
    ssq_b = sum([(df2[df2.Treatment ==l].Yield.mean()-grand_mean)**2 for l in df2.Treatment])
    ssq_w=0

    for group in  df2.Treatment.unique():
        vc = df2[df2.Treatment == group]
        vc_dose_means = [vc[vc.Domain_Name == d].Yield.mean() for d in vc.Domain_Name]
        ssq_w = ssq_w+sum(
            (vc.Yield - vc_dose_means) ** 2)
    ssq_axb = ssq_t-ssq_a-ssq_b-ssq_w

    print("anova here 2")
    ms_a = ssq_a / float(df1)

    ms_b = ssq_b / float(df_2)
    ms_axb = ssq_axb/df_axb
    ms_w = ssq_w/df_w

    f_a = ms_a/ms_w
    f_b = ms_b/ms_w
    f_axb = ms_axb/ms_w

    p_a = stats.f.sf(f_a, df1, df_w)
    p_b = stats.f.sf(f_b, df_2, df_w)
    print("anova here 3")
    p_axb = stats.f.sf(f_axb, df_axb, df_w)

    DF_subjects = len(df2.Individual.unique())*len(df2.Domain_Name.unique())-len(df2.Domain_Name.unique())
    DF_Resid =df_w-DF_subjects
    pre_SS_Subj = 0

    for l in df2.Domain_Name.unique():
        df1 = df2[df2.Domain_Name == l]
        AVG = df1.Yield.mean()
        for l2 in df1.Individual.unique():
            AVG2 = df1[df1.Individual == l2].Yield.mean()
            value= (AVG2-AVG)**2
            pre_SS_Subj=pre_SS_Subj+value

    Subject = SS_Subj =len(df1.Treatment.unique())*pre_SS_Subj
    Residual = SS_Resid = ssq_w-SS_Subj
    print("anova here 4")

    '''SS values'''
    Row_Factor_x_Time = ssq_axb
    Row_Factor = ssq_a
    Time=ssq_b
    MS_Subj = SS_Subj/DF_subjects
    MS_Resid = SS_Resid/DF_Resid


    if not paired:
        print("not paired")
        MS_Resid = ssq_w / df_w
        DF_Resid= df_w

    '''DF values'''
    # DF_Row_Factor_x_Time
    # DF_Row_Factor
    # DF_Time
    # DF_MS_Subj
    # DF_MS_Resid
    #
    # '''MS values'''
    # MS_Row_Factor_x_Time=ms_axb
    # MS_Row_Factor=ms_a
    # MS_Time=ms_b
    # MS_MS_Subj=MS_Subj
    # MS_MS_Resid=

    results = {'sum_sq': [ssq_a, ssq_b, ssq_axb, ssq_w,SS_Subj,SS_Resid],
               'DF': [df1, df_2, df_axb, df_w,DF_subjects,DF_Resid],
                'MS': [ms_a, ms_b, ms_axb, ms_w,MS_Subj,MS_Resid],
               'F': [f_a, f_b, f_axb, 'NaN'],
               'PR(>F)': [p_a, p_b, p_axb, 'NaN']}
    print("anova here 6")
    p_values = calculate_P_values(MS_Resid,DF_Resid,df)
    p_values_adjusted = p_values.copy()
    from statsmodels.sandbox.stats.multicomp import multipletests
    # Now feed the p-values in a bonferoni correction
    if(p_values.columns.__len__()>2):
        for i in range(0,p_values.__len__()):
            p_adjusted = multipletests(p_values_adjusted.iloc[i,1:].values.tolist(), alpha=0.05, method='bonferroni')
            p_values_adjusted.iloc[i, 1:] = p_adjusted[1]
    else:
        p_adjusted = multipletests(p_values_adjusted.iloc[:,1], alpha=0.05, method='bonferroni')
        d = p_adjusted[1]
        a = np.array(d, dtype=float)
        p_values_adjusted.iloc[:, 1]=a
    return p_values_adjusted



