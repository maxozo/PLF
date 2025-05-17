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
from statsmodels.sandbox.stats.multicomp import multipletests
import string

#https://www.marsja.se/three-ways-to-carry-out-2-way-anova-with-python/
#http://jpktd.blogspot.com/2013/03/multiple-comparison-and-tukey-hsd-or_25.html

def data_set_concerter(df):
    
    names = list(string.ascii_lowercase)
    key = list(df.keys())
    df_returned = pd.DataFrame()
    df_returned['Individual']=""
    df_returned['Treatment']=""
    df_returned['Domain_Name']=""
    df_returned['Yield']=""
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
    return df_returned


def calculate_P_values(MS_Resid, DF_Resid, df):
    condition_names = list(df.keys())
    num_conditions = len(condition_names)
    num_domains = len(df[condition_names[0]])   # added this for number of domains (both conditions should be same), so can use in loop (below)
    results = []                                # more efficient to start with a list, and just convert to pd at the end
    
    for i in range(num_domains):        # loop through each domain
        # Extract domain name from the first condition
        domain_name = df[condition_names[0]].iloc[i]["Domain_Name"]     # get current domain name
        row_data = {"Domain_Name": domain_name}                         # row_data is a dict - put in the current domain name

        domain_values = {       # a dictionary comprehension that gets all sample values for this domain, for each condition
            cond: df[cond].iloc[i].drop("Domain_Name").to_numpy()       # get df for a condition (df[cond]), values for a domain (.iloc[i]), without domain_name, all converted to a numpy array (.to_numpy())
            for cond in condition_names     # for each condition
        }

        for i1 in range(num_conditions):                  # loop through each pair of conditions
            for i2 in range(i1 + 1, num_conditions):
                group1 = domain_values[condition_names[i1]]         # get values from the domain_values dict ^, for each condition being compared in this pair
                group2 = domain_values[condition_names[i2]]

                mean_diff = np.mean(group1) - np.mean(group2)       # difference between mean of each condition's data
                SE = math.sqrt(MS_Resid * (1.0 / len(group1) + 1.0 / len(group2)))    # standard error calculated using MS_resid (mean square within groups, input to function) and sample size)
                t_value = mean_diff / SE                                              # t value calc
                pval = stats.t.sf(abs(t_value), DF_Resid) * 2                         # using stats package to calculate two-tailed p-value from the t-val and residual degrees of freedom (input to function)

                colname = f'p: {condition_names[i2]} vs {condition_names[i1]}'     # dynamic col name for p-vals condition1 vs condition2
                row_data[colname] = round(pval, 4)                                 # add to the row_data dict - colname as above, and value being p-val rounded to 4dp

        results.append(row_data)        # add the data for this domain to the empty results list made at start
    return pd.DataFrame(results)        # once all domains done (in all condition combos), convert to pd df

def Two_Way_mixed_Anova(df,paired=True):

    df2 = data_set_concerter(df)

    N = len(df2)
    df1=len(df2["Domain_Name"].unique())-1
    df_2=len(df2['Treatment'].unique())-1
    df_axb = df1*df_2

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

    ms_a = ssq_a / float(df1)

    ms_b = ssq_b / float(df_2)
    ms_axb = ssq_axb/df_axb
    ms_w = ssq_w/df_w

    try:
        f_a = ms_a/ms_w
        f_b = ms_b/ms_w
        f_axb = ms_axb/ms_w
    except:
        f_a = ms_a/(ms_w+0.0000000000000000001)
        f_b = ms_b/(ms_w+0.0000000000000000001)
        f_axb = ms_axb/(ms_w+0.0000000000000000001)

    p_a = stats.f.sf(f_a, df1, df_w)
    p_b = stats.f.sf(f_b, df_2, df_w)

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


    '''SS values'''
    Row_Factor_x_Time = ssq_axb
    Row_Factor = ssq_a
    Time=ssq_b
    MS_Subj = SS_Subj/DF_subjects
    MS_Resid = SS_Resid/DF_Resid


    if not paired:
        MS_Resid = ssq_w / df_w
        DF_Resid= df_w

    '''DF values'''
    results = {'sum_sq': [ssq_a, ssq_b, ssq_axb, ssq_w,SS_Subj,SS_Resid],
               'DF': [df1, df_2, df_axb, df_w,DF_subjects,DF_Resid],
                'MS': [ms_a, ms_b, ms_axb, ms_w,MS_Subj,MS_Resid],
               'F': [f_a, f_b, f_axb, 'NaN'],
               'PR(>F)': [p_a, p_b, p_axb, 'NaN']}
   
    p_values = calculate_P_values(MS_Resid,DF_Resid,df)
    p_values_adjusted = p_values.copy()
   
    
 
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

    return p_values_adjusted.fillna(1)



