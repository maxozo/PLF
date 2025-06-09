import pandas as pd
import string
from scipy import stats
import numpy as np
import math
from statsmodels.stats.multitest import multipletests




def data_set_converter(df):
    names = list(string.ascii_lowercase)
    key = list(df.keys())
    data_list = []

    for x in range(len(df)):
        Working_Experimental_Condition = df[key[x]]
        Treatment = key[x]

        for row in range(len(Working_Experimental_Condition)):
            Data = Working_Experimental_Condition.iloc[row]
            Domain = Data['Domain_Name']
            Data = Data.drop(["Domain_Name"])

            for line_index, Yield in enumerate(Data):
                if line_index >= len(names):
                    Individual = f"sample_{line_index}"
                else:
                    Individual = names[line_index]
                data_list.append({
                    'Individual': Individual,
                    'Treatment': Treatment,
                    'Domain_Name': Domain,
                    'Yield': Yield
                })
    df_returned = pd.DataFrame(data_list)
    return df_returned





def calculate_P_values(MS_Resid, DF_Resid, df):
    condition_names = list(df.keys())
    num_conditions = len(condition_names)
    num_domains = len(df[condition_names[0]])
    results = []
    for i in range(num_domains):
        domain_name = df[condition_names[0]].iloc[i]["Domain_Name"]
        row_data = {"Domain_Name": domain_name}
        domain_values = {
            cond: df[cond].iloc[i].drop("Domain_Name").to_numpy()
            for cond in condition_names
        }
        for i1 in range(num_conditions):
            for i2 in range(i1 + 1, num_conditions):
                group1 = domain_values[condition_names[i1]]
                group2 = domain_values[condition_names[i2]]
                mean_diff = np.mean(group1) - np.mean(group2)
                SE = math.sqrt(MS_Resid * (1.0 / len(group1) + 1.0 / len(group2)))
                t_value = mean_diff / SE
                pval = stats.t.sf(abs(t_value), DF_Resid) * 2
                colname = f'p: {condition_names[i2]} vs {condition_names[i1]}'
                row_data[colname] = round(pval, 4)
        results.append(row_data)
    return pd.DataFrame(results)




def Two_Way_mixed_Anova(df, paired=True):
    long_df = data_set_converter(df)
    n_rows = len(long_df)
    domains = long_df["Domain_Name"].unique()
    treatments = long_df["Treatment"].unique()
    individuals = long_df["Individual"].unique()
    df_A = len(domains) - 1
    df_B = len(treatments) - 1
    df_axb = df_A * df_B
    df_w = n_rows - len(treatments) * len(domains)
    grand_mean = long_df['Yield'].mean()
    ssq_t = np.sum((long_df['Yield'] - grand_mean) ** 2)
    ssq_a = sum([
        (long_df[long_df.Domain_Name == l].Yield.mean() - grand_mean)**2
        for l in long_df.Domain_Name
    ]) 
    ssq_b = sum([
            (long_df[long_df.Treatment == l].Yield.mean() - grand_mean)**2
            for l in long_df.Treatment
    ]) 
    ssq_w = 0
    for group in treatments:
        vc = long_df[long_df['Treatment'] == group]
        vc_domain_means = vc.groupby('Domain_Name')['Yield'].transform('mean')
        ssq_w += np.sum((vc['Yield'] - vc_domain_means) ** 2)
    ssq_axb = ssq_t - ssq_a - ssq_b - ssq_w
    ms_a = ssq_a / df_A
    ms_b = ssq_b / df_B
    ms_axb = ssq_axb / df_axb
    ms_w = ssq_w / df_w
    try:
        f_a = ms_a / ms_w
        f_b = ms_b / ms_w
        f_axb = ms_axb / ms_w
    except ZeroDivisionError:
        ms_w += 1e-19
        f_a = ms_a / ms_w
        f_b = ms_b / ms_w
        f_axb = ms_axb / ms_w
    p_a = stats.f.sf(f_a, df_A, df_w)
    p_b = stats.f.sf(f_b, df_B, df_w)
    p_axb = stats.f.sf(f_axb, df_axb, df_w)
    DF_subjects = len(individuals) * len(domains) - len(domains)
    DF_Resid = df_w - DF_subjects
    domain_means = long_df.groupby('Domain_Name')['Yield'].mean()
    subject_means = long_df.groupby(['Domain_Name', 'Individual'])['Yield'].mean()
    pre_SS_Subj = np.sum((subject_means - subject_means.index.get_level_values(0).map(domain_means)) ** 2)
    SS_Subj = len(treatments) * pre_SS_Subj
    SS_Resid = ssq_w - SS_Subj
    '''SS values'''
    MS_Subj = SS_Subj / DF_subjects
    MS_Resid = SS_Resid / DF_Resid
    if not paired:
        MS_Resid = ms_w
        DF_Resid = df_w
    '''DF values'''
    results = {
        'sum_sq': [ssq_a, ssq_b, ssq_axb, ssq_w, SS_Subj, SS_Resid],
        'DF': [df_A, df_B, df_axb, df_w, DF_subjects, DF_Resid],
        'MS': [ms_a, ms_b, ms_axb, ms_w, MS_Subj, MS_Resid],
        'F': [f_a, f_b, f_axb, 'NaN'],
        'PR(>F)': [p_a, p_b, p_axb, 'NaN']
    }
    p_values = calculate_P_values(MS_Resid, DF_Resid, df)
    p_values_adjusted = p_values.copy()
    if p_values.shape[1] > 2:
        for i in range(len(p_values)):
            adj = multipletests(
                p_values_adjusted.iloc[i, 1:].values.tolist(),
                alpha=0.05,
                method='bonferroni'
            )[1]
            p_values_adjusted.iloc[i, 1:] = adj
    else:           
        adj = multipletests(
            p_values_adjusted.iloc[:, 1],
            alpha=0.05,
            method='bonferroni'
        )[1]
        p_values_adjusted.iloc[:, 1] = np.array(adj, dtype=float)
    return p_values_adjusted.fillna(1) 