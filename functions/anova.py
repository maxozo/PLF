
import pandas as pd
import string
from scipy import stats
import numpy as np
import math
from statsmodels.stats.multitest import multipletests


def data_set_converter(df_dict):
    data_list = [
        {'Individual': col, 'Treatment': treatment, 'Domain_Name': row.Domain_Name, 'Yield': getattr(row, col)}
        for treatment, condition_df in df_dict.items()
        for row in condition_df.itertuples(index=False)
        for col in condition_df.columns if col != 'Domain_Name'
    ]
    return pd.DataFrame(data_list)


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



def Two_Way_mixed_Anova_2(df, paired=True):

    long_df = data_set_converter(df)  # use the above function to make data into long format, with cols for Individual, Treatment, Domain_Name, & Yield
    # stuff that is used a few times, so more efficient to get here:
    n_rows = len(long_df)  # num rows in the long-format data
    domains = long_df["Domain_Name"].unique()     # domains
    treatments = long_df["Treatment"].unique()    # experimental conditions
    individuals = long_df["Individual"].unique()  # samples

    df_A = len(domains) - 1  # degrees of freedom for Domain_Name (factor A) = number of domains in protein - 1
    df_B = len(treatments) - 1  # degrees of freedom for Treatment (factor B) = number of conditions - 1
    df_axb = df_A * df_B                # these two ^^ multiplied - gives the degrees of freedom for interaction effect (df_axb)
    df_w = n_rows - len(treatments) * len(domains)  # Degrees of freedom for within-group (residual error) = num rows in long data - (num conditions * num domains)

    grand_mean = long_df['Yield'].mean()  # mean of all yield/abundance values (all domains, all samples, all conditions)

    ssq_t = np.sum((long_df['Yield'] - grand_mean) ** 2)  # Total Sum of Squares (overall variation from grand mean)
    ssq_a = sum([  # Sum of Squares for Domain_Name (Factor A: domain effect)
        (long_df[long_df.Domain_Name == l].Yield.mean() - grand_mean)**2    # per row in long_df: get the mean yield of each domain (across samples/conditions), minus grand mean, square it. Add up values from all rows (there will be the same value on multiple rows - yes add)
        for l in long_df.Domain_Name
    ]) 
    ssq_b = sum([  # Sum of Squares for Treatment (Factor B: treatment/condition effect) 
            (long_df[long_df.Treatment == l].Yield.mean() - grand_mean)**2   # same as above, but per treatment (condition)
            for l in long_df.Treatment        # for each treatment (condition) in df2
    ]) 
    
    ssq_w = 0   # Initialize the Within-Group Sum of Squares (SSW), which represents unexplained variability
    for group in treatments:  # for each condition
        vc = long_df[long_df['Treatment'] == group]  # get all rows in long_df for current condition
        vc_domain_means = vc.groupby('Domain_Name')['Yield'].transform('mean')  # the mean yield per domain, for the current condition
        ssq_w += np.sum((vc['Yield'] - vc_domain_means) ** 2)   # add on to ssq_w (started at 0 outside of this loop), per row the sum of: yield - mean domain yield, squared - done per condition

    ssq_axb = ssq_t - ssq_a - ssq_b - ssq_w  # Interaction Sum of Squares: how much variation is left after removing main effects and error

    # Now calculate Mean Squares (MS = SS / DF) for each factor and interaction
    ms_a = ssq_a / df_A  # MS for Domain_Name (Factor A)
    ms_b = ssq_b / df_B  # MS for Treatment (Factor B)
    ms_axb = ssq_axb / df_axb  # MS for Interaction (A x B)
    ms_w = ssq_w / df_w  # MS for within-group error (residual)

    try:
        # Calculate F-statistics: MS of each factor divided by MS of the error
        f_a = ms_a / ms_w
        f_b = ms_b / ms_w
        f_axb = ms_axb / ms_w
    except ZeroDivisionError:
        # Handle edge case where ms_w might be zero: add on a small number, then calculate mean squares
        ms_w += 1e-19
        f_a = ms_a / ms_w
        f_b = ms_b / ms_w
        f_axb = ms_axb / ms_w

    p_a = stats.f.sf(f_a, df_A, df_w)  # use stats package to calculate p-vals 
    p_b = stats.f.sf(f_b, df_B, df_w)
    p_axb = stats.f.sf(f_axb, df_axb, df_w)

    DF_subjects = len(individuals) * len(domains) - len(domains)  # degrees of freedom for samples, within each domain  (num unique samples * num of unique domains - num unique domains)
    DF_Resid = df_w - DF_subjects  # Degrees of freedom for residuals = Degrees of freedom for within-group (residual error) - degrees of freedom for samples

    # useful to store:
    domain_means = long_df.groupby('Domain_Name')['Yield'].mean()     # mean yield per domain (across conditions/samples)
    subject_means = long_df.groupby(['Domain_Name', 'Individual'])['Yield'].mean()      # mean yield per domain per sample (across all conditions)

    pre_SS_Subj = np.sum((subject_means - subject_means.index.get_level_values(0).map(domain_means)) ** 2)  # pre- sum of squares for samples. For each sample, do subject_mean - domain mean, square it, then all sample values added together

    SS_Subj = len(treatments) * pre_SS_Subj  # Final sum of squares for samples: multiply pre-value (^) by number of conditions
    SS_Resid = ssq_w - SS_Subj  # Final residual sum of squares:  take the sample SS (^) away from Within-Group Sum of Squares (total unexplained variation, accounting for samples)

    '''SS values'''
    MS_Subj = SS_Subj / DF_subjects  # mean square for samples = sum of squares for samples / degrees of freedom for samples
    MS_Resid = SS_Resid / DF_Resid  # mean square for residuals = residual sum of squares / Degrees of freedom for residuals
    if not paired:  # If the data is not paired (repeated measures), revert to the basic within-group MS and DF
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

    p_values = calculate_P_values(MS_Resid, DF_Resid, df)   # use p-vals function, put results in 'p_values'
    p_values_adjusted = p_values.copy()                     # make a copy of 'p_values'

    # Now make adjusted p-values with Bonferroni correction
    if p_values.shape[1] > 2:  # if there is more than 1 comparison of conditions (i.e. if data had >2 conditions) (p-val results table has 2 cols if there is just 1 comparison)
        for i in range(len(p_values)):  # for each row in p-val table (rows are domains - so for each domain):
            adj = multipletests(
                p_values_adjusted.iloc[i, 1:].values.tolist(),  # for each column after the first (so each comparison) (in a list)
                alpha=0.05,
                method='bonferroni'
            )[1]
            p_values_adjusted.iloc[i, 1:] = adj  # replace og p-val with adj-p value
    else:           
        adj = multipletests(                     # (same as above, but for when there is just 1 comparison col)
            p_values_adjusted.iloc[:, 1],
            alpha=0.05,
            method='bonferroni'
        )[1]
        p_values_adjusted.iloc[:, 1] = np.array(adj, dtype=float)

    return p_values_adjusted.fillna(1)  # return the final adjusted p-vals, with any missing p-vals put as 1 (ns)



