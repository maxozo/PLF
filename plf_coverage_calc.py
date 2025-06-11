import argparse
import pandas as pd
import sys
import os

# Set up argument parser
parser = argparse.ArgumentParser(description="Calculate segment coverage from PLF results.")

parser.add_argument('--input_file', required=True, help='Input PLF results .tsv file')
parser.add_argument('--condition1', required=True, help='First condition column name')
parser.add_argument('--condition2', required=True, help='Second condition column name')

# Parse the arguments
args = parser.parse_args()

# Access arguments like this:
plf_results_file = args.input_file
condition1 = args.condition1
condition2 = args.condition2

input_data_name = os.path.splitext(os.path.basename(plf_results_file))[0]

print(f"Reading input file: {plf_results_file}")

# Try reading the file
try:
    plf_results = pd.read_csv(plf_results_file, sep="\t")
except Exception as e:
    print(f"Error reading input file: {e}")
    sys.exit(1)

# Segment coverage calculation function
def calculate_segment_coverage(df):
    print("Calculating segment coverage...")
   
    result_rows = []
    proportion_covered_dict = {}  # stores geneAC and coverage

    # Check if the provided condition names exist
    try:
        df[[condition1, condition2]]
    except KeyError as e:
        print(f"Error: One or both of the specified conditions were not found in the data. \n Please enter the conditions as they appear in the columns.")
        # print(f"-{', '.join(df.columns)}")
        print(f"Details: {e}")
        sys.exit(1)

    grouped = df.groupby('GeneAC')

    for gene, group in grouped:
        num_segments = len(group)
        num_segments_present = (group[condition1] != 0).sum()

        if num_segments_present == 0:
            num_segments_present = (group[condition2] != 0).sum()

        proportion_covered = num_segments_present / num_segments if num_segments > 0 else None

        result_rows.append({
            'GeneAC': gene,
            'num_segments': num_segments,
            'num_segments_present': num_segments_present,
            'proportion_covered': proportion_covered
        })

        proportion_covered_dict[gene] = proportion_covered

    result_df = pd.DataFrame(result_rows)
    result_df.to_csv(f"{input_data_name}_segment_coverage_summary_table.csv", index=False)
    
    df['protein_segment_coverage'] = df['GeneAC'].map(proportion_covered_dict)
    df.to_csv(f"{input_data_name}_with_segment_coverage.csv", index=False)
    
    return result_df


# Run the function
calculate_segment_coverage(plf_results)

print(f"Segment coverage calculation finished. Results table saved to {input_data_name}_segment_coverage_summary_table.csv")
print(f"Input data with segment coverage added saved to {input_data_name}_with_segment_coverage.csv")




# example input:      
#   python plf_coverage_calc.py --input_file Aim2_EDTA_PLF_BMprots_PLFresults.tsv --condition1 Epidermis+ --condition2 Epidermis-


