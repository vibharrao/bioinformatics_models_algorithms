import os
import sys
import csv
import numpy as np
from statistics import mean, median

def normalize_counts(data):
    """Normalize the gene expression counts within each replicate."""
    normalized_data = []
    for row in data:
        total = sum(row)
        if total == 0:
            # Avoid division by zero by normalizing to zeros if the total is zero
            normalized_row = [0 for value in row]
        else:
            normalized_row = [value / total for value in row]
        normalized_data.append(normalized_row)
    return normalized_data

def compute_stats(normalized_data):
    """Compute mean and median normalized expression for each gene."""
    num_genes = len(normalized_data[0])
    means = [mean([row[i] for row in normalized_data]) for i in range(num_genes)]
    medians = [median([row[i] for row in normalized_data]) for i in range(num_genes)]
    return means, medians

def calculate_log2_fold_change(mean_treatment, mean_control):
    """Compute log2 Fold Change for each gene, handling division by zero."""
    log2_fc = []
    for treat, ctrl in zip(mean_treatment, mean_control):
        if ctrl == 0:
            # If control is zero, set fold change to a large positive value if treatment is > 0
            log2_fc.append(float('inf') if treat > 0 else 0)
        else:
            log2_fc.append(np.log2(treat / ctrl) if treat > 0 else -float('inf'))
    return log2_fc

def mann_whitney_u_test(control, treatment):
    """Perform Mann-Whitney U Test on two sets of gene expression data (control and treatment)."""
    n1 = len(control)
    n2 = len(treatment)

    # Rank all the values combined
    all_values = sorted(control + treatment)
    ranks = {val: rank+1 for rank, val in enumerate(all_values)}

    # Calculate rank sums
    R1 = sum(ranks[val] for val in control)
    R2 = sum(ranks[val] for val in treatment)

    # Calculate U1 and U2
    U1 = R1 - (n1 * (n1 + 1)) / 2
    U2 = R2 - (n2 * (n2 + 1)) / 2

    # Return the smaller U value (because we want the one with smaller significance)
    U = min(U1, U2)

    # Calculate p-value approximation for large samples (use normal approximation of U statistic)
    mean_U = n1 * n2 / 2
    std_U = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z_score = (U - mean_U) / std_U
    p_value = 2 * min(1, 1 - abs(z_score))  # Two-tailed test approximation

    return p_value

def read_csv_files(file_paths):
    """Read multiple CSV files and concatenate the data."""
    data = []
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip header
            row_values = []
            for row in reader:
                row_values.append(float(row[1]))  # Extract the expression level (second column)
            data.append(row_values)
    return data

def main(control_path, treatment_path):
    # Load all control and treatment files
    control_files = []
    for f in os.listdir(control_path):
        if f.endswith('.csv'):
            control_files.append(os.path.join(control_path, f))

    treatment_files = []
    for f in os.listdir(treatment_path):
        if f.endswith('.csv'):
            treatment_files.append(os.path.join(treatment_path, f))
    
    # Read and concatenate all control and treatment data
    control_data = read_csv_files(control_files)
    treatment_data = read_csv_files(treatment_files)
    
    # Normalize the gene expression counts
    normalized_control = normalize_counts(control_data)
    normalized_treatment = normalize_counts(treatment_data)
    
    # Compute statistics
    mean_control, median_control = compute_stats(normalized_control)
    mean_treatment, median_treatment = compute_stats(normalized_treatment)
    
    # Compute log2 Fold Change
    log2_fold_change = calculate_log2_fold_change(mean_treatment, mean_control)
    
    # Compute p-values using Mann-Whitney U Test
    p_values = []
    for i in range(len(normalized_control[0])):
        control_gene_values = [row[i] for row in normalized_control]
        treatment_gene_values = [row[i] for row in normalized_treatment]
        p_value = mann_whitney_u_test(control_gene_values, treatment_gene_values)
        p_values.append(p_value)
    
    # Create output and sort by p-value (lowest first)
    genes = [f'Gene{i+1}' for i in range(len(mean_control))]
    output_data = sorted(zip(
        genes, mean_control, median_control, mean_treatment, median_treatment, log2_fold_change, p_values
    ), key=lambda x: x[6])
    
    # Write output to a tab-delimited text file
    with open('output_Tier3.txt', 'w') as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerow(['#gene', 'mean normalized control expression', 'median normalized control expression', 
                         'mean normalized treatment expression', 'median normalized treatment expression', 
                         'logFoldChange', 'p value'])
        writer.writerows(output_data)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python your_assignment.py path/to/controls path/to/treatments")
        sys.exit(1)
    
    control_path = sys.argv[1]
    treatment_path = sys.argv[2]
    main(control_path, treatment_path)
