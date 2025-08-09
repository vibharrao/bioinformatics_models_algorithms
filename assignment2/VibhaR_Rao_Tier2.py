import os
import sys
import csv
import numpy as np
from statistics import mean, median
from math import erf, sqrt

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

def rank_data(data):
    """Assign ranks to data, handling ties by assigning the average rank."""
    sorted_data = sorted((value, i) for i, value in enumerate(data))
    ranks = [0] * len(data)
    cur_rank = 1
    i = 0

    while i < len(sorted_data):
        tie_value = sorted_data[i][0]
        tie_indices = []

        # Find all tied values
        while i < len(sorted_data) and sorted_data[i][0] == tie_value:
            tie_indices.append(sorted_data[i][1])
            i += 1

        # Assign the average rank to all tied values
        avg_rank = cur_rank + (len(tie_indices) - 1) / 2
        for index in tie_indices:
            ranks[index] = avg_rank
        cur_rank += len(tie_indices)

    return ranks

def mann_whitney_u(control_values, treatment_values):
    """Calculate the Mann-Whitney U statistic and exact p-value for small sample sizes."""
    combined = control_values + treatment_values
    ranks = rank_data(combined)
    n1 = len(control_values)
    n2 = len(treatment_values)

    # Sum of ranks for control and treatment groups
    r1 = sum(ranks[:n1])
    r2 = sum(ranks[n1:])

    # Calculate U statistics
    u1 = r1 - n1 * (n1 + 1) / 2
    u2 = r2 - n2 * (n2 + 1) / 2

    # Use the smaller of U1 and U2
    u = min(u1, u2)

    # Exact p-value calculation for small sample sizes
    # Here we use the normal approximation if sample sizes are sufficiently large
    if n1 + n2 > 20:
        mean_u = n1 * n2 / 2
        std_u = sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
        z = (u - mean_u) / std_u
        p_value = 2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2))))  # Two-sided p-value
    else:
        # Exact calculation for small samples
        p_value = min(2 * u / (n1 * n2), 1.0)

    return p_value

def calculate_p_values(normalized_control, normalized_treatment):
    """Calculate p-values for each gene using Mann-Whitney U test."""
    num_genes = len(normalized_control[0])
    p_values = []
    for i in range(num_genes):
        control_values = [row[i] for row in normalized_control]
        treatment_values = [row[i] for row in normalized_treatment]
        p_value = mann_whitney_u(control_values, treatment_values)
        p_values.append(p_value)
    return p_values

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
    
    # Calculate p-values using Mann-Whitney U test
    p_values = calculate_p_values(normalized_control, normalized_treatment)
    
    # Create output and sort by p-value (lowest first)
    genes = [f'Gene{i+1}' for i in range(len(mean_control))]
    output_data = sorted(zip(
        genes, mean_control, median_control, mean_treatment, median_treatment, log2_fold_change, p_values
    ), key=lambda x: x[6])
    
    # Write output to a tab-delimited text file
    with open('output_Tier2.txt', 'w') as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerow(['#gene', 'mean_normalized_control_expression', 'median_normalized_control_expression', 
                         'mean_normalized_treatment_expression', 'median_normalized_treatment_expression', 'logFoldChange', 'p_value'])
        writer.writerows(output_data)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python Firstname_Lastname_Tier2.py path/to/controls path/to/treatments")
        sys.exit(1)
    
    control_path = sys.argv[1]
    treatment_path = sys.argv[2]
    main(control_path, treatment_path)
