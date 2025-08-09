import sys
import numpy as np

# Function to parse a BED file and return a list of genomic ranges (chromosome, start, end).
def parse_bed_file(bed_file):
    """Parse the BED file and return a list of tuples (chrom, start, end)."""
    ranges = []
    with open(bed_file) as f:
        # Read each line, extract chromosome, start, and end, and append them to the ranges list.
        for line in f:
            chrom, start, end = line.strip().split()[:3]
            ranges.append((chrom, int(start), int(end)))
    return ranges

# Function to parse a .fai file and return a dictionary with chromosome lengths.
def parse_fai_file(fai_file):
    """Parse the .fai file and return a dictionary of chromosome lengths."""
    chrom_lengths = {}
    with open(fai_file) as f:
        # Read each line, extract chromosome name and length, and add them to the dictionary.
        for line in f:
            chrom, length = line.strip().split()[:2]
            chrom_lengths[chrom] = int(length)
    return chrom_lengths

# Function to merge overlapping or adjacent ranges.
def merge_ranges(ranges):
    """Merge overlapping or adjacent ranges."""
    merged = []
    # Sort the ranges by chromosome, start, and end.
    for chrom, start, end in sorted(ranges):
        # If the merged list is empty or the current range doesn't overlap with the last range, add it as a new range.
        if not merged or merged[-1][0] != chrom or merged[-1][2] < start:
            merged.append([chrom, start, end])
        else:
            # If the current range overlaps with the last one, merge them by updating the end position.
            merged[-1][2] = max(merged[-1][2], end)
    return merged

# Function to calculate the number of overlapping bases between two sets of merged ranges.
def calculate_overlap(merged_a, merged_b):
    """Calculate the number of overlapping bases between two merged sets of ranges."""
    overlap = 0
    i, j = 0, 0
    # Iterate through both sets of ranges.
    while i < len(merged_a) and j < len(merged_b):
        a_chrom, a_start, a_end = merged_a[i]
        b_chrom, b_start, b_end = merged_b[j]
        
        # Skip ranges on different chromosomes.
        if a_chrom < b_chrom:
            i += 1
        elif a_chrom > b_chrom:
            j += 1
        else:  # Same chromosome
            # If the ranges don't overlap, move to the next range in one of the sets.
            if a_end < b_start:
                i += 1
            elif b_end < a_start:
                j += 1
            else:
                # Calculate the overlap between the two ranges and add it to the total overlap.
                overlap += min(a_end, b_end) - max(a_start, b_start)
                # Move to the next range in one of the sets.
                if a_end < b_end:
                    i += 1
                else:
                    j += 1
    return overlap

# Function to permute the positions of the ranges randomly while ensuring no overlap within each set.
def positions_permuted(ranges, chrom_lengths):
    """Randomly permute the ranges by assigning non-overlapping start positions."""
    assigned_ranges = []
    
    # Iterate through each chromosome in the genome.
    for chrom in chrom_lengths:
        # Filter the ranges that belong to the current chromosome.
        chrom_ranges = [r for r in ranges if r[0] == chrom]
        if not chrom_ranges:
            continue  # Skip if no ranges for this chromosome
        
        chrom_length = chrom_lengths[chrom]
        # Calculate the total length of all ranges on this chromosome.
        total_range_length = sum(end - start for _, start, end in chrom_ranges)
        available_space = chrom_length - total_range_length
        
        if available_space < 0:
            raise ValueError(f"Not enough space on chromosome {chrom} to fit all ranges.")
        
        # Sort the ranges by length in descending order to place longer ranges first.
        chrom_ranges.sort(key=lambda x: x[2] - x[1], reverse=True)
        
        # Divide the available space into gaps between the ranges.
        num_gaps = len(chrom_ranges) + 1
        gap_probabilities = [1 / num_gaps] * num_gaps  # Equal probabilities for all gaps
        gaps = np.random.multinomial(available_space, gap_probabilities)
        
        # Place the ranges sequentially with the random gaps.
        current_position = gaps[0]  # Start after the first gap
        for i, (_, start, end) in enumerate(chrom_ranges):
            length = end - start
            new_start = current_position
            new_end = new_start + length
            assigned_ranges.append((chrom, new_start, new_end))
            current_position = new_end + gaps[i + 1]  # Add the next gap after this range
    
    return assigned_ranges

# Function to shuffle the ranges of Set A and Set B together, then assign non-overlapping positions.
def randomize_ranges_together(set_a, set_b, chrom_lengths):
    """Randomly shuffle Set A and Set B together, then assign non-overlapping positions to each set."""
    # Concatenate Set A and Set B, then shuffle the combined set.
    all_ranges = set_a + set_b
    np.random.shuffle(all_ranges)
    
    # Split the shuffled set back into Set A and Set B.
    num_ranges_a = len(set_a)
    shuffled_set_a = all_ranges[:num_ranges_a]
    shuffled_set_b = all_ranges[num_ranges_a:]
    
    # Assign non-overlapping positions for Set A and Set B using the `positions_permuted` function.
    randomized_set_a = positions_permuted(shuffled_set_a, chrom_lengths)
    randomized_set_b = positions_permuted(shuffled_set_b, chrom_lengths)
    
    return randomized_set_a, randomized_set_b

# Function to perform a permutation test and calculate the p-value.
def permutation_test(set_a, set_b, chrom_lengths, num_permutations=10000):
    """Perform a permutation test to calculate the significance of the observed overlap."""
    # Merge the ranges in Set A and Set B.
    merged_a = merge_ranges(set_a)
    merged_b = merge_ranges(set_b)
    # Calculate the observed overlap between the two sets.
    observed_overlap = calculate_overlap(merged_a, merged_b)
    
    permuted_overlaps = []
    
    # Perform the specified number of permutations.
    for _ in range(num_permutations):
        # Randomly shuffle the ranges and assign new positions.
        permuted_a, permuted_b = randomize_ranges_together(set_a, set_b, chrom_lengths)
        merged_permuted_a = merge_ranges(permuted_a)
        merged_permuted_b = merge_ranges(permuted_b)
        # Calculate the overlap for the permuted sets.
        permuted_overlap = calculate_overlap(merged_permuted_a, merged_permuted_b)
        permuted_overlaps.append(permuted_overlap)
        
    permuted_overlaps = np.array(permuted_overlaps)
    
    # Calculate the p-value by determining how many permuted overlaps are greater than or equal to the observed overlap.
    p_value = np.sum(permuted_overlaps >= observed_overlap) / num_permutations
    
    return observed_overlap, p_value

# Main function to handle input and output.
def main():
    # Parse command line arguments.
    set_a_file = sys.argv[1]
    set_b_file = sys.argv[2]
    genome_fai_file = sys.argv[3]
    num_permutations = int(sys.argv[4]) if len(sys.argv) > 4 else 10000
    
    # Parse the BED files for Set A and Set B, and the genome .fai file for chromosome lengths.
    set_a = parse_bed_file(set_a_file)
    set_b = parse_bed_file(set_b_file)
    chrom_lengths = parse_fai_file(genome_fai_file)
    
    # Perform the permutation test.
    observed_overlap, p_value = permutation_test(set_a, set_b, chrom_lengths, num_permutations)
    
    # Output the results.
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

# Run the main function if this script is executed.
if __name__ == "__main__":
    main()
