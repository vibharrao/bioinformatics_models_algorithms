import sys
import numpy as np

def parse_bed_file(bed_file):
    """Parse the BED file and return a list of tuples (chrom, start, end)."""
    ranges = []  # List to store parsed ranges (chromosome, start, end).
    with open(bed_file) as f:
        for line in f:
            # Split each line by spaces or tabs, and extract the first three fields (chrom, start, end).
            chrom, start, end = line.strip().split()[:3]
            # Append the chromosome and the start/end positions as a tuple to the list of ranges.
            ranges.append((chrom, int(start), int(end)))
    return ranges  # Return the list of parsed ranges.

def parse_fai_file(fai_file):
    """Parse the .fai file and return a dictionary of chromosome lengths."""
    chrom_lengths = {}  # Dictionary to store chromosome lengths.
    with open(fai_file) as f:
        for line in f:
            # Split each line to get chromosome name and its length.
            chrom, length = line.strip().split()[:2]
            chrom_lengths[chrom] = int(length)  # Store chromosome and its length in the dictionary.
    return chrom_lengths  # Return the dictionary of chromosome lengths.

def merge_ranges(ranges):
    """Merge overlapping or adjacent ranges."""
    merged = []  # List to store merged ranges.
    # Sort ranges by chromosome and start position to ensure they are processed in order.
    for chrom, start, end in sorted(ranges):
        # If merged list is empty or there is no overlap with the last merged range, append the new range.
        if not merged or merged[-1][0] != chrom or merged[-1][2] < start:
            merged.append([chrom, start, end])
        else:
            # If there is overlap or adjacency, merge by updating the end position of the last range.
            merged[-1][2] = max(merged[-1][2], end)
    return merged  # Return the list of merged ranges.

def calculate_overlap(merged_a, merged_b):
    """Calculate the number of overlapping bases between two merged sets of ranges."""
    overlap = 0  # Variable to store the total overlap.
    i, j = 0, 0  # Pointers to iterate through both merged range lists.
    
    # While there are ranges to compare in both sets:
    while i < len(merged_a) and j < len(merged_b):
        a_chrom, a_start, a_end = merged_a[i]  # Current range from Set A.
        b_chrom, b_start, b_end = merged_b[j]  # Current range from Set B.
        
        # If chromosomes are different, move the pointer in the smaller chromosome.
        if a_chrom < b_chrom:
            i += 1
        elif a_chrom > b_chrom:
            j += 1
        else:  # If the chromosomes are the same, compare the ranges.
            if a_end < b_start:  # No overlap, move Set A pointer.
                i += 1
            elif b_end < a_start:  # No overlap, move Set B pointer.
                j += 1
            else:
                # Overlap exists; calculate the overlap size and add to total.
                overlap += min(a_end, b_end) - max(a_start, b_start)
                # Move the pointer in the set where the range ends first.
                if a_end < b_end:
                    i += 1
                else:
                    j += 1
    return overlap  # Return the total overlap between the two sets.

def positions_permuted(ranges, chrom_lengths):
    """Randomly permute the ranges by assigning non-overlapping start positions."""
    assigned_ranges = []  # List to store the permuted ranges.

    # Loop over each chromosome in the chromosome lengths dictionary.
    for chrom in chrom_lengths:
        # Filter for ranges belonging to the current chromosome.
        chrom_ranges = [r for r in ranges if r[0] == chrom]
        if not chrom_ranges:
            continue  # Skip if no ranges are present for this chromosome.

        chrom_length = chrom_lengths[chrom]  # Get the total length of the chromosome.
        total_range_length = sum(end - start for _, start, end in chrom_ranges)  # Total length of all ranges on this chromosome.
        available_space = chrom_length - total_range_length  # Calculate the available space for random gaps.

        # If there isn't enough space on the chromosome for all ranges, raise an error.
        if available_space < 0:
            raise ValueError(f"Not enough space on chromosome {chrom} to fit all ranges.")
        
        # Sort ranges by their length in descending order to handle larger ranges first.
        chrom_ranges.sort(key=lambda x: x[2] - x[1], reverse=True)
        
        # Generate random gaps between the ranges, such that the total gaps fit the available space.
        num_gaps = len(chrom_ranges) + 1  # Number of gaps is one more than the number of ranges.
        gap_probabilities = [1 / num_gaps] * num_gaps  # Equal probabilities for each gap.
        gaps = np.random.multinomial(available_space, gap_probabilities)  # Randomly assign the gaps.

        # Place each range sequentially with the generated random gaps.
        current_position = gaps[0]  # Start position after the first gap.
        for i, (_, start, end) in enumerate(chrom_ranges):
            length = end - start  # Length of the current range.
            new_start = current_position  # New start position for the permuted range.
            new_end = new_start + length  # New end position.
            assigned_ranges.append((chrom, new_start, new_end))  # Add the new range to the assigned list.
            current_position = new_end + gaps[i + 1]  # Move to the next start position with a gap.
    
    return assigned_ranges  # Return the list of permuted ranges.

def randomize_ranges_together(set_a, set_b, chrom_lengths):
    """Randomly shuffle Set A and Set B together, then assign non-overlapping positions to each set."""
    all_ranges = set_a + set_b  # Combine ranges from both Set A and Set B.
    np.random.shuffle(all_ranges)  # Shuffle the combined ranges randomly.
    
    # Split the shuffled ranges back into Set A and Set B.
    num_ranges_a = len(set_a)  # Number of ranges in Set A.
    shuffled_set_a = all_ranges[:num_ranges_a]  # First portion goes to Set A.
    shuffled_set_b = all_ranges[num_ranges_a:]  # Second portion goes to Set B.
    
    # Assign non-overlapping positions to both sets.
    randomized_set_a = positions_permuted(shuffled_set_a, chrom_lengths)
    randomized_set_b = positions_permuted(shuffled_set_b, chrom_lengths)
    
    return randomized_set_a, randomized_set_b  # Return the permuted sets A and B.

def permutation_test(set_a, set_b, chrom_lengths, num_permutations=10000):
    """Perform a two-tailed permutation test to calculate the significance of the observed overlap."""
    merged_a = merge_ranges(set_a)  # Merge overlapping/adjacent ranges in Set A.
    merged_b = merge_ranges(set_b)  # Merge overlapping/adjacent ranges in Set B.
    observed_overlap = calculate_overlap(merged_a, merged_b)  # Calculate observed overlap between Set A and Set B.
    
    permuted_overlaps = []  # List to store overlaps from each permutation.
    
    # Run the permutation test for a given number of iterations.
    for _ in range(num_permutations):
        permuted_a, permuted_b = randomize_ranges_together(set_a, set_b, chrom_lengths)  # Randomize positions.
        merged_permuted_a = merge_ranges(permuted_a)  # Merge permuted Set A ranges.
        merged_permuted_b = merge_ranges(permuted_b)  # Merge permuted Set B ranges.
        permuted_overlap = calculate_overlap(merged_permuted_a, merged_permuted_b)  # Calculate overlap for permuted sets.
        permuted_overlaps.append(permuted_overlap)  # Store the permuted overlap.
        
    permuted_overlaps = np.array(permuted_overlaps)  # Convert the list to a numpy array.
    
    # Calculate the p-value: proportion of permutations with overlap greater than or equal to the observed overlap.
    p_value = np.sum(permuted_overlaps >= observed_overlap) / num_permutations
    
    return observed_overlap, p_value  # Return the observed overlap and the calculated p-value.

def main():
    # Parse command line arguments
    set_a_file = sys.argv[1]  # BED file for Set A.
    set_b_file = sys.argv[2]  # BED file for Set B.
    genome_fai_file = sys.argv[3]  # FAI file with chromosome lengths.
    num_permutations = int(sys.argv[4]) if len(sys.argv) > 4 else 10000  # Number of permutations (default 10000).
    
    # Parse the input files
    set_a = parse_bed_file(set_a_file)  # Parse Set A BED file.
    set_b = parse_bed_file(set_b_file)  # Parse Set B BED file.
    chrom_lengths = parse_fai_file(genome_fai_file)  # Parse FAI file for chromosome lengths.
    
    # Perform the permutation test
    observed_overlap, p_value = permutation_test(set_a, set_b, chrom_lengths, num_permutations)
    
    # Output the result
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

if __name__ == "__main__":
    main()  # Run the main function when executed as a script.
