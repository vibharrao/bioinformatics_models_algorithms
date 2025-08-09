import sys
import numpy as np

def parse_bed(file_path):
    ranges = []
    with open(file_path) as f:
        for line in f:
            chrom, start, end = line.strip().split()[:3]
            ranges.append((chrom, int(start), int(end)))
    return ranges

def parse_fai(file_path):
    chrom_lengths = {}
    with open(file_path) as f:
        for line in f:
            parts = line.strip().split()
            chrom = parts[0]
            length = int(parts[1])
            chrom_lengths[chrom] = length
    return chrom_lengths

def calculate_overlap(ranges_a, ranges_b):
    overlap_bases = 0
    # Sort ranges by chromosome and start position
    ranges_a = sorted(ranges_a, key=lambda x: (x[0], x[1]))
    ranges_b = sorted(ranges_b, key=lambda x: (x[0], x[1]))

    i, j = 0, 0
    while i < len(ranges_a) and j < len(ranges_b):
        chrom_a, start_a, end_a = ranges_a[i]
        chrom_b, start_b, end_b = ranges_b[j]
        
        if chrom_a == chrom_b:
            # Check for overlap
            if start_a < end_b and end_a > start_b:
                # Calculate overlap length
                overlap_start = max(start_a, start_b)
                overlap_end = min(end_a, end_b)
                overlap_bases += overlap_end - overlap_start
            
            # Move to the next range in Set A or Set B
            if end_a < end_b:
                i += 1
            else:
                j += 1
        elif chrom_a < chrom_b:
            i += 1
        else:
            j += 1

    return overlap_bases

def permute_ranges(ranges, chrom_lengths):
    permuted_ranges = []
    for chrom, start, end in ranges:
        chrom_length = chrom_lengths[chrom]
        range_length = end - start
        new_start = np.random.randint(0, chrom_length - range_length)
        permuted_ranges.append((chrom, new_start, new_start + range_length))
    return permuted_ranges

def permutation_test(ranges_a, ranges_b, chrom_lengths, num_permutations=10000):
    observed_overlap = calculate_overlap(ranges_a, ranges_b)
    permuted_overlaps = []
    
    for _ in range(num_permutations):
        permuted_a = permute_ranges(ranges_a, chrom_lengths)
        overlap = calculate_overlap(permuted_a, ranges_b)
        permuted_overlaps.append(overlap)
    
    p_value = np.sum(np.array(permuted_overlaps) >= observed_overlap) / num_permutations
    return observed_overlap, p_value

def main():
    set_a_file = sys.argv[1]
    set_b_file = sys.argv[2]
    fai_file = sys.argv[3]
    num_permutations = int(sys.argv[4]) if len(sys.argv) > 4 else 10000    
    ranges_a = parse_bed(set_a_file)
    ranges_b = parse_bed(set_b_file)
    chrom_lengths = parse_fai(fai_file)
    observed_overlap, p_value = permutation_test(ranges_a, ranges_b, chrom_lengths, num_permutations)
    
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

if __name__ == "__main__":
    main()
