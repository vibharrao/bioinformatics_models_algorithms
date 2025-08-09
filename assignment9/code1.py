import sys
import numpy as np

# Constants
e = 1 / 1000  # Sequencing error rate
TRANS_INBRED_TO_OUTBRED = 1 / (1.5 * 10**6)
TRANS_OUTBRED_TO_INBRED = 1 / (4 * 10**6)
MIN_REGION_LENGTH = 10  # Minimum length for inbred regions to be included

def parse_vcf(vcf_file):
    positions = []  # List of positions
    genotypes = {}  # Genotype information for each individual
    allele_counts = {}  # Reference and alternate allele counts at each position

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):  # Skip header lines
                if line.startswith("#CHROM"):
                    samples = line.strip().split("\t")[9:]  # Extract sample names
                continue

            fields = line.strip().split("\t")
            pos = int(fields[1])
            positions.append(pos)

            # Initialize allele counts
            ref_count = 0
            alt_count = 0

            # Extract genotypes for all samples
            for i, genotype in enumerate(fields[9:]):
                if samples[i] not in genotypes:
                    genotypes[samples[i]] = []
                genotypes[samples[i]].append(genotype.split(":")[0])  # Store only the genotype

                # Count alleles
                if "0" in genotype:
                    ref_count += genotype.count("0")
                if "1" in genotype:
                    alt_count += genotype.count("1")

            allele_counts[pos] = (ref_count, alt_count)  # Store allele counts for the position

    return positions, genotypes, allele_counts, samples

def calculate_emission_probabilities(genotype, p, q):
    if p < 0 or q < 0 or p + q != 1:
        return None, None

    if genotype in {"0|0", "0/0"}:  # Homozygous reference
        return 1 - e, 1 - 2 * p * q
    elif genotype in {"1|1", "1/1"}:  # Homozygous alternate
        return 1 - e, 1 - 2 * p * q
    elif genotype in {"0|1", "1|0", "0/1", "1/0"}:  # Heterozygous
        return e, 2 * p * q
    elif genotype in {"./.", ".|."}:  # Missing data
        return None, None
    else:
        return None, None

def viterbi_algorithm(positions, genotypes, allele_counts, samples):
    results = []
    for sample in samples:
        sample_genotypes = genotypes[sample]
        n = len(sample_genotypes)
        
        # Initialize Viterbi tables
        dp_inbred = np.zeros(n)
        dp_outbred = np.zeros(n)
        prev_inbred = np.zeros(n, dtype=int)
        prev_outbred = np.zeros(n, dtype=int)

        # Start probabilities
        dp_inbred[0] = 0.5
        dp_outbred[0] = 0.5

        # Fill DP tables
        for i in range(1, n):
            ref_count, alt_count = allele_counts[positions[i]]
            total = ref_count + alt_count
            if total == 0:
                continue  # Skip positions without allele counts
            p = ref_count / total
            q = 1 - p
            em_inbred, em_outbred = calculate_emission_probabilities(sample_genotypes[i], p, q)
            
            # Skip positions with invalid emission probabilities
            if em_inbred is None or em_outbred is None:
                continue

            dp_inbred[i] = max(
                dp_inbred[i - 1] * (1 - TRANS_INBRED_TO_OUTBRED) * em_inbred,
                dp_outbred[i - 1] * TRANS_OUTBRED_TO_INBRED * em_inbred
            )
            dp_outbred[i] = max(
                dp_outbred[i - 1] * (1 - TRANS_OUTBRED_TO_INBRED) * em_outbred,
                dp_inbred[i - 1] * TRANS_INBRED_TO_OUTBRED * em_outbred
            )
            prev_inbred[i] = 1 if dp_inbred[i - 1] * (1 - TRANS_INBRED_TO_OUTBRED) > \
                                 dp_outbred[i - 1] * TRANS_OUTBRED_TO_INBRED else 0
            prev_outbred[i] = 1 if dp_outbred[i - 1] * (1 - TRANS_OUTBRED_TO_INBRED) > \
                                  dp_inbred[i - 1] * TRANS_INBRED_TO_OUTBRED else 0

        # Traceback
        in_inbred = dp_inbred[-1] > dp_outbred[-1]
        state = [in_inbred]
        for i in range(n - 2, -1, -1):
            state.append(prev_inbred[i + 1] if state[-1] else prev_outbred[i + 1])
        state = state[::-1]

        # Identify and merge consecutive inbred regions
        start = None
        for i, s in enumerate(state):
            if s:
                if start is None:
                    start = positions[i]
            else:
                if start is not None:
                    results.append((sample, start, positions[i - 1]))
                    start = None
        if start is not None:
            results.append((sample, start, positions[-1]))

    # Filter results by minimum region length
    filtered_results = [
        (sample, start, end) for sample, start, end in results if end - start + 1 >= MIN_REGION_LENGTH
    ]

    return filtered_results

def main():
    if len(sys.argv) != 2:
        print("Usage: python code.py input.vcf")
        return

    vcf_file = sys.argv[1]
    positions, genotypes, allele_counts, samples = parse_vcf(vcf_file)
    results = viterbi_algorithm(positions, genotypes, allele_counts, samples)

    # Output results
    print("individual\tstart_position\tstop_position")
    for result in results:
        print("\t".join(map(str, result)))

if __name__ == "__main__":
    main()
