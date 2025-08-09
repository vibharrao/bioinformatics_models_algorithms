import sys
import numpy as np

# Constants for transition probabilities
P_INBRED_TO_OUTBRED = 1 / (1.5 * 10**6)
P_OUTBRED_TO_INBRED = 1 / (4 * 10**6)

def parse_vcf(file_path):
    """Parses the VCF file to extract positions and genotypes."""
    positions = []
    individuals = []
    genotypes = []

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("##"):
                continue  # Skip metadata lines
            elif line.startswith("#CHROM"):
                # Extract individual names from header line
                individuals = line.strip().split("\t")[9:]
            else:
                fields = line.strip().split("\t")
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                genotype_data = fields[9:]
                positions.append(pos)
                genotypes.append((ref, alt, genotype_data))
    
    return positions, individuals, genotypes

def compute_emission_probabilities(ref_freq, genotype, inbred=True):
    """Computes emission probabilities for inbred or outbred states."""
    alt_freq = 1 - ref_freq
    SEQUENCING_ERROR_RATE = 1 / 1000  # Sequencing error rate
    if genotype == "1|1" or genotype == "0|0":  # Homozygous
        if inbred:
            return 1 - SEQUENCING_ERROR_RATE
        else:
            return 1 - 2 * ref_freq * alt_freq
    elif genotype == "1|0" or genotype == "0|1":  # Heterozygous
        if inbred:
            return SEQUENCING_ERROR_RATE
        else:
            return 2 * ref_freq * alt_freq
    else:
        return 0.0  # Invalid genotype

def viterbi_algorithm(positions, genotypes, ref_freqs):
    """Applies the Viterbi algorithm to determine inbred regions."""
    num_positions = len(positions)
    num_states = 2  # Inbred, Outbred

    # Initialize Viterbi variables
    viterbi_log_probs = np.full((num_states, num_positions), -np.inf)
    backtrack = np.zeros((num_states, num_positions), dtype=int)

    # Initial probabilities
    viterbi_log_probs[0, 0] = np.log(0.5)  # Inbred
    viterbi_log_probs[1, 0] = np.log(0.5)  # Outbred

    # Transition probabilities
    log_p_inbred_to_inbred = np.log(1 - P_INBRED_TO_OUTBRED)
    log_p_inbred_to_outbred = np.log(P_INBRED_TO_OUTBRED)
    log_p_outbred_to_inbred = np.log(P_OUTBRED_TO_INBRED)
    log_p_outbred_to_outbred = np.log(1 - P_OUTBRED_TO_INBRED)

    for t in range(1, num_positions):
        for current_state in range(num_states):
            max_log_prob = -np.inf
            max_state = None
            for prev_state in range(num_states):
                log_transition_prob = (
                    log_p_outbred_to_inbred if prev_state == 1 and current_state == 0
                    else log_p_inbred_to_outbred if prev_state == 0 and current_state == 1
                    else log_p_inbred_to_inbred if prev_state == 0 and current_state == 0
                    else log_p_outbred_to_outbred
                )
                log_prob = viterbi_log_probs[prev_state, t-1] + log_transition_prob
                if log_prob > max_log_prob:
                    max_log_prob = log_prob
                    max_state = prev_state

            genotype = genotypes[t]
            ref_freq = ref_freqs[t]
            emission_prob = compute_emission_probabilities(ref_freq, genotype, inbred=(current_state == 0))
            log_emission_prob = np.log(emission_prob) if emission_prob > 0 else -np.inf
            viterbi_log_probs[current_state, t] = max_log_prob + log_emission_prob
            backtrack[current_state, t] = max_state

    # Backtrack to find most likely states
    most_likely_states = []
    current_state = np.argmax(viterbi_log_probs[:, -1])  # Last position's most likely state
    for t in range(num_positions-1, -1, -1):
        most_likely_states.append(current_state)
        current_state = backtrack[current_state, t]

    return most_likely_states[::-1]

def find_inbred_segments(individual, positions, states):
    """Finds consecutive inbred segments based on state path."""
    inbred_segments = []
    start = None

    for i, state in enumerate(states):
        if state == 0:  # Inbred state
            if start is None:
                start = positions[i]
        else:
            if start is not None:
                inbred_segments.append((individual, start, positions[i-1]))
                start = None

    if start is not None:
        inbred_segments.append((individual, start, positions[-1]))

    return inbred_segments

def main(vcf_file,output_file):
    positions, individuals, genotypes = parse_vcf(vcf_file)
    ref_freqs = [0.5] * len(positions)

    results = []
    for i, individual in enumerate(individuals):
        indiv_genotypes = [g[2][i] for g in genotypes]
        states = viterbi_algorithm(positions, indiv_genotypes, ref_freqs)
        inbred_segments = find_inbred_segments(individual, positions, states)
        results.extend(inbred_segments)

#     print("individual\tstart_position\tstop_position")
#     for individual, start, stop in results:
#         print(f"{individual}\t{start}\t{stop}")

# if __name__ == "__main__":
#     if len(sys.argv) != 2:
#         print("Usage: python script.py input.vcf")
#         sys.exit(1)
#     main(sys.argv[1])

    with open(output_file, 'w') as f:
        f.write("individual\tstart_position\tstop_position\n")
        for individual, start, stop in results:
            f.write(f"{individual}\t{start}\t{stop}\n")
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.vcf output.txt")
        sys.exit(1)
    vcf_file = sys.argv[1]
    output_file = sys.argv[2]
    main(vcf_file, output_file)