import numpy as np
import sys

# Constants
SEQUENCING_ERROR_RATE = 1 / 1000  # e
INITIAL_P_OUTBRED_TO_INBRED = 1 / (4 * 10**6)
INITIAL_P_INBRED_TO_OUTBRED = 1 / (1.5 * 10**6)
PERTURBATION = 1e-7  # Perturbation for simplex initialization
MAX_ITER = 500  # Maximum iterations for Nelder-Mead
TOL = 1e-6  # Convergence tolerance

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

def calculate_reference_allele_frequencies(genotypes, individuals):
    """Calculate reference allele frequencies for each position."""
    ref_freqs = []
    for genotype_data in genotypes:
        ref_count = 0
        total_count = 0
        for genotype in genotype_data[2]:  # Iterate over all individuals
            alleles = genotype.split("|")
            ref_count += alleles.count("0")  # Count reference alleles
            total_count += len(alleles)  # Total alleles
        ref_freqs.append(ref_count / total_count if total_count > 0 else 0.5)  # Avoid division by zero
    return ref_freqs

def compute_emission_probabilities(ref_freq, genotype, inbred=True):
    """Computes emission probabilities for inbred or outbred states."""
    alt_freq = 1 - ref_freq
    if genotype == "1|1" or genotype == "0|0":  # Homozygous
        if inbred:
            return max(1 - SEQUENCING_ERROR_RATE, 1e-10)
        else:
            return max(1 - 2 * ref_freq * alt_freq, 1e-10)
    elif genotype == "1|0" or genotype == "0|1":  # Heterozygous
        if inbred:
            return max(SEQUENCING_ERROR_RATE, 1e-10)
        else:
            return max(2 * ref_freq * alt_freq, 1e-10)
    else:
        return 1e-10  # Small fallback probability

def forward_algorithm(positions, genotypes, ref_freqs, p_outbred_to_inbred, p_inbred_to_outbred):
    """Calculate the likelihood of the model using the forward algorithm."""
    num_positions = len(positions)
    alpha = np.zeros((2, num_positions))  # Forward probabilities for inbred (0) and outbred (1)

    # Initial probabilities
    alpha[0, 0] = 0.5 * compute_emission_probabilities(ref_freqs[0], genotypes[0], inbred=True)  # Inbred
    alpha[1, 0] = 0.5 * compute_emission_probabilities(ref_freqs[0], genotypes[0], inbred=False)  # Outbred

    for t in range(1, num_positions):
        distance = positions[t] - positions[t-1]
        p_inbred_to_inbred = 1 - np.exp(-distance / (1 / p_outbred_to_inbred))
        p_outbred_to_outbred = 1 - np.exp(-distance / (1 / p_inbred_to_outbred))

        # Emission probabilities
        e_inbred = compute_emission_probabilities(ref_freqs[t], genotypes[t], inbred=True)
        e_outbred = compute_emission_probabilities(ref_freqs[t], genotypes[t], inbred=False)

        # Update forward probabilities
        alpha[0, t] = (alpha[0, t-1] * p_inbred_to_inbred + alpha[1, t-1] * p_outbred_to_inbred) * e_inbred
        alpha[1, t] = (alpha[0, t-1] * p_inbred_to_outbred + alpha[1, t-1] * p_outbred_to_outbred) * e_outbred

    # Ensure sum is greater than zero before applying log
    final_sum = np.sum(alpha[:, -1])
    if final_sum <= 0:
        return -np.inf  # Log of zero or negative probability

    return np.log(final_sum)  # Log likelihood


def nelder_mead(objective_function, initial_guess):
    """Manual implementation of the Nelder-Mead algorithm."""
    simplex = np.array([
        initial_guess,
        initial_guess + np.array([PERTURBATION, 0]),
        initial_guess + np.array([0, PERTURBATION])
    ])
    scores = np.array([objective_function(p) for p in simplex])

    for iteration in range(MAX_ITER):
        order = np.argsort(scores)
        simplex = simplex[order]
        scores = scores[order]

        # Centroid of the best points
        centroid = np.mean(simplex[:-1], axis=0)

        # Reflection
        reflected = centroid + (centroid - simplex[-1])
        reflected_score = objective_function(reflected)
        if scores[0] <= reflected_score < scores[-2]:
            simplex[-1] = reflected
            scores[-1] = reflected_score
            continue

        # Expansion
        if reflected_score < scores[0]:
            expanded = centroid + 2 * (centroid - simplex[-1])
            expanded_score = objective_function(expanded)
            if expanded_score < reflected_score:
                simplex[-1] = expanded
                scores[-1] = expanded_score
            else:
                simplex[-1] = reflected
                scores[-1] = reflected_score
            continue

        # Contraction
        contracted = centroid + 0.5 * (simplex[-1] - centroid)
        contracted_score = objective_function(contracted)
        if contracted_score < scores[-1]:
            simplex[-1] = contracted
            scores[-1] = contracted_score
            continue

        # Shrink
        simplex[1:] = simplex[0] + 0.5 * (simplex[1:] - simplex[0])
        scores[1:] = [objective_function(p) for p in simplex[1:]]

        # Convergence check
        if np.max(np.abs(simplex[1:] - simplex[0])) < TOL:
            break

    return simplex[0]

def optimize_with_nelder_mead(positions, genotypes, ref_freqs):
    """Optimize transition probabilities using manual Nelder-Mead."""
    def objective_function(params):
        p_outbred_to_inbred, p_inbred_to_outbred = params
        if p_outbred_to_inbred <= 0 or p_inbred_to_outbred <= 0:
            return np.inf  # Penalize invalid probabilities
        return -forward_algorithm(positions, genotypes, ref_freqs, p_outbred_to_inbred, p_inbred_to_outbred)
    
    initial_guess = [INITIAL_P_OUTBRED_TO_INBRED, INITIAL_P_INBRED_TO_OUTBRED]
    return nelder_mead(objective_function, initial_guess)

def main(vcf_file):
    positions, individuals, genotypes = parse_vcf(vcf_file)
    ref_freqs = calculate_reference_allele_frequencies(genotypes, individuals)

    # Select one individual for simplicity
    individual_idx = 0
    indiv_genotypes = [g[2][individual_idx] for g in genotypes]

    # Optimize transition probabilities
    p_outbred_to_inbred, p_inbred_to_outbred = optimize_with_nelder_mead(positions, indiv_genotypes, ref_freqs)

    # Output results
    print(f"P(transition outbred>inbred): {p_outbred_to_inbred}")
    print(f"P(transition inbred>outbred): {p_inbred_to_outbred}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.vcf")
        sys.exit(1) 
    main(sys.argv[1])