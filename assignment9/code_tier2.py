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
        ref_freqs.append(max(ref_count / total_count, 1e-10) if total_count > 0 else 0.5)  # Avoid division by zero
    return ref_freqs

def compute_emission_probabilities(ref_freq, genotype, inbred=True):
    """Computes emission probabilities for inbred or outbred states with added numerical stability."""
    alt_freq = 1 - ref_freq
    SEQUENCING_ERROR_RATE_STABLE = max(SEQUENCING_ERROR_RATE, 1e-10)
    
    if genotype in {"1|1", "0|0"}:  # Homozygous
        if inbred:
            return max(1 - SEQUENCING_ERROR_RATE_STABLE, 1e-10)
        else:
            return max(1 - 2 * ref_freq * alt_freq, 1e-10)
    elif genotype in {"1|0", "0|1"}:  # Heterozygous
        if inbred:
            return max(SEQUENCING_ERROR_RATE_STABLE, 1e-10)
        else:
            return max(2 * ref_freq * alt_freq, 1e-10)
    else:
        return 1e-10  # Default fallback

def forward_algorithm(positions, genotypes, ref_freqs, p_outbred_to_inbred, p_inbred_to_outbred):
    """Calculate the likelihood of the model using a log-scaled forward algorithm."""
    num_positions = len(positions)
    if num_positions == 1:
        return np.log(compute_emission_probabilities(ref_freqs[0], genotypes[0]))

    log_alpha = np.full((2, num_positions), -np.inf)  # Log-scaled forward probabilities (in log space)

    # Initial probabilities (in log-space)
    log_alpha[0, 0] = np.log(0.5) + np.log(compute_emission_probabilities(ref_freqs[0], genotypes[0], inbred=True))
    log_alpha[1, 0] = np.log(0.5) + np.log(compute_emission_probabilities(ref_freqs[0], genotypes[0], inbred=False))

    for t in range(1, num_positions):
        distance = positions[t] - positions[t - 1]
        p_inbred_to_inbred = 1 - np.exp(-distance / (1 / p_outbred_to_inbred))
        p_outbred_to_outbred = 1 - np.exp(-distance / (1 / p_inbred_to_outbred))
        
        # Add numerical stability to transition probabilities
        p_inbred_to_inbred = np.clip(p_inbred_to_inbred, 1e-10, 1)
        p_outbred_to_outbred = np.clip(p_outbred_to_outbred, 1e-10, 1)
        p_outbred_to_inbred = np.clip(1 - p_inbred_to_inbred, 1e-10, 1)
        p_inbred_to_outbred = np.clip(1 - p_outbred_to_outbred, 1e-10, 1)

        # Log transition probabilities
        log_p_inbred_to_inbred = np.log(p_inbred_to_inbred)
        log_p_outbred_to_outbred = np.log(p_outbred_to_outbred)
        log_p_outbred_to_inbred = np.log(p_outbred_to_inbred)
        log_p_inbred_to_outbred = np.log(p_inbred_to_outbred)

        # Emission probabilities in log-space
        log_e_inbred = np.log(compute_emission_probabilities(ref_freqs[t], genotypes[t], inbred=True))
        log_e_outbred = np.log(compute_emission_probabilities(ref_freqs[t], genotypes[t], inbred=False))

        # Update log-alpha
        log_alpha[0, t] = log_e_inbred + np.logaddexp(
            log_alpha[0, t - 1] + log_p_inbred_to_inbred,
            log_alpha[1, t - 1] + log_p_outbred_to_inbred
        )
        log_alpha[1, t] = log_e_outbred + np.logaddexp(
            log_alpha[0, t - 1] + log_p_inbred_to_outbred,
            log_alpha[1, t - 1] + log_p_outbred_to_outbred
        )

    # Final log-sum to compute log-likelihood
    log_likelihood = np.logaddexp(log_alpha[0, -1], log_alpha[1, -1])
    return log_likelihood


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
