import numpy as np
import argparse

def wright_fisher_simulation(allele_freq, pop_size, fitness):
    """Simulates one Wright-Fisher generation with selection."""
    # Calculate adjusted fitness-based probability for allele inheritance
    prob_allele = allele_freq * fitness / (allele_freq * fitness + (1 - allele_freq))
    # Generate the next generation's allele frequency by binomial sampling
    next_gen_count = np.random.binomial(pop_size, prob_allele)
    return next_gen_count / pop_size

def simulate_fixation_loss(allele_freq, pop_size, fitness, replicates):
    fixation_generations = []
    loss_generations = []
    
    for _ in range(replicates):
        freq = allele_freq
        generations = 0
        
        # Run simulation until allele fixation (freq = 1) or loss (freq = 0)
        while 0 < freq < 1:
            freq = wright_fisher_simulation(freq, pop_size, fitness)
            generations += 1
            
        
        # Record results depending on fixation or loss
        if freq == 1:
            fixation_generations.append(generations)
        else:
            loss_generations.append(generations)
    
    return fixation_generations, loss_generations

def calculate_statistics(generations):
    """Calculate mean and variance for the generations list."""
    mean_gen = np.mean(generations) if generations else 0
    variance_gen = np.var(generations) if generations else 0
    return mean_gen, variance_gen

def main():
    parser = argparse.ArgumentParser(description="Wright-Fisher Model Simulation")
    parser.add_argument("--allele_freq", type=float, required=True, help="Initial allele frequency (0 to 1)")
    parser.add_argument("--pop_size", type=int, required=True, help="Population size (integer)")
    parser.add_argument("--fitness", type=float, required=True, help="Relative fitness of the allele")
    parser.add_argument("--replicates", type=int, required=True, help="Number of simulation replicates")
    args = parser.parse_args()
    
    # Run simulations
    fixation_generations, loss_generations = simulate_fixation_loss(
        args.allele_freq, args.pop_size, args.fitness, args.replicates
    )
    
    # Calculate statistics
    if fixation_generations:
        mean_fix, var_fix = calculate_statistics(fixation_generations)
        print(f"Allele was fixed in {mean_fix:.2f} generations. Variance: {var_fix:.2f}")
    
    if loss_generations:
        mean_loss, var_loss = calculate_statistics(loss_generations)
        print(f"Allele was lost in {mean_loss:.2f} generations. Variance: {var_loss:.2f}")

if __name__ == "__main__":
    main()
