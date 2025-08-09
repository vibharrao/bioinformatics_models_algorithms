import numpy as np
import argparse

def load_population_schedule(filename):
    """Loads population schedule from a TSV file."""
    pop_schedule = []
    with open(filename, 'r') as file:
        for line in file:
            generation, pop_size = line.strip().split('\t')
            pop_schedule.append((int(generation), int(pop_size)))
    return pop_schedule

def get_population_size(generation, pop_schedule):
    """Gets the population size for the current generation based on the schedule."""
    current_size = pop_schedule[0][1]  # Default to first population size
    for gen, pop_size in pop_schedule:
        if generation < gen:
            break
        current_size = pop_size
    return current_size

def wright_fisher_simulation(allele_freq, pop_size, fitness):
    """Simulates one Wright-Fisher generation with selection."""
    prob_allele = allele_freq * fitness / (allele_freq * fitness + (1 - allele_freq))
    next_gen_count = np.random.binomial(pop_size, prob_allele)
    return next_gen_count / pop_size

def simulate_fixation_loss(allele_freq, pop_schedule, fitness, replicates):
    fixation_generations = []
    loss_generations = []
    
    for _ in range(replicates):
        freq = allele_freq
        generations = 0
        
        while 0 < freq < 1:
            pop_size = get_population_size(generations, pop_schedule)
            freq = wright_fisher_simulation(freq, pop_size, fitness)
            generations += 1
        
        if freq == 1:
            fixation_generations.append(generations)
        else:
            loss_generations.append(generations)
    
    return fixation_generations, loss_generations

def calculate_statistics(generations):
    mean_gen = np.mean(generations) if generations else 0
    variance_gen = np.var(generations) if generations else 0
    return mean_gen, variance_gen

def main():
    parser = argparse.ArgumentParser(description="Wright-Fisher Model Simulation with Complex Demography")
    parser.add_argument("--allele_freq", type=float, required=True, help="Initial allele frequency (0 to 1)")
    parser.add_argument("--pop_size_file", type=str, required=True, help="Path to TSV file with population size schedule")
    parser.add_argument("--fitness", type=float, required=True, help="Relative fitness of the allele")
    parser.add_argument("--replicates", type=int, required=True, help="Number of simulation replicates")
    args = parser.parse_args()
    
    # Load population schedule from TSV file
    pop_schedule = load_population_schedule(args.pop_size_file)
    
    # Run simulations
    fixation_generations, loss_generations = simulate_fixation_loss(
        args.allele_freq, pop_schedule, args.fitness, args.replicates
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
