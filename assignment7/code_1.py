import argparse
import numpy as np

def simulate_wright_fisher(allele_freq, pop_size, fitness):
    generations = 0
    while 0 < allele_freq < 1:
        generations += 1
        # Calculate the expected frequency with selection
        p = allele_freq * fitness / ((allele_freq * fitness) + (1 - allele_freq))
        
        # Sample the number of alleles in the next generation
        allele_count = np.random.binomial(pop_size, p)
        allele_freq = allele_count / pop_size

    return generations, allele_freq

def main():
    parser = argparse.ArgumentParser(description="Wright-Fisher model with selection")
    parser.add_argument("--allele_freq", type=float, required=True, help="Initial allele frequency")
    parser.add_argument("--pop_size", type=int, required=True, help="Population size (number of haploid individuals)")
    parser.add_argument("--fitness", type=float, required=True, help="Relative fitness of the allele")
    parser.add_argument("--replicates", type=int, required=True, help="Number of simulation replicates")

    args = parser.parse_args()

    fixation_times = []
    loss_times = []

    for _ in range(args.replicates):
        generations, final_freq = simulate_wright_fisher(args.allele_freq, args.pop_size, args.fitness)
        
        if final_freq == 1:
            fixation_times.append(generations)
        elif final_freq == 0:
            loss_times.append(generations)

    if fixation_times:
        mean_fixation = np.mean(fixation_times)
        var_fixation = np.var(fixation_times)
        print(f"Allele was fixed in {mean_fixation:.2f} generations. Variance: {var_fixation:.2f}")

    if loss_times:
        mean_loss = np.mean(loss_times)
        var_loss = np.var(loss_times)
        print(f"Allele was lost in {mean_loss:.2f} generations. Variance: {var_loss:.2f}")

if __name__ == "__main__":
    main()
