import numpy as np
import argparse

def simulate_coalescent_time(pop_size, sample_size, event_number):
    """
    Simulate the time to a specific coalescent event for a sample size within a population.
    Returns the time in generations to reach the specified coalescent event.
    """
    k = sample_size
    total_time = 0
    coalescent_events = 0
    
    while coalescent_events < event_number and k > 1:
        # Mean time to the next coalescent event
        time_to_next_event = np.random.exponential(2 * pop_size / (k * (k - 1)))
        total_time += time_to_next_event
        coalescent_events += 1
        k -= 1  # Reduce the sample size as lineages coalesce
    
    return total_time if coalescent_events == event_number else None

def simulate_replicates(pop_size, sample_size, event_number, replicates):
    """
    Perform multiple coalescent simulations and collect the time to the specified coalescent event.
    Returns a list of times for each replicate.
    """
    coalescent_times = []
    
    for _ in range(replicates):
        time = simulate_coalescent_time(pop_size, sample_size, event_number)
        if time is not None:
            coalescent_times.append(time)
    
    return coalescent_times

def calculate_statistics(times):
    mean_time = np.mean(times) if times else 0
    variance_time = np.var(times) if times else 0
    return mean_time, variance_time

def main():
    parser = argparse.ArgumentParser(description="Coalescent Simulation to the Eighth Event")
    parser.add_argument("--pop_size", type=int, required=True, help="Population size (integer)")
    parser.add_argument("--sample_size", type=int, required=True, help="Sample size (integer)")
    parser.add_argument("--replicates", type=int, required=True, help="Number of simulation replicates")
    args = parser.parse_args()
    
    # Set the event number for the eighth coalescent event
    event_number = 8
    
    # Run coalescent simulations for the specified number of replicates
    coalescent_times = simulate_replicates(args.pop_size, args.sample_size, event_number, args.replicates)
    
    # Calculate statistics and print the results
    if coalescent_times:
        mean_time, variance_time = calculate_statistics(coalescent_times)
        print(f"Time to eighth coalescent event: {mean_time:.2f} generations. Variance: {variance_time:.2f}")
    else:
        print("Eighth coalescent event not reached in any replicate.")

if __name__ == "__main__":
    main()
