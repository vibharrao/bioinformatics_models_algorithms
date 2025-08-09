# Assignment 7: Monte Carlo Simulation of Allele Fixation in a Wright-Fisher Model

## Background

In this assignment, you will create a Python script to model allele frequency dynamics in a Wright-Fisher population of haploids. The goal is to estimate the expected number of generations until an allele either fixes or is lost from the population.

This assignment uses Monte Carlo simulations to estimate fixation time under different initial conditions, population sizes, and selection pressures.

## Description
You will design a Python script called that simulates allele fixation dynamics in a haploid Wright-Fisher population. The script will estimate the expected number of generations until allele fixation and the variance in this estimation across multiple replicate simulations.

## Program Requirements
For all tiers, your program must run in the conda env provided: `BME205_Assn7.yml`. This env only has Python 3.12.4 installed and numpy 2.1.1, so your program must run using only standard Python libraries as well as numpy.

## Tier 1

### Simulation Model:
Implement the Wright-Fisher model with selection. Each generation, simulate the allele frequencies based on the current allele frequency and fitness values, with the following assumptions:

1. Haploid organisms: Each organism contributes one allele.
2. Generations: Assume non-overlapping generations.
3. Selection: Adjust the probability of allele transmission based on the specified fitness value. Use this to simulate natural selection, affecting the probability of allele inheritance in the next generation.

The simulation should track generations until:

The allele either reaches fixation (frequency of 1) or loss (frequency of 0).
Each replicate simulation will yield the number of generations required for either fixation or loss.

### Command-Line Arguments:
Your script should accept the following command-line arguments:

1. `--allele_freq`: The initial frequency of the allele (float between 0 and 1).
2. `--pop_size`: The population size (integer). Since it is a haploid population, this is the number of individuals, each contributing one allelee.
3. `--fitness`: The relative fitness of the allele (float, where 1 represents neutral fitness, values >1 indicate positive selection, and values <1 indicate negative selection).
4. `--replicates`: The number of Monte Carlo simulation replicates (integer).

*Use argparse to handle the command-line arguments.*

### Output:
Your script should output:

1. Expected number of generations until fixation or loss (average across replicates, separate for each outcome).
2. Variance in fixation time across replicates.

Output Format: Print the following to the console:

**If the allele was fixed:**

`Allele was fixed in <mean_generations>. Variance: <variance_generations>`

**OR, if the allele was lost:**

`Allele was lost in <mean_generations>. Variance: <variance_generations>`

Where `<mean_generations>` and `<variance_generations>` are the values your script computed from the simulations.


### Expected Execution

`python Firstname_Lastname_assignment7_Tier1.py --allele_freq 0.1 --pop_size 100 --fitness 1.05 --replicates 1000`

## Tier 2: Complex Demography

What if the population size is not constant through time but instead changes instantaneously at particular generations?  
Add a feature to your scriot that can accept a tsv that indicates a generation and population size.  
The population will remain that size until the next generation indicated or until one allele goes to fixation.

### Population sizes
Where popsize.tsv is a tsv that has two elements per line: generation and population size. For example:
```
0 100
50 1000
100 200
```
This population file indicates that from generation 0 to 50, the population size is 100. Then at generation 50 the population size instantaneously increases to 1000. Then at generation 100 and onwards the population size is reduced to 200.

We have included a popsize.tsv to use. You can expect that this file will be present in the working directory that your submitted script is run.

### Output
Same as tier 1.

### Expected Execution

```python Firstname_Lastname_assignment7_Tier2.py --allele_freq 0.1 --pop_size_file popsize.tsv --fitness 1.05 --replicates 1000 ```

 
## Tier 3: Do it, but backwards. 

It is also possible to record only the backward-in-time processes that determine the ancestry of a current sample using Coalescent theory.

Use a coalescent simulation that can determine the expected time and variance of the **eighth** coalescent event for a sample size (`--sample_size`) to be defined via the command line within a larger population, (`--pop_size`). Time should be in units generations. The following command line is needed:

### Output
The output should be formatted like so:

`Time to eighth coalescent event: <mean_generations>. Variance: <variance_generations>`


### Expected execution
```python Firstname_Lastname_assignment7_Tier3.py --pop_size 100000 --sample_size 10 --replicates 1000```


