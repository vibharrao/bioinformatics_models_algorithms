# Assignment 9

## Program Requirements
For all tiers, your program must run in the conda env provided: `BME205_Assn9.yml`. This env only has Python 3.12.4 installed and numpy 2.1.1, so your program must run using only standard Python libraries as well as numpy.

## Tier 1: Implement Viterbi Decoding

**Inbred regions** refer to segments of a genome where an individual has inherited identical alleles from both parents due to mating between closely related individuals. These regions can result from a process known as **inbreeding**, where genetic diversity is reduced due to repeated crosses within a limited gene pool. Identifying inbred regions is crucial for understanding genetic disorders, breeding programs, and evolutionary dynamics, as these regions are more likely to harbor deleterious recessive mutations.

Write a script that accepts a provided VCF file and determines the **most likely set of inbred regions** for each individual using the **Viterbi algorithm**. To do this, you will need to determine the emission probabilities and use the following provided transition probabilities:  

- **Transition from inbred to outbred**: `1 / (1.5 × 10^6)` (per base pair)  
- **Transition from outbred to inbred**: `1 / (4 × 10^6)` (per base pair)  

### Emission Probabilities

| Genotype     | **Inbred State** | **Outbred State** |
|--------------|------------------|-------------------|
| Homozygous   | `1 - e`          | `1 - 2pq`         |
| Heterozygous | `e`              | `2pq`             |

- `e`: Sequencing error rate `e = 1 / 1000`.  
- `p`: Reference allele frequency.  
- `q`: Alternative allele frequency `q = 1 - p`.  

### Input

There is an included VCF file for you to use: `synthetic_population.vcf`

Your script should accept a VCF file as an argument. We will run your script as so:


`python FirstName_LastName_Tier1.py input.vcf`


### Output

After running the Viterbi algorithm, track consecutive inbred states and record the **start** and **stop positions** of each inbred segment in each sample in the VCF, sorted by sample name and start position.

The output should be written to **stdout** in a **tab-delimited format** with the following columns: 
1. indiviual
2. start_position
3. stop_position

For example:
<pre><code>individual    start_position    stop_position
sample1       10               100
sample1       150              300
sample2       50               200
sample2       250              400
sample3       100              250
</code></pre>

## Tier 2: Model Optimization Using Forward Equations

Using the **forward equation** to compute the model likelihoods given the starting parameters above, use the **Amoeba direct search algorithm** to optimize the expected transition rates.  

### Input

Your script should accept a VCF file as an argument. We will run your script as so:


`python FirstName_LastName_Tier2.py input.vcf`


### Output

Your program should output the results to **stdout** in the following format:
<pre><code>P(transition outbred>inbred): [estimate] 
P(transition inbred>outbred): [estimate]</code></pre>

## Tier 3: 

Is there anything unusual about the parts of the genome that are inbred or outbred in many samples? Why might that be the case? Provide you answer in the canvas submission.