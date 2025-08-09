# Bioinformatics Models & Algorithms

This repository contains implementations for various assignments from **BME 205 — Bioinformatics Models & Algorithms**, completed in Python (with NumPy where allowed).  
Each assignment models a core computational biology concept, from sequence alignment and phylogenetics to stochastic population models and Hidden Markov Models.

---

## Table of Contents
- [Assignment 2 — Sequence Alignment with Affine Gap Penalties](#assignment-2--sequence-alignment-with-affine-gap-penalties)
- [Assignment 3 — Profile HMM Construction & Sequence Scoring](#assignment-3--profile-hmm-construction--sequence-scoring)
- [Assignment 4 — Phylogenetic Tree Reconstruction](#assignment-4--phylogenetic-tree-reconstruction)
- [Assignment 5 — Gibbs Sampling for Motif Discovery](#assignment-5--gibbs-sampling-for-motif-discovery)
- [Assignment 6 — Expectation Maximization for HMM Parameters](#assignment-6--expectation-maximization-for-hmm-parameters)
- [Assignment 7 — Monte Carlo Simulation of Allele Fixation](#assignment-7--monte-carlo-simulation-of-allele-fixation)
- [Assignment 9 — HMM for Inbred Region Detection](#assignment-9--hmm-for-inbred-region-detection)

---

## Assignment 2 — Sequence Alignment with Affine Gap Penalties

**Goal:**  
Implement **global and local sequence alignment** algorithms supporting **affine gap penalties** using dynamic programming.

**Key Features:**
- **Affine gap model**: Gap opening and gap extension penalties treated separately.
- Supports both **Needleman–Wunsch** (global) and **Smith–Waterman** (local) modes.
- Outputs:
  - Alignment score
  - Aligned sequences
  - Traceback path

**Run:**
```bash
python assignment2.py seq1.fasta seq2.fasta --gap_open -5 --gap_extend -1 --mode global
```
## Assignment 3 — Profile HMM Construction & Sequence Scoring
**Goal:**
Construct a Profile Hidden Markov Model (Profile-HMM) from a multiple sequence alignment and score new sequences against it.

**Key Features:**

- Estimate match, insert, delete state probabilities from the MSA.

- Handle pseudocounts for zero-frequency correction.

- Viterbi decoding to find most likely state path.

- Sequence scoring to determine likelihood under the profile.

**Run:**

```bash
python assignment3.py alignment.fasta test_sequences.fasta
```
## Assignment 4 — Phylogenetic Tree Reconstruction

**Goal:**  
Infer a phylogenetic tree from sequence data using distance-based and/or maximum likelihood methods.

**Key Features:**
- Pairwise distance matrix computation (p-distance, Jukes–Cantor).
- Neighbor-Joining algorithm for tree construction.
- Outputs Newick-formatted tree.
- Optional tree visualization using `matplotlib` or `ete3`.

**Run Example:**
```bash
python assignment4.py sequences.fasta --method neighbor_joining
```

## Assignment 5 — Gibbs Sampling for Motif Discovery

**Goal:**  
Implement Gibbs sampling to identify shared sequence motifs across multiple DNA sequences.

**Features:**
- Random motif initialization.
- Iterative sampling to refine motif position in each sequence.
- Motif profile updates excluding current sequence.

**Outputs:**
- Discovered motif consensus.
- Motif position in each sequence.
- Final Position Weight Matrix (PWM).

**Run Example:**
```bash
python assignment5.py sequences.fasta --motif_length 8 --iterations 1000
```

**Example Output:**
```
Consensus Motif: ACGTAGCT
Motif Positions:
seq1: position 5
seq2: position 3
seq3: position 7

Final PWM:
A: 0.80 0.10 0.00 0.70 0.20 0.60 0.10 0.90
C: 0.10 0.70 0.20 0.10 0.60 0.20 0.70 0.05
G: 0.05 0.10 0.70 0.10 0.10 0.10 0.10 0.02
T: 0.05 0.10 0.10 0.10 0.10 0.10 0.10 0.03
```

---

## Assignment 6 — Expectation Maximization for HMM Parameters

**Goal:**  
Use the Baum–Welch algorithm (EM for HMMs) to learn model parameters from observed sequences.

**Features:**
- Implements forward–backward probability calculation.
- Updates transition and emission probabilities until convergence.
- Handles multiple sequences.

**Outputs:**
- Final transition matrix.
- Final emission matrix.
- Log-likelihood progression.

**Run Example:**
```bash
python assignment6.py sequences.txt --states 3 --alphabet A C G T
```

---

## Assignment 7 — Monte Carlo Simulation of Allele Fixation

**Goal:**  
Simulate allele frequency dynamics in a haploid Wright–Fisher population to estimate fixation/loss times under selection.

**Tiers:**
- Tier 1: Basic simulation with constant population size.
- Tier 2: Variable demography from a population-size change file.
- Tier 3: Backward-time coalescent simulation to estimate time to 8th coalescent event.

**Run Examples:**
```bash
# Tier 1
python Firstname_Lastname_assignment7_Tier1.py --allele_freq 0.1 --pop_size 100 --fitness 1.05 --replicates 1000

# Tier 2
python Firstname_Lastname_assignment7_Tier2.py --allele_freq 0.1 --pop_size_file popsize.tsv --fitness 1.05 --replicates 1000

# Tier 3
python Firstname_Lastname_assignment7_Tier3.py --pop_size 100000 --sample_size 10 --replicates 1000
```

**Example Output:**
```
Allele was fixed in 250.3. Variance: 15.6
```

---

## Assignment 9 — HMM for Inbred Region Detection

**Goal:**  
Identify inbred genomic segments per individual using a Hidden Markov Model with Viterbi decoding (Tier 1) and optimize transition rates with the forward algorithm + Amoeba search (Tier 2).

**Model:**
- States: Inbred, Outbred.
- Emissions: Based on genotype (homozygous/heterozygous) and allele frequency.
- Transitions:
  - P(I→O) = 1 / (1.5 × 10^6) per bp.
  - P(O→I) = 1 / (4 × 10^6) per bp.

**Run Examples:**
```bash
# Tier 1 — Viterbi decoding
python FirstName_LastName_Tier1.py synthetic_population.vcf

# Tier 2 — Optimize transitions
python FirstName_LastName_Tier2.py synthetic_population.vcf
```

**Example Output (Tier 1):**
```
individual    start_position    stop_position
sample1       10               100
sample1       150              300
sample2       50               200
```

**Example Output (Tier 2):**
```
P(transition outbred>inbred): 2.5e-07
P(transition inbred>outbred): 6.7e-07
```
