# Assignment 3: Permutation Test for Overlaps

## Assignment Description

In this assignment, you will implement a permutation test to analyze the overlap between two sets of ranges. This will help you understand the concept of permutation tests and how they can be applied in practice.

## Learning Objectives

1. Understand the principles of permutation tests.
2. Implement a permutation test algorithm in Python.
3. Analyze and interpret the results of permutation tests.

## Scenario

In your research lab, you are studying the genomic regions associated with a particular trait in a plant species. You have two sets of genomic ranges: one set corresponds to regions of the genome that are associated with high expression levels of a particular gene under drought conditions (Set A), and the other set corresponds to regions associated with the binding sites of a transcription factor known to be involved in drought response (Set B).

You suspect that the transcription factor binding sites (Set B) might play a role in regulating the genes that are highly expressed under drought conditions (Set A). To test this hypothesis, you need to determine whether the number of overlapping bases between these two sets of ranges is statistically significant.

*To achieve this, you will perform a permutation test, which will help you understand if the observed overlap is greater than what would be expected by chance.*

## Instructions

### Tier 1 (requried for full credit)

You are provided with two sets of genomic ranges in bed files: `SetA.bed` and `SetB.bed`. Additionally, there is a fasta index `genome.fa.fai` that describes the length of the chromosome. These files can be found in the directory `Tier_1_Files`.  
Write a python program that implements a permutation test to determine if the *number of overlapping bases* between the genomic ranges is statistically significant.

### Tier 2 (extra credit)

You are provided with two sets of genomic ranges in bed files: `SetA.bed` and `SetB.bed`. Additionally, there is a fasta index `genome.fa.fai` that describes the lengths of each of the chromosomes. These files can be found in the directory `Tier_2_Files`. You will notice that there are *multiple chromosomes* represented in the bed files. If ignored, and you redistribute the ranges at random, you suspect it will cause an inflation in the p-value. Implement a permutation test that accounts for this structure in your data.  

## Program Requirements

1. You must name your python file as `Firstname_Lastname_TierX.py` for each tier you submit. As the programs outputs to standard out, you do not need to submit an output file, just your python scripts.

2. Your program should take 3 required command line arguments, and 1 optional argument in this order:  
  `$ python submission.py path/to/SetA.bed path/to/SetB.bed path/to/genome.fa.fai`
    - Path to SetA.bed (str)
    - Path to SetB.bed (str)
    - Path to genome.fa.fai (str)
    - Number of permutations to perform (int, default=10,000)

3. Your program should output the following line to *stdout*:  
  `Number of overlapping bases observed: <ANSWER (int)>, p value: <ANSWER (float)>`  
  You MUST follow this print format, replacing the sections marked with `<...>` with your results. Here is an example output with numbers I made up:  
  `Number of overlapping bases observed: 345279, p value: 0.5432`

4. Your program must run in the following conda environment:  
  `conda create -c conda-forge -n bme205-assn3 "python==3.12.4" "numpy==2.1.1"`
    - This is the same as Assignment 2. You can reuse that environment if you would like.
    - *DO NOT INSTALL ANY ADDITIONAL PACKAGES. YOU CAN ONLY USE PYTHON STDLIB AND NUMPY*

## FAQ

1. What are .bed files?  
    BED files store genomic ranges, in this case it is three columns with the chromosome, start position, and stop position respectively.

2. What are .fa.fai files?  
    FA.FAI files are fasta index files which store the chromosome name, the length of the sequence, the offset from the the start of the sequence, the number of bases per line, and the number of bytes per line.

3. How do I permute genomic ranges fairly?  
    To randomly permute a set of ranges, keep the number of ranges and their lengths, and randomly reassign their starting positions. Keep in mind as you permute the ranges, should two ranges *within* a set be overlapping? How can that be avoided?

4. How many permututations should my program perform?  
    We ask for about 10,000 random permutations. **IF YOUR CODE TAKES MORE THAN 5 MINUTES TO RUN IT WILL BE CONSIDERED FAILING!**
