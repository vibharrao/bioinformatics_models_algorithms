# Assignment 2
## Description
You are a bioinformatician tasked with analyzing gene expression data. 

These data have been assembled with large-scale replication for 10 genes.

In the `control_files` directory, you will find 100 csv files, each representing an experimental unit with gene expression data for the 10 genes under control conditions.

In the `treatment_files` directory, you will find 100 csv files, each representing an experimental unit with gene expression data for the 10 genes under treatment conditions.

Assignment (Tier 1):
1. Normalize the gene expression counts within each replicate experimental measurements. To do this, you can simply divide each gene expression count by the total expression measured from a given replicate.
2. Compute the mean and median normalized expression for each gene in the experimental groups.
3. Compute the log_2 Fold Change for each gene such that a positive fold change indicates increased expression in the treatment group.
4. Output gene IDs ranked by the log_2 Fold Change (lowest first) in a *tab-delimited* table with the following header:
#gene   (mean normalized control expression)    (median normalized control expression)  (mean noramlized treatment expression)   (median normalized treatment expression) (logFoldChange)


We will run your script as follows:
`$ python your_assignment.py path/to/controls path/to/treatments`


Extra Credit (Tier 2):
[do this if you finished the assignment in less than 4 hours]

As a diligent statistician you recognize that the fold change of two genes is not comparable as they are expressed at different levels. To determine the significance of the change in expression, you will calculate the p-value of the change in the median of the treatments using a Mann-Whitney U test.
Using the experiment and treatment files as input, write a script that can:

1. normalize the gene expression counts for experiment and treatment datasets. 
2. compare the normalized gene expression estimates using a mann-whitney U test to determine if the median expression is different between experiment and treatment experiments.
3. output a differential expression matrix, sorted by lowest p-value, with the following format and including a header:
#gene   (mean normalized control expression)    (median normalized control expression)  (mean noramlized treatment expression)   (median normalized treatment expression) (logFoldChange) (p value)

We will run your script as follows:
`$ python your_assignment.py path/to/controls path/to/treatments`


Extra Credit 2 (Tier 3):
[do this if you finish all of the above in less than 4 hours]

Changes to one gene's expression does slightly affect normalized expression estimates for all other genes because it's a measure of *relative* expression. 
With most gene expression analysis, only a few genes out of thousands might be significantly differentially expressed. 
So the relative effect on normalized expression estimates for other genes is spread across thousands and relatively minimal. 
In this experiment, we only measured expression at 10 genes. 
Explain in one paragraph or less what mildly clever thing I did when creating this experimental dataset such that not all 10 genes are apparently differentially expressed when comparing their normalized expression. 

For this tier, please turn in your script for Tier 2 as well as the paragraph as a seperate file, called `paragraph.txt`

## Assignment Requirements:

### Filepaths:
Your script will be interacting with files and your TAs will not be running the code on your computer.
We will be passing the paths to the data files using the command line, as reflected above. 
We will compare the output you turn in to the output your script produces, as well as the expected output.

Submission Requirements:
1. You **must** submit your solution for Tier 1 with the following naming scheme:

        Firstname_Lastname_Tier1.py, output_Tier1.txt

2. If you attempted Tier 2 for extra credit, **you must submit all of the above, plus the following**:

        Firstname_Lastname_Tier2.py, output_Tier2.txt

3. If you attempted Tier 3, **you must submit all of the above, plus the following**:

       paragraph.txt


Hints:

        1. command line arguments are accessible via the `argv` variable in the builtin `sys` library.
        2. python's builtin `os` library has a listdir function which lists all files in a directory.
        3. python 3.10 has a builtin `statistics` library with a NormalDist function.
