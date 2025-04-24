# Overview
This analysis was performed to determine if there was concordance between our ancestries-specific GWAS summary statistics. This repository was created by AS on 4/24/2023, with code written by DS. For any questions, please email: alyssa.c.scartozzi@vanderbilt.edu. 

# Concordance Analysis
Concordance analysis compares summary statistics between our genetic ancestries-specific GWAS to determine whether the concordance rate between two summary statistics is greater than expected. The concordance rate was calculated as the proportion of contiguous p-value < 0.005 variants that shared the same direction of effect (defined as a minimum distance of 10kb), divided by the total number of genomic regions. 

To establish the expected rate of concordance of effects under the null, 25 permutations were performed with random assignments of effect direction for each genomic region. Each simulation randomly shuffles the direction of effects of entire blocks (blocks defined by variants that are within 10kb and have the same direction of effect) for both genetic ancestry groups and calculates concordance based on the shuffled gwas summary statistics. The rate of concordance of effects between permuted datasets was used as the expected null concordance rate for the binomial t-test.

Please see Shaw et al. 2021 for example use of this method: 10.1016/j.ajhg.2021.11.004

# Input File Formatting for the Scripts to Work Without Editing:
Files must be tab separated with the following columns

File 1: 
MARKERNAME	CHR	POS_b37	EA	NEA	EAF	EFFECT	STDERR	PVAL

File 2:
MARKERNAME	CHR	POS_b37	EA	NEA	EAF	EFFECT	STDERR	PVAL
