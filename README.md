# Overview
This analysis was performed to determine if there was concordance between our ancestries-specific GWAS summary statistics. This repository was created by AS on 4/24/2023, with code written by DS. For any questions, please email: alyssa.c.scartozzi@vanderbilt.edu. 

# Concordance Analysis
Concordance analysis compares summary statistics between our genetic ancestries-specific GWAS to determine whether the concordance rate between two summary statistics is greater than expected. The concordance rate is calculated by the proportion of variants that had the same direction of effect over the total variants present in both GWAS. p-value = .05 was used as our threshold. 

Since we were comparing summary statistics within two different genetic ancestries, we wanted to take into account LD. To do this, the script simulates concordance analyses X times and creates an expected concordance based on the average concordance of those simulations. Each simulation randomly shuffles the direction of effects of entire blocks (blocks defined by variants that are within 10kb and have the same direction of effect) for both genetic ancestry groups and calculates concordance based on the shuffled gwas summary statistics.

Please see Shaw et al. 2021 for example use of this method: 10.1016/j.ajhg.2021.11.004
