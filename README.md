# BEST
BEST is a free phylogenetics package of programs to estimate the posterior distribution of species trees using multilocus molecular data that accounts for deep coalescence but not for other issues such as horizontal transfer. It is intended to implement the Bayesian hierarchical model proposed by Liang Liu and Dennis Pearl and further developed in collaboration with Scott Edwards.

Estimating the posterior distribution of species trees involves two consecutive Markov Chain Monte Carlo (MCMC) procedures. The first MCMC is performed in BEST part 1 (a modified version of MrBayes 3.1.2) in which a new function is added to approximate the joint probability of gene trees. The output gene trees (multilocus) then form the input file of the second MCMC program BEST part 2. Beginning with version 1.6, the weights for the importance sampling component (details in Liang's Thesis) are also calulated within BEST.


Citation:
1. Liu, L. BEST: Bayesian estimation of species trees under the coalescent model. Bioinformatics 2008, 24(21):2542:2543.
2. Liu, L. and D.K. Pearl. Species trees from gene trees: reconstructing Bayesian posterior distributions of a species phylogeny using estimated gene tree distributions. Syst Biol 2007, 56:504-514.
