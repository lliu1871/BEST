# BEST
BEST is a free phylogenetics package of programs to estimate the posterior distribution of species trees using multilocus molecular data that accounts for deep coalescence but not for other issues such as horizontal transfer. It is intended to implement the Bayesian hierarchical model proposed by Liang Liu and Dennis Pearl and further developed in collaboration with Scott Edwards. Details of the Bayesian hierarchical model are available in Liang's Thesis (pdf format).

Estimating the posterior distribution of species trees involves two consecutive Markov Chain Monte Carlo (MCMC) procedures. The first MCMC is performed in BEST part 1 (a modified version of MrBayes 3.1.2) in which a new function is added to approximate the joint probability of gene trees. The output gene trees (multilocus) then form the input file of the second MCMC program BEST part 2. Beginning with version 1.6, the weights for the importance sampling component (details in Liang's Thesis) are also calulated within BEST.



The Bayesian Hierarchical model has been applied to several mutilocus datasets including finch, yeast, and macaques (NEXUS files). Interesting results from these analyses can be found in this technical report and in publications by Liu and Pearl (Systematic Biology, 56:504-514, 2007) and by Edwards, Liu, and Pearl (PNAS, 104:5936-5941, 2007).
