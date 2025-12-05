## Social loops in gut microbiome strain-sharing networks

This github repo describes workflow and codes used in Honduras socio-microbial cycle study. Specifically, it provides codebase for constructing loops and cross-analyzing network cycles across multiple modalities.

### Abstract:

The human microbiome is transmitted across diverse social ties.  But broader network-level motifs, such as geodesic cycles, or loops of individuals, in which people might be repeatedly exposed, remain unexplored.  Social loops of connected individuals can reinforce re-infection by microbes, eventually leading to higher probability of colonization of microbes in new hosts. In addition, microbes also form their own microbial networks, and some microbial species may also demonstrate similar cyclic patterns of propagation. Here, we examine these geodesic cycles in microbial and social networks among 1,787 rural Hondurans in 18 isolated villages. We find that cycles are an over-represented motifs in all of these microbial species (validated using null villages). In addition, a histogram of the cycle distribution per person showed a log-normal distribution, which suggests multiplicative growth. We also find significant association between microbial cycle centrality and social cycle centrality in four species: SGB8640 (Candidatus Gastranaeiphilales bacterium), SGB1667 (Prevotellaceae family), Prevotella Pectinovera (SGB1662), and Clostridia bacterium (SGB14005), showing existence of social microbes. Moreover, increasing overlap between social and microbial transmission pattern was positively associated with increased social cycle centrality in individuals, which suggests that social ties may be used as underlying paths for microbial transmission as well. 

### Code content:
# kcycle_strains_final_code.R
Code for (i) Pre-processing microbiome data and preparing networks (ii) Constructing cycles (iii) Performing association across cross-modal cycles
(iv) Measuring overlap across cross-modal cycles

# kcycle_strains_random_permute tests.R
Code for randomizing vilalges for assessing statistical viability of cycle formation beyond random chance. (Null hypothesis/robustness test)
