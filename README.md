Coalescent-Theory.

The goal of my project was to study the continuum between the island and the one-dimensional, circular stepping stone model concerning the expected mean coalescence time. To do this, I used the general migration model (which was derived in (Nagylaki, 1998, see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1460216/pdf/9649546.pdf) and adapted the migration rates. I did this by inducing a parameter a, which corresponds to the percentage of local migration in the model. When a is set to 0, then the migrants would -  in the case of the island model - only go to d-3 out of d demes, because then two arbitrary demes and the focal deme are excluded from migration. In the cas eof the stepping stone model, the migrants would then only be going to the d-3 non-neighborung demes (as well excluding the focal deme).

The file "Formula" is a representation of this adapted formula together with the formulas for the expected mean coalescence time of the island model and the stepping stone model (which both can be found in (Slatkin, 1991, see: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/FCC418CBC6F021B741C83FDE6A0E7558/S0016672300029827a.pdf/inbreeding_coefficients_and_coalescence_times.pdf).


"Code_msprime" works as follows: 
First of all, it is good to know that the goal was to compare the performance of the formulas for the expected mean coalescence time of the island model, the stepping stone model from before with the results for the expected mean coalescence time of the island model and the one-dimensional, circular stepping stone model according to the population genetics simulator msprime. Therefore, I included the code for the formulas from the file "Formula" and added a code with msprime, which works as follows: A tuple named tmrca is created which will later contain the time to the mrca of each replicate. The demography of the respective model is created and then msprime is called to simulate num_replicate such coalescent trees with the respective demography. In order to do mimic randomness, it needs a random seed, which I, w.l.o.g., just set to 1. For each of these replicates, each pair of lineages is taken, first within each deme (to get te expected within-deme coalescence time) and then between demes (to get the expected between-deme cosalescence time). Then it returns the averaged mean of all these coalescence times that were obtained. 


_____________________________________________________________
How can the code be used?
