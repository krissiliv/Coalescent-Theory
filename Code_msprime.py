import msprime
import matplotlib.pyplot as plt
import numpy as np
import sympy

d = 5 #number of demes
a = 2/(d-1) #percentage of local migration
N = 1000/d #number of individuals per deme
num_replicates = 1000 #number of replicates
#_______________________________________________________


Lst= [[i,j] for i in np.arange(0,d) for j in np.arange(0,d)] #this is a list of all possible deme pairs
def give_coal_timeGMM_original(m, d, N, a): #then the function for the coalescent time for the general migration model is defined
    def Migration(p):#this is the migration matrix
        Migration=np.zeros((p,p)) #first we need a matrix (which is "empty" in the sense that it just contains zeros)
        for [i,j] in Lst: #then for each deme pair
            if i==j: 
                Migration[i,j]=(1-m) #we set the rate of non-miration 
            elif i==j-1 or i==j+1:
                Migration[i,j]=a*0.5*m #and the rates of migrating to one of the other demes
            elif [i,j]==[0,d-1] or [i,j]==[d-1,0]:
                Migration[i,j]=a*0.5*m #this rate should also be used if we look at the two edge demes (as I want to built a circular model => I want to connect the edges)
            else:
                Migration[i,j]=(1-a)*(1/(d-3))*m #and this is the rate of global migration (to the non-neighboring demes)
        return Migration
    Migration=Migration(d) #m should be d-dimensional
    def TimeVariable_h(p): #now I want the code to automatically define the set of variables depending on the number of demes d (which is called p here)
        TimeVariable_h=np.zeros((p,p), dtype=dict)
        for [i,j] in Lst:
            if i>j: #I want to make sure that only entries like T_(12) exist as they are equal to T_(21) entries (that holds because the model is symmetric)
                TimeVariable_h[i,j]=sympy.symbols(str('T_' + str(j) + str(i)))
            else:
                TimeVariable_h[i,j]=sympy.symbols(str('T_' + str(i) + str(j)))
        return TimeVariable_h   
    TimeVariable=TimeVariable_h(d) #this is our set of variables. In the ecological sense, the entry T_(12) of this matrix is the coalescence time of two randomly chosen lineages, where one is in deme 1 and the other in deme 2
    def SumTimes(p): #now I define sum, which connects all these coalescent times, like in (Nagylaki, 1998): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1460216/pdf/9649546.pdf
        SumTimes=np.zeros((p,p), dtype=dict)
        for [i,j] in Lst:
            SumTimes[i,j]=sum(Migration[i,k]*Migration[j,l]*TimeVariable[k,l] for k in np.arange(p) for l in np.arange(p)) - (sum(Migration[i,k]*Migration[j,k]*TimeVariable[k,k]/(2*N) for k in np.arange(p)))
        return SumTimes
    SumTimes=SumTimes(d) #this matrix contains the explicit sums for all deme pairs
    def HelpingFunction_1(p):
        HelpingFunction_1=np.zeros((p,p), dtype=dict)
        for [i,j] in Lst:
            HelpingFunction_1[i,j]=(1+SumTimes[i,j])-TimeVariable[i,j] #now I set the above defined sums to 0, so that I can later use sympy.solve to solve them after T_(ij) for all deme pairs (ij)
        return HelpingFunction_1 
    HelpingFunction=list(HelpingFunction_1(d).reshape(d*d)) 
    TimeVariables_list=list(TimeVariable.reshape(d*d))
    TimeVariables_liste1=set(TimeVariables_list)
    TimeVariables_liste=list(TimeVariables_liste1) #with the last few lines I brought TimeVariables_list and HelpFunction into the right shape so that sympy.solve can work with them
    HelpingFunction_result=sympy.solve(HelpingFunction,TimeVariables_liste) #now I can let the variables and the sums run through the "solver", which gives me the appropriate value for each T_(ij), such that all the sums hold true (it solves the equation system)
    def TimesSolution_t(p): 
        TimesSolution_t=np.zeros((p,p))
        for i in np.arange(0,p):
            for j in np.arange(0,p):
                TimesSolution_t[i,j]=HelpingFunction_result[TimeVariable[i,j]] #here I bring the solution of the equation system in the right shape suhc that I can work with them (I put them all into a matrix)
        return TimesSolution_t    
    TimesSolution_1=TimesSolution_t(d)  #Now we have the d-dim matrix of the solutions (it contains each T_(ij) value)
    Sum_generalModel=abs(sum(TimesSolution_1[i,j]/(d*d) for i in np.arange(0,d) for j in np.arange(0,d))) #as my goal is to compare the expected mean coalescence times with each other, I calculate the mean of the coalescence times T_(ij)
    H=np.zeros(d)
    for h in np.arange(d):
        H[h]=2*d*N+((d-h)*h)/(2*m) #this formula and the formula in S_3 can be found in (Slatkin, 1991) https://www.cambridge.org/core/services/aop-cambridge-core/content/view/FCC418CBC6F021B741C83FDE6A0E7558/S0016672300029827a.pdf/inbreeding_coefficients_and_coalescence_times.pdf
    Sum_STSTmodel=np.mean(H) #This has nothing to do with the above, it is a different, less complex formula (for the expected mean coalescence time) of a very restricted model: the Stepping Stone model
    Sum_Imodel=(2*d*N + (d-1)*(d-1)/(2*m*d)) #This also has nothing to do with the above, it is a different, less complex formula (for the expected mean coalescence time) of a very restricted model: the Island model
    return [int(Sum_generalModel), int(Sum_STSTmodel), int(Sum_Imodel)]  #GMM,STST,IM



#_______________________________________________________

def give_coal_timeIM(m, d, N, num_replicates):
    samples = dict() 
    for i in range(d):
        samples[i] = 1 # 1 sample (lineage) is taken out of every deme
    values_per_tree = (2*d)*(2*d-1)/2
    tmrca = np.zeros(num_replicates) #this will later contain the time until the MRCA (most recent common ancestor) of the lineage pairs
    demographyIM = msprime.Demography.island_model([N]*d, migration_rate=m/(d-1)) #now I set up the island model, with symmetric migration rates. Migration goes to each other deme at equal rate
    replicates = msprime.sim_ancestry(samples=samples, num_replicates=num_replicates, demography=demographyIM, random_seed=1) #then I define what a replicate is (it can be seen as an ancestral tree of the model)
    for replicate_index, ts in enumerate(replicates):
        tree = ts.first() #tree should be the ancestral tree within that replicate
        acc_within = 0.  #first set the within-deme coalescence time as 0 (coalescence time of a pair of lineages, where both are in the same deme)
        acc_between = 0.  #as well as the between-deme coalescence time (coalescence time of a pair of lineages, where both are in different demes)
        for i in range(2*d): 
            for j in range(i+1, 2*d): 
                if i%2==0 and j==i+1: #the lineages are numbered, such that for example 4 (4 can be divided by 2 without a rest) and 5 (=i+1) are in the same deme
                    acc_within += tree.time(tree.mrca(i, j)) #therefore I calculate the within deme coalescence time
                else:
                    acc_between += tree.time(tree.mrca(i, j)) # else they are in two different demes and we need the between-deme coalescence time
        tmrca[replicate_index] = (acc_within+acc_between) / values_per_tree #of all the solutions (all the coalescence times) for this replicate we now want to know the mean of them
    return(np.mean(tmrca)) #after we got the mean of each replicate, we calculate the mean over them so that we get a good representor of the expected mean coalescence time of the model
#_______________________________________________________
#this code is analogeous to the above but with the Stepping Stone model
def give_coal_timeSTST(m, d, N, num_replicates):
    samples = dict() 
    for i in range(d):
        samples[i] = 1 
    values_per_tree = (2*d)*(2*d-1)/2 
    tmrca = np.zeros(num_replicates)
    demographyIM = msprime.Demography.stepping_stone_model([N]*d, migration_rate=m/2, boundaries=False)
    replicates = msprime.sim_ancestry(samples=samples, num_replicates=num_replicates, demography=demographyIM, random_seed=1)
    for replicate_index, ts in enumerate(replicates):
        tree = ts.first()  
        acc_within = 0.  
        acc_between = 0.  
        for i in range(2*d):
            for j in range(i+1, 2*d): 
                if i%2==0 and j==i+1:
                    acc_within += tree.time(tree.mrca(i, j))
                else:
                    acc_between += tree.time(tree.mrca(i, j))
        tmrca[replicate_index] = (acc_within+acc_between) / values_per_tree
    return(np.mean(tmrca)) 

#__________plotting_____________________________________________


M_emp = [10**i for i in np.linspace(-10, -0, 10)] #this is the range of the parameter m, which is the migration rate of the model


means_GMM_original = list()
for m in M_emp:
    output = give_coal_timeGMM_original(m=m, d=d, N=N, a=a)[0]  #I evaluate the above defined function with the above defined parameters d,N and a for m in the range defined with M_emp
    means_GMM_original.append(output) #and I append all the outputs to the list (which will later be the "y-axis" in the plot)

#Analogeously I get the values for the "vertical axis" for the Island model and the Stepping Stone model function defined above:
means_IM = list()
for m in M_emp:
    output = give_coal_timeIM(m=m, d=d, N=N, num_replicates=num_replicates) 
    means_IM.append(output)

means_STST = list()
for m in M_emp:
    output = give_coal_timeSTST(m=m, d=d, N=N, num_replicates=num_replicates) 
    means_STST.append(output)    

#and this formula here represents the value of the expected mean coalescence time for the simplest, most restricted biomathematical model: the Wright-Fisher model:
means_m = list()
for m in M_emp: #I iterate over m because I want to have exactly as many values as in the above "vertical axis' lists" in order to be able to compare it with them in the plot
    output = 2*d*N 
    means_m.append(output) 

#nor I tell the plotter what the "horizontal axis" and "vertical axis" should be as well as what the color and shape of the plotted dots/line should be and the name it should have in the plot legend
plt.plot(M_emp, means_GMM_original, 'ks', label='GMM Formula')
plt.plot(M_emp, means_IM, 'b-', label='IM Code')
plt.plot(M_emp, means_STST, 'c-', label='STST Code')


plt.legend()
plt.xlabel("Migration rate"); #name of the horizontal axis
plt.ylabel("Mean coalescence time"); #name of the vertical axis
plt.xscale('log'); #scaling of the axis
plt.yscale('log');
plt.savefig('C1') #name, with which the plot will be saved (in the same folder as the code-file)
