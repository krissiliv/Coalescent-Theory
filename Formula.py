import sympy
import numpy as np
import matplotlib.pyplot as plt


d= 5 #demes
a= 1 #a% local migration
N= 1000/d #number of individuals per deme
#_______________________________________________________



Lst= [[i,j] for i in np.arange(0,d) for j in np.arange(0,d)] #this is a list of all possible deme pairs
def give_coal_timeGMM_original(m): #then the function for the coalescent time for the general migration model is defined
    def M(p):#this is the migration matrix
        M=np.zeros((p,p)) #first we need a matrix (which is "empty" in the sense that it just contains zeros)
        for [i,j] in Lst: #then for each deme pair
            if i==j: 
                M[i,j]=(1-m) #we set the rate of non-miration 
            elif i==j-1 or i==j+1:
                M[i,j]=a*0.5*m #and the rates of migrating to one of the other demes
            elif [i,j]==[0,d-1] or [i,j]==[d-1,0]:
                M[i,j]=a*0.5*m #this rate should also be used if we look at the two edge demes (as I want to built a circular model => I want to connect the edges)
            else:
                M[i,j]=(1-a)*(1/(d-3))*m #and this is the rate of global migration (to the non-neighboring demes)
        return M
    M=M(d) #m should be d-dimensional
    def T_h(p): #now I want the code to automatically define the set of variables depending on the number of demes d (which is called p here)
        T_h=np.zeros((p,p), dtype=dict)
        for [i,j] in Lst:
            if i>j: #I want to make sure that only entries like T_(12) exist as they are equal to T_(21) entries (that holds because the model is symmetric)
                T_h[i,j]=sympy.symbols(str('T_' + str(j) + str(i)))
            else:
                T_h[i,j]=sympy.symbols(str('T_' + str(i) + str(j)))
        return T_h   
    T=T_h(d) #this is our set of variables. In the ecological sense, the entry T_(12) of this matrix is the coalescence time of two randomly chosen lineages, where one is in deme 1 and the other in deme 2
    def SumT(p): #now I define sum, which connects all these coalescent times, like in (Nagylaki, 1998): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1460216/pdf/9649546.pdf
        SumT=np.zeros((p,p), dtype=dict)
        for [i,j] in Lst:
            SumT[i,j]=sum(M[i,k]*M[j,l]*T[k,l] for k in np.arange(p) for l in np.arange(p)) - (sum(M[i,k]*M[j,k]*T[k,k]/(2*N) for k in np.arange(p)))
        return SumT
    SumT=SumT(d) #this matrix contains the explicit sums for all deme pairs
    def F_1(p):
        F_1=np.zeros((p,p), dtype=dict)
        for [i,j] in Lst:
            F_1[i,j]=(1+SumT[i,j])-T[i,j] #now I set the above defined sums to 0, so that I can later use sympy.solve to solve them after T_(ij) for all deme pairs (ij)
        return F_1 
    F=list(F_1(d).reshape(d*d)) 
    T_list=list(T.reshape(d*d))
    T_liste1=set(T_list)
    T_liste=list(T_liste1) #with the last few lines I brought T and F into the right shape so that sympy.solve can work with them
    F_res=sympy.solve(F,T_liste) #now I can let the variables and the sums run through the "solver", which gives me the appropriate value for each T_(ij), such that all the sums hold true (it solves the equation system)
    def T_t(p): 
        T_t=np.zeros((p,p))
        for i in np.arange(0,p):
            for j in np.arange(0,p):
                T_t[i,j]=F_res[T[i,j]] #here I bring the solution of the equation system in the right shape suhc that I can work with them (I put them all into a matrix)
        return T_t    
    T_1=T_t(d)  #Now we have the d-dim matrix of the solutions (it contains each T_(ij) value)
    S_1=abs(sum(T_1[i,j]/(d*d) for i in np.arange(0,d) for j in np.arange(0,d))) #as my goal is to compare the expected mean coalescence times with each other, I calculate the mean of the coalescence times T_(ij)
    H=np.zeros(d)
    for h in np.arange(d):
        H[h]=2*d*N+((d-h)*h)/(2*m) #this formula and the formula in S_3 can be found in (Slatkin, 1991) https://www.cambridge.org/core/services/aop-cambridge-core/content/view/FCC418CBC6F021B741C83FDE6A0E7558/S0016672300029827a.pdf/inbreeding_coefficients_and_coalescence_times.pdf
    S_2=np.mean(H) #This has nothing to do with the above, it is a different, less complex formula (for the expected mean coalescence time) of a very restricted model: the Stepping Stone model
    S_3=(2*d*N + (d-1)*(d-1)/(2*m*d)) #This also has nothing to do with the above, it is a different, less complex formula (for the expected mean coalescence time) of a very restricted model: the Island model
#    S_w = 2*N #this is the within deme coalescence time of the models (I put it here in case I would like to plot it as well, but it is not plotted below)
#    S_b = 2*d*N + (d-1)/(2*m) #this is the between deme coalescence time of the models (I put it here in case I would like to plot it as well, but it is not plotted below)
    return [S_1, S_2, S_3] #GMM,STST,IM  #, S_w, S_b]  can also be added - if needed

#---------plotting------------------------------------------------------------------------
M_emp = [10**i for i in np.linspace(-13, -1, 20)]

means_GMM_original = list()
means_IM_formula = list()
means_STST_formula = list()

for m in M_emp:
    output_GMM_original = give_coal_timeGMM_original(m=m)
    means_GMM_original.append(output_GMM_original[0]) #this is the first entry of line 68 (which is the output of give_coal_timeGMM_original(m))
    means_IM_formula.append(output_GMM_original[2])  #this is the third entry of the output of give_coal_timeGMM_original(m)
    means_STST_formula.append(output_GMM_original[1]) #and this is the second one
    
#-----------------------------------------------------------------------------------------
    
#here I specify the horizontal and vertical axis, the color and shape of the dots/line and the name of it
plt.plot(M_emp, means_IM_formula, 'b--', label='IM Formula')
plt.plot(M_emp, means_STST_formula, 'c--', label='STST Formula')
plt.plot(M_emp, means_GMM_original, 'kD', label='GMM Formula')

plt.legend() #I want to have a legend in the plot
plt.xlabel("Migration rate"); #this is the name of the horizontal axis
plt.ylabel("Mean coalescence time"); # this is the name of the vertical axis
plt.xscale('log'); #this is the scaling of the axis
plt.yscale('log');
plt.savefig('F1') #this is the name, with which the plot will be saved in the folder of the code-file
