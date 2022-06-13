import msprime
import matplotlib.pyplot as plt
import numpy as np
import sympy

d = 5
a = 2/(d-1)
N = 1000/d
num_replicates = 1000
#_______________________________________________________


Lst= [[i,j] for i in np.arange(0,d) for j in np.arange(0,d)]
def give_coal_timeGMM_original(m, d, N, a): 
    def M(p): 
        M=np.zeros((p,p))
        for [i,j] in Lst:
            if i==j:
                M[i,j]=(1-m)
            elif i==j-1 or i==j+1:
                M[i,j]=a*0.5*m
            elif [i,j]==[0,d-1] or [i,j]==[d-1,0]:
                M[i,j]=a*0.5*m
            else:
                M[i,j]=(1-a)*(1/(d-3))*m
        return M
    M=M(d)
    def T_h(p):
        T_h=np.zeros((p,p), dtype=dict)
        for [i,j] in Lst:
            if i>j:
                T_h[i,j]=sympy.symbols(str('T_' + str(j) + str(i)))
            else:
                T_h[i,j]=sympy.symbols(str('T_' + str(i) + str(j)))
        return T_h   
    T=T_h(d)
    def SumT(p):
        SumT=np.zeros((p,p), dtype=dict)
        for [i,j] in Lst:
            SumT[i,j]=sum(M[i,k]*M[j,l]*T[k,l] for k in np.arange(p) for l in np.arange(p)) - (sum(M[i,k]*M[j,k]*T[k,k]/(2*N) for k in np.arange(p)))
        return SumT
    SumT=SumT(d)
    def F_1(p):
        F_1=np.zeros((p,p), dtype=dict)
        for [i,j] in Lst:
            F_1[i,j]=(1+SumT[i,j])-T[i,j]
        return F_1 
    F=list(F_1(d).reshape(d*d)) 
    T_list=list(T.reshape(d*d))
    T_liste1=set(T_list)
    T_liste=list(T_liste1)
    F_res=sympy.solve(F,T_liste) #kann ich nicht kontrollieren..
    def T_t(p): 
        T_t=np.zeros((p,p))
        for i in np.arange(0,p):
            for j in np.arange(0,p):
                T_t[i,j]=F_res[T[i,j]]
        return T_t    
    T_1=T_t(d)  #bis hier her getestet
    S_1=abs(sum(T_1[i,j]/(d*d) for i in np.arange(0,d) for j in np.arange(0,d)))
    H=np.zeros(d)
    for h in np.arange(d):
        H[h]=2*d*N+((d-h)*h)/(2*m)
    S_2=np.mean(H)
    S_3=(2*d*N + (d-1)*(d-1)/(2*m*d))
    return [int(S_1), int(S_2), int(S_3)]  #GMM,STST,IM



#_______________________________________________________

def give_coal_timeIM(m, d, N, num_replicates):
    samples = dict() 
    for i in range(d):
        samples[i] = 1 
    values_per_tree = (2*d)*(2*d-1)/2 
    tmrca = np.zeros(num_replicates) 
    demographyIM = msprime.Demography.island_model([N]*d, migration_rate=m/(d-1)) 
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
#_______________________________________________________

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


M_emp = [10**i for i in np.linspace(-10, -0, 10)]


means_GMM_original = list()
for m in M_emp:
    output = give_coal_timeGMM_original(m=m, d=d, N=N, a=a)  
    means_GMM_original.append(output)

 
means_IM = list()
for m in M_emp:
    output = give_coal_timeIM(m=m, d=d, N=N, num_replicates=num_replicates) 
    means_IM.append(output)

means_STST = list()
for m in M_emp:
    output = give_coal_timeSTST(m=m, d=d, N=N, num_replicates=num_replicates) 
    means_STST.append(output)    
   
means_m = list()
for m in M_emp:
    output = 2*d*N 
    means_m.append(output) 

plt.plot(M_emp, means_GMM_original, 'ks', label='GMM Formula')
plt.plot(M_emp, means_IM, 'b-', label='IM Code')
plt.plot(M_emp, means_STST, 'c-', label='STST Code')


plt.legend()
plt.xlabel("Migration rate");
plt.ylabel("Mean coalescence time");
plt.xscale('log');
plt.yscale('log');
plt.savefig('C1')
