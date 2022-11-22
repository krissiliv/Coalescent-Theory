import sympy
import numpy as np
import matplotlib.pyplot as plt


d= 5 #demes
a= 1 #a% local migration
N= 1000/d #number of individuals per deme
#_______________________________________________________


#for this function (the steps until line 65) please refer to https://github.com/krissiliv/Coalescent-Theory/blob/main/Code_msprime.py
Lst= [[i,j] for i in np.arange(0,d) for j in np.arange(0,d)]
def give_coal_timeGMM_original(m): 
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
    F_res=sympy.solve(F,T_liste) 
    def T_t(p): 
        T_t=np.zeros((p,p))
        for i in np.arange(0,p):
            for j in np.arange(0,p):
                T_t[i,j]=F_res[T[i,j]]
        return T_t    
    T_1=T_t(d)  
    S_1=abs(sum(T_1[i,j]/(d*d) for i in np.arange(0,d) for j in np.arange(0,d)))
    H=np.zeros(d)
    for h in np.arange(d):
        H[h]=2*d*N+((d-h)*h)/(2*m) #this formula can be found in (Slatkin, 1991)
    S_2=np.mean(H)
    S_3=(2*d*N + (d-1)*(d-1)/(2*m*d)) #this formula can be found in (Slatkin, 1991)
#    S_w = 2*N #this is the within deme coalescence time of the models (I put it here in case I would like to plot it as well, but it is not plotted below)
#    S_b = 2*d*N + (d-1)/(2*m) #this is the between deme coalescence time of the models (I put it here in case I would like to plot it as well, but it is not plotted below)
    return [S_1, S_2, S_3] #GMM,STST,IM  #, S_w, S_b]  can also be added - if needed

#---------plotting------------------------------------------------------------------------
M_emp = [10**i for i in np.linspace(-13, -1, 20)] #this is the range of the parameter m, which represents the migration rate

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
