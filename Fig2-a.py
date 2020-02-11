# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 10:24:16 2019

@author: Jin Xu
"""

import matplotlib.pyplot as plt
import numpy as np

S_tot=10e3 #unit change from microM (10e-6) to nanoM (10e-9)
E_tot=2.8e3
F_tot=2.8e3

K_E_i=[1.43e1,1.00e2,1.43e1,1.00e2]
K_F_j=[1.00e2,8.33e0,1.00e2,5.00e1] #j=i+1
lambda_i=[6.38e-2,5.05e1,1.01e-2,3.67e1]

grid=1200 #grid*grid mesh, try120,1200; 1200 is good
phi_1_data_x=[]
phi_1_data_y=[]
phi_2_data_x=[]
phi_2_data_y=[]
phi_cross_data_x=[]
phi_cross_data_y=[]
stable=[]
E1=0.
F1=0.
E2=0.
F2=0.
E3=0.
F3=0.
count1=0
count2=0
count3=0

step=12/grid #log_10(E) and log_10(F) in the range of (-6,6)
#step is also the accuracy of the mesh
equal=step/2 # assume within half of the accuracy as equal

for x in range (grid):
    E=10**(-6+x*step) #log_10(E) (-6,6)
    for y in range (grid):
        F=10**(-6+y*step) #log_10(F) (-6,6)
        t=E/F
    
        ri=[1]
        phi1=1 #initial with only r0(t)=1
        for a in range (len(lambda_i)):
            ri_add=lambda_i[a]*t**(a+1)
            for i in range (a):
                ri_add=ri_add*lambda_i[i]

            ri.append(ri_add)
            phi1=phi1+ri_add
            
        phi2=ri[0]/K_E_i[0]
        for b in range (1,len(K_E_i)):
            phi2=phi2+ri[b]/K_E_i[b]
        
        phi3=ri[1]/K_F_j[0]
        for c in range (1,len(K_F_j)):
            phi3=phi3+ri[c+1]/K_F_j[c]     
            
        S0=S_tot/(phi1+E*phi2+F*phi3)         
        phi_1_EF=E*(1.+S0*phi2)
        phi_2_EF=F*(1.+S0*phi3)
     
        if abs(np.log10(phi_1_EF)-np.log10(E_tot))<equal:
            phi_1_data_x.append(np.log10(E))
            phi_1_data_y.append(np.log10(F))
            
        if abs(np.log10(phi_2_EF)-np.log10(F_tot))<equal:
            phi_2_data_x.append(np.log10(E))
            phi_2_data_y.append(np.log10(F))
        
# prepare for Fig2(a) illustration

        if abs(np.log10(phi_1_EF)-np.log10(phi_2_EF))<equal and abs(np.log10(phi_1_EF)-np.log10(E_tot))<equal and abs(np.log10(phi_2_EF)-np.log10(F_tot))<equal:
            phi_cross_data_x.append(np.log10(E))
            phi_cross_data_y.append(np.log10(F))
            
            if (np.log10(F)>3) and (np.log10(F)<4):
                E1=E1+E
                F1=F1+F
                count1=count1+1
            else:
                if np.log10(E)>3:
                    E3=E3+E
                    F3=F3+F
                    count3=count3+1
                if (1.2<np.log10(E)) and (np.log10(E)<1.8):
                    E2=E2+E
                    F2=F2+F
                    count2=count2+1
                   
    
plt.plot(phi_1_data_x,phi_1_data_y,'g^')
plt.plot(phi_2_data_x,phi_2_data_y,'bs')
plt.plot(phi_cross_data_x,phi_cross_data_y,'ro')
plt.xlabel('log$_{10}$[E(nM)]')
plt.ylabel('log$_{10}$[F(nM)]')
plt.axis([0,5,0,5])
plt.savefig('Fig2a-1.pdf')
plt.show()


#fig2(a) illustration
E1=E1/count1
F1=F1/count1
stable=[[E1,F1]]
E2=E2/count2
F2=F2/count2
stable.append([E2,F2])
E3=E3/count3
F3=F3/count3
stable.append([E3,F3])


Si=[]
for j in range (len(stable)):
    E_stat=stable[j][0]
    F_stat=stable[j][1]
    t=E_stat/F_stat
    
    ri=[1]
    phi1=1 #initial with only r0(t)=1
    for a in range (len(lambda_i)):
        ri_add=lambda_i[a]*t**(a+1)
        for i in range (a):
            ri_add=ri_add*lambda_i[i]

        ri.append(ri_add)
        phi1=phi1+ri_add

            
    phi2=ri[0]/K_E_i[0]
    for b in range (1,len(K_E_i)):
        phi2=phi2+ri[b]/K_E_i[b]
        
    phi3=ri[1]/K_F_j[0]
    for c in range (1,len(K_F_j)):
        phi3=phi3+ri[c+1]/K_F_j[c]                
            
    S0=S_tot/(phi1+E_stat*phi2+F_stat*phi3)
        
    for k in range (len(ri)):
        Si.append(S0*ri[k]*0.001) #0.001: nM->microM

m=np.arange(15)
plt.bar(m,height=Si)
plt.xticks(m,['0','1','2','3','4','0','1','2','3','4','0','1','2','3','4'])
plt.xlabel('stable point #')
plt.ylabel('Si($\mu$M)')
plt.savefig('Fig2a-2.pdf')
plt.show()
