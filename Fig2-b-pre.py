# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 10:24:16 2019

@author: Jin Xu
"""

import numpy as np
#from decimal import Decimal #overflow

S_tot=10e3 #unit change from microM (10e-6) to nanoM (10e-9)
E_tot=2.8e3
F_tot=2.8e3


K_E_i=[1.43e1,1.00e2,1.43e1,1.00e2]
K_F_j=[1.00e2,8.33e0,1.00e2,5.00e1] #j=i+1
lambda_i=[6.38e-2,5.05e1,1.01e-2,3.67e1]


#grid=1200 #grid*grid mesh, try12,120,1200; 1200 is good
grid=30
phi_1_data_xy=[]
phi_2_data_xy=[]


step=5/grid #log_10(E) and log_10(F) in the range of (-6,6)
#step is also the accuracy of the mesh
equal=step/10 # assume within 1/10 of the accuracy as equal

# step one have the values on the curve of E_tot and F_tot
for x in range (1,grid):
    E=10**(0+x*step) 
    for y in range (grid):
        F=10**(0+y*step) 
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
            #phi_1_data_x.append(np.log10(E))
            #phi_1_data_y.append(np.log10(F))
            phi_1_data_xy.append([np.log10(E),np.log10(F)])
            
        if abs(np.log10(phi_2_EF)-np.log10(F_tot))<equal:
            #phi_2_data_x.append(np.log10(E))
            #phi_2_data_y.append(np.log10(F))
            phi_2_data_xy.append([np.log10(E),np.log10(F)])
        
#print(len(phi_1_data_xy))
#print(len(phi_2_data_xy))
print(phi_1_data_xy)
#print(phi_2_data_xy)


for j in range (len(phi_1_data_xy)):
    Si=[]
    E=10**phi_1_data_xy[j][0]
    F=10**phi_1_data_xy[j][1]
    print(E,F)

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
        
    for k in range (len(ri)):
        Si.append(S0*ri[k]) #0.001: nM->microM
        
    print(Si)
