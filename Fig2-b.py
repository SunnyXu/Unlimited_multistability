# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 11:03:38 2019

@author: Jin Xu
"""

import tellurium as te
import roadrunner
import random
import numpy as np
import pylab as plt

# to load from SBML
# r = roadrunner.RoadRuner(r'C:\path\to\sbml.xml')

r = te.loada("""
             
#Reactions:
    E+S0 -> ES0;  a_E_0*E*S0-b_E_0*ES0;
    #S0 -> ES0;    a_E_0*E*S0-b_E_0*ES0;
    ES0  -> E+S1; c_E_0_1*ES0;
    E+S1 -> ES1;  a_E_1*E*S1-b_E_1*ES1;
    ES1  -> E+S2; c_E_1_2*ES1;
    E+S2 -> ES2;  a_E_2*E*S2-b_E_2*ES2;
    ES2  -> E+S3; c_E_2_3*ES2;
    E+S3 -> ES3;  a_E_3*E*S3-b_E_3*ES3;
    #ES3 -> S4;   c_E_3_4*ES3;
    ES3 ->  E+S4; c_E_3_4*ES3;
   
    S4+F -> FS4;  a_F_4*S4*F-b_F_4*FS4;
    FS4  -> S3+F; c_F_4_3*FS4;
    S3+F -> FS3;  a_F_3*S3*F-b_F_3*FS3;
    FS3  -> S2+F; c_F_3_2*FS3;
    S2+F -> FS2;  a_F_2*S2*F-b_F_2*FS2;
    FS2  -> S1+F; c_F_2_1*FS2;
    S1+F -> FS1;  a_F_1*S1*F-b_F_1*FS1;
    #FS1 -> S0;    c_F_1_0*FS1;
    FS1  -> S0+F; c_F_1_0*FS1;
         
#Rules: keep update as part of the reactions 
#    ES0:=E*S0/K_E_0; ES1=E*S1/K_E_1; ES2=E*S2/K_E_2; ES3=E*S3/K_E_3
#    FS1=F*S1/K_F_1; FS2=F*S2/K_F_2; FS3=F*S3/K_F_3; FS4=F*S4/K_F_4

#Species initializations:     
    # X=E+F, X=0, so E=F=0

    # select the point on the curve of E_tot (F_tot), refer to paper supplementary 2.2
    
    
    #E = 5.6 # stable point 1, curled 100,1000, E:20
    #F = 1000
    
    #E = 5.6 #0.74
    #F = 562 #2.75
    
    
    # three stable points E:30, culled 100,1000; phi_2
    E = 14.68 
    F = 14.68
    
    
    #E=14.68 #1.17 #phi_1
    #F=21.54 #1.33 


    S1 = 1; #S0+S4=S_tot, so S1=S2=S3=0
    S2 = 2;
    S3 = 3;

 
    S0 = 2e3; #set initial alpha=0.2
    S4 = 8e3;

#Parameters: only give the value once as am initial condition    
    a_E_0=8.12e-3; a_E_1=1.02e-1;  a_E_2=8.12e-3;  a_E_3=1.02e-1
    b_E_0=1.6e-2;  b_E_1=2.04e-1;  b_E_2=1.6e-2;   b_E_3=2.04e-1
    c_E_0_1=1e-1;  c_E_1_2=1e1;    c_E_2_3=1e-1;   c_E_3_4=1e1
    a_F_1=1.12e-1; a_F_2=2.64e-3;  a_F_3=6.51e-1;  a_F_4=2.85e-3
    b_F_1=2.24e-1; b_F_2=5e-3;     b_F_3=1.3;      b_F_4=6e-3
    c_F_1_0=1.1e1; c_F_2_1=1.7e-2; c_F_3_2=6.39e1; c_F_4_3=1.36e-1
    
#    K_E_0=1.43e1;  K_E_1=1.00e2;   K_E_2=1.43e1;   K_E_3=1.00e2
#    K_F_1=1.00e2;  K_F_2=8.33e0;   K_F_3=1.00e2;   K_F_4=5.00e1
    

""")



S_tot=10e3 #unit change from microM (10e-6) to nanoM (10e-9)
#Si=r.getFloatingSpeciesIds()
#print(Si)

for x in range (100): #100 samples will be tested
    alpha=random.uniform(0,1)
    r.model["init(S0)"]=alpha*S_tot
    r.model["init(S4)"]=(1-alpha)*S_tot
    #r.time=np.log10(r.time)
    #r.S4=np.log10(r.S4)
    m1 = r.simulate (0.1, 100,  20000, ['time','S4'])
    m2 = r.simulate (100, 1000, 20000, ['time','S4'])
    m3 = r.simulate (1000, 1000000, 20000, ['time','S4'])
    
    m = np.vstack((m1,m2,m3))
    #r.plot(m, logx=True, xlim=[10E-2,10E3], ylim=[10E-9,10E3]);

#    plt.plot(m[:,0],m[:,1:])
    plt.plot(np.log10(m[:,0]),np.log10(m[:,1:]))
    #plt.plot(np.log10(m[:,0]),m[:,1:])
    
# =============================================================================
# for x in range (100): #100 samples will be tested
#     alpha=random.uniform(0,1)
#     r.model["init(S0)"]=alpha*S_tot
#     r.model["init(S4)"]=(1-alpha)*S_tot
#     r.model["init(E)"]=5.6
#     r.model["init(F)"]=562
#     #r.time=np.log10(r.time)
#     #r.S4=np.log10(r.S4)
#     m1 = r.simulate (0.1, 100,  20000, ['time','S4'])
#     m2 = r.simulate (100, 1000, 20000, ['time','S4'])
#     m3 = r.simulate (1000, 1000000, 20000, ['time','S4'])
#     
#     m = np.vstack((m1,m2,m3))
#     #r.plot(m, logx=True, xlim=[10E-2,10E3], ylim=[10E-9,10E3]);
# 
#     #plt.plot(m[:,0],m[:,1:])
#     #plt.plot(np.log10(m[:,0]),np.log10(m[:,1:]))
#     #plt.plot(np.log10(m[:,0]),m[:,1:])
# =============================================================================

    
# =============================================================================
#     r.conservedMoietyAnalysis = True
# 
#     #This step is for steady analysis,  
#     #to check if all are real stable steady states
#     #obtain the Jacobian matrix
#     Jac=r.getFullJacobian()    
#     #print(Jac)
#     #print(Jac.shape)
# 
#     #print the eigenvalues of the full Jacobian matrix
#     Eigen=r.getFullEigenValues()
#     print(Eigen.real)
#     print(len(Eigen))
#     for i in range (len(Eigen)):
#         if Eigen.real[i] >0:
#             print("hi")
# =============================================================================
    
plt.xlabel('log$_{10}$[Time(s)]')
plt.ylabel('log$_{10}$[S4(nM)]')
plt.savefig('S4.pdf')
plt.show()



# =============================================================================
# #This part is for stability analysis
# #Evaluation of the steady state
# r.steadyStateSelections=['S4']
# val=r.getSteadyStateValues()
# print(val)
# 
# #a value representing the sum of square of the rates of change.
# #10e-6 represents a steady value
# 
# import random
# 
# S_tot=10e3 #unit change from microM (10e-6) to nanoM (10e-9)
# 
# a_E_i=[8.12e-3,1.02e-1,8.12e-3,1.02e-1]
# b_E_i=[1.6e-2,2.04e-1,1.6e-2,2.04e-1]
# c_E_i_j=[1e-1,1e1,1e-1,1e1]
# a_F_j=[1.12e-1,2.64e-3,6.51e-1,2.85e-3]
# b_F_j=[2.24e-1,5e-3,1.3,6e-3]
# c_F_j_i=[1.1e1,1.7e-2,6.39e1,1.36e-1]
# 
# K_E_i=[1.43e1,1.00e2,1.43e1,1.00e2]
# K_F_j=[1.00e2,8.33e0,1.00e2,5.00e1] #j=i+1
# 
# #initial condition
# E=F=0 #X=0, so E=F=0
# S1=S2=S3=0 #S0+S4=S_tot, so S1=S2=S3=0
# =============================================================================


