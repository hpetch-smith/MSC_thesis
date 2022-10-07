# -*- coding: utf-8 -*-

import math;
import numpy as np;
from scipy.linalg import svdvals

pi = math.pi;



#--------------dnorm_log-----------------------------------------------
# Density function
def dnorm_log(x,mu,sigma):
    p = math.log(1.0/math.sqrt(2 * pi * sigma**2)) -0.5 * ((x-mu)/sigma)**2;
    return p;


#--------------C_ii-------------------------------------------------
# Calculates autocovariance within each model
def C_ii(par, d):
    
    A = math.exp(par[2]);
    B = math.exp(par[3]);
    v = math.exp(par[0]);
    w = math.exp(par[1]);
    
    C_ii_u = pi**0.5 * v**2 * 1.0/math.sqrt(A) * math.exp(-A*d**2*0.25);
    C_ii_v = pi**0.5 * w**2 * 1.0/math.sqrt(B) * math.exp(-B*d**2*0.25);
    
    C_ii = C_ii_u + C_ii_v;
    return C_ii;


# #--------------C_12-------------------------------------------------
# Calculates cross covariance between model 1 and 2

def C_12(par,d):
    A_1 = math.exp(par[2]);
    A_2 = math.exp(par[3]);
    v_1 = par[0];
    v_2 = par[1];
    miu = par[4];
    Sigma = A_1*(A_1 + A_2)**(-1)*A_2;
    C_12_u = (2*pi)**0.5*v_1*v_2*1.0/math.sqrt(A_1 + A_2) * math.exp(-Sigma*(d-miu)**2*0.5);
    return C_12_u


# #--------------C_21-------------------------------------------------
# Calculates cross covariance between model 2 and 1

def C_21(par, d):
    A_1 = math.exp(par[2]);
    A_2 = math.exp(par[3]);
    v_1 = par[0];
    v_2 = par[1];
    miu = par[4];
    
    Sigma = A_1*(A_1 + A_2)**(-1)*A_2;
    C_21_u = (2*pi)**0.5*v_1*v_2*1.0/math.sqrt(A_1+A_2) * math.exp(-Sigma*(d+miu)**2*0.5);
    return C_21_u;
    
#--------------Covariance-------------------------------------------------   
# Construct covariance matrices using the above covariance functions  
# parameters:; [v1 w1 f1 g1 Beta1 v2 w2 f2 g2 Beta2 mu];

def Covariance(par, t1, t2):
    # Parameters for independant 1
    par1 = par[0:5]; 
    
    # Parameters for independant 2
    par2 = par[5:10]; 
    # Shared parameters
    par3 = [par[0], par[5], par[2], par[7], par[10]];  
    # Number of samples in model one 
    N1 = len(t1);
    # Number of samples in model 2
    N2 = len(t2);
    # Total dimensions for matrices
    dim = N1 + N2; 
    
    # Set up blank covariance matrices
    I = np.ndarray(shape = (dim,dim));
    Cv = np.ndarray(shape = (dim, dim)); 
    C11 = np.ndarray(shape = (N1, N1));
    C22 = np.ndarray(shape = (N2, N2));
    C12 = np.ndarray(shape = (N1, N2));
    C21 = np.ndarray(shape = (N2, N1));
    count = 0;
    # Model one autocovariance
    for i in range(N1):
        for j in range(N1):
            C11[i,j] = C_ii(par1, t1[i] - t1[j])
            if i == j:
                I[i,j] = math.exp(par1[4])**2;

    # Model two auto covariance
    for i in range(N2):
        for j in range(N2): 
            C22[i,j] = C_ii(par2, t2[i]-t2[j])
            if i == j:
                I[(N1 + i),(N1 + j)] = math.exp(par2[4])**2;
                
    # Cross covariance for 1 and 2
    for i in range(N1):
        for j in range(N2): 
            C12[i,j] = C_12(par3, t1[i]-t2[j])
            
    for i in range(N2):
        for j in range(N1):
            C21[i,j] = C_21(par3, t2[i]-t1[j])
            
    # Construct final matrix
    Cv[0:N1,0:N1] = C11; 
    Cv[(N1):dim + 1,(N1):dim] = C22; 

    Cv[(N1):dim + 1,0:N1] = C21       
    Cv[0:N1,(N1):dim + 1] = C12      
    Cv = Cv + I
    return Cv
 
#--------------log_lik-------------------------------------------------
def log_lik(par,Y,t1,t2,prior):
    
    N = len(Y)
    C = Covariance(par, t1, t2)
    C_I = np.linalg.inv(C)

    
    prior_v1 = dnorm_log(par[0], prior[0], prior[1])
    prior_v2 = dnorm_log(par[5], prior[0], prior[1])
    
    prior_w1 = dnorm_log(par[1], prior[2], prior[3])
    prior_w2 = dnorm_log(par[6], prior[2], prior[3])
    
    prior_f1 = dnorm_log(par[2], prior[4], prior[5])
    prior_f2 = dnorm_log(par[7], prior[4], prior[5])

    prior_g1 = dnorm_log(par[3], prior[6], prior[7])
    prior_g2 = dnorm_log(par[8], prior[6], prior[7])

    prior_b1 = dnorm_log(par[4], prior[8], prior[9])
    prior_b2 = dnorm_log(par[9], prior[8], prior[9])
    
    prior_mu = dnorm_log(par[10], prior[10], prior[11])
    
    ty = np.transpose(Y)
    M = ty.dot(C).dot(C_I)

    v = prior_v1 + prior_v2
    w = prior_w1 + prior_w2
    f = prior_f1 + prior_f2
    g = prior_g1 + prior_g2
    beta = prior_b1 + prior_b2

    SVD = svdvals(C) 
    A = np.log10(SVD)
    A = np.sum(A)
    print(A)
    log_lkh = -1*(-0.5 * A - 0.5 * M - 0.5 * N * math.log(2*pi) + v + w + f + g + beta + prior_mu)
    
    return log_lkh;