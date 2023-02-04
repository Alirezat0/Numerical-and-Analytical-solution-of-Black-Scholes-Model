#!/usr/bin/env python
# coding: utf-8

# In[21]:


import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

def C_BS(S, K, T, r, q, sigma):
    d1 = (np.log(S/K) + (r - q + sigma**2/2)*T)/(sigma*np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    C = S*np.exp(-q*T)*stats.norm.cdf(d1) - K*np.exp(-r*T)*stats.norm.cdf(d2)
    return C

def dC(C, S, t, r, q, sigma, K):
    dS = S*sigma*np.sqrt(t)
    C_new = C_BS(S + dS, K, T-t, r, q, sigma)
    
    return (C_new - C)/dS
# Parameter
K = 100
r = 0.12
q = 0
sigma = 0.1
T = 1

# Numerical Solution
dt = T/100
S = np.linspace(70, 130, 1001)
C_num = np.zeros(len(S))
C_num[0] = C_BS(S[0], K, T, r, q, sigma)
for i in range(1, len(S)):
    C_new = C_num[i-1] + dC(C_num[i-1], S[i-1], T-i*dt, r, q, sigma, K)*dt
    C_num[i] = np.array([C_new])
C_num = C_num[::-1]


# Analytic Solution
C_BS_vec = np.vectorize(C_BS)
C_ana = C_BS_vec(S, K, T, r, q, sigma)

# 2D Plot
plt.plot(S, C_ana, 'r', label='Analytic Solution')
#plt.plot(S, C_num, 'g--', label='Numerical Solution')
plt.legend()
plt.xlabel('Stock Price')
plt.ylabel('Option Price')
plt.title('Black-Scholes Option Pricing')
plt.show()




