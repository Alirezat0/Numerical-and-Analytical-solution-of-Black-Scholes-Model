#!/usr/bin/env python
# coding: utf-8

# In[21]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm

def European_Option_Price(S, K, T, r, sigma, q, option_type, div_schedule):
    """
    Returns the price of a European option using the Black-Scholes model.
    S: initial stock price
    K: strike price
    T: time to expiration (in years)
    r: risk-free interest rate
    sigma: volatility
    q: dividend rate
    option_type: 'call' or 'put'
    div_schedule: a list of tuples, where each tuple is the date and dividend amount on that date
    """
    # Constants
    d1 = (np.log(S / K) + (r - q + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    
    # Option Price
    if option_type == 'call':
        price = S * np.exp(-q * T) * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    elif option_type == 'put':
        price = K * np.exp(-r * T) * norm.cdf(-d2) - S * np.exp(-q * T) * norm.cdf(-d1)
    
    return price
# Parameters
S = 100 # initial stock price
K = 100 # strike price
r = 0.03 # risk-free interest rate
q = 0.03 # dividend rate
sigma = 0.3 # volatility
Tav = 1 # time to expiration (in years)
option_type = 'call' # option type
div_schedule = [] # dividend schedule

# Stock Price
S_min = 50
S_max = 150
n = 100
S = np.linspace(S_min, S_max, n)

# Time to Expiration
tao_min = 0.5
tao_max = 2
n = 100
Tav = np.linspace(tao_min, tao_max, n)
T, S = np.meshgrid(Tav, S)

# Option Price
C_total = European_Option_Price(S, K, T, r, sigma, q, option_type, div_schedule)

# 3D Plot
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(121, projection='3d')
ax.plot_surface(T, S, C_total, linestyle=':')
ax.set_xlabel('Time Until Expiration (Years)')
ax.set_ylabel('Stock Price ($)')
ax.set_zlabel('Option Price ($)')

# 2D Plot
ax = fig.add_subplot(122)
ax.plot(S, C_total[:, 0], color='blue', linestyle=':')
ax.set_xlabel('Stock Price ($)')
ax.set_ylabel('Option Price ($)')



plt.show()



# In[ ]:





# In[ ]:




