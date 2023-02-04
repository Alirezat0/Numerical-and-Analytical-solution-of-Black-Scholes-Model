# Numerical-and-Analytical-solution-of-Black-Scholes-Model
there is 2 code : first one is for analytic solution 2d (also it contains a 3d plot). second one is numerical solution in 2d where x is stock price($) and y is option price it calculates 2 mentioned solution where: K = 100 strike price r = 0.12 risk-free interest rate q = 0 dividend yield sigma = 0.1 volatility of the stock price T = 1 time to maturity

This function C_BS calculates the Black-Scholes formula for the price of a European call option

d1 and d2 are intermediate variables used in the Black-Scholes formula. C is the final option price, which is calculated as the sum of two terms representing the present value of the expected stock price and the present value of the strike price

This function dC calculates the derivative of the option price with respect to the stock price

dS is an intermediate variable representing a small change in the stock price. C_new is the option price calculated using the C_BS function with the stock price increased by dS. The final output of the function is the derivative of the option price with respect to the stock price, calculated as the difference between C_new and C divided by dS.
