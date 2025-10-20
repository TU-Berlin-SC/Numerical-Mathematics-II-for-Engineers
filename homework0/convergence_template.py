"""
Numerical Mathematics for Engineering II WS 25/26
Homework 00 Exercise 0.3
Revision on convergence curves and tables 
"""

# Import necessary libraries
import numpy as np
import scipy.integrate # library to calculate integrate
import matplotlib.pyplot as plt

# Methods for Integrating Functions given fixed samples.

#    trapezoid            -- Use trapezoidal rule to compute integral.
#    cumulative_trapezoid -- Use trapezoidal rule to cumulatively compute integral.
#    simpson              -- Use Simpson's rule to compute integral from samples.

def f(x):
    return 2 + np.cos(x)

# interval
a = 0
b = np.pi/2

def approximate_integral(n,rule):
    h = (b-a)/n # step size
    x = np.linspace(a, b, n+1) # range [a,b] => chop with n sections

    if(rule == "trapezoid"):
        integral = scipy.integrate.trapezoid(f(x), x)
    elif(rule == "simpson"):
        integral = scipy.integrate.simpson(f(x), x)

    # i got pi + 1, when i evaluated the integral
    exact = np.pi + 1
    relative_error = np.abs((integral - exact)/exact)     # E_rel = (|I_approx - I_exact|) / I_exact

    return relative_error 



if __name__ == '__main__':
    N = 5
    error_t = np.zeros(N)
    error_s = np.zeros(N)
    for k in range(N):
        error_t[k] = approximate_integral(2**(k+1), "trapezoid")
        error_s[k] = approximate_integral(2**(k+1), "simpson")

    h = (b-a) / (2**np.arange(1,N+1)) # N+1 due to dimension error
    
    plt.plot(h, error_t, "o-", label="trapezoid")  # using dot line
    plt.plot(h, error_s, "s-", label="Simpson")     # using square line
    plt.xlabel("$h$")
    plt.ylabel("relative error")
    plt.legend()
    plt.show()
    # loglog slope = p
    plt.loglog(h, error_t, "x-", label="trapezoid")
    plt.loglog(h, error_s, "x-", label="Simpson")
    plt.loglog(h, h**2, "--", label="$h^2$")   # Trapezoid rule -> O(h^2)
    plt.loglog(h, h**4, "--", label="$h^4$")   # Simpson rule -> O(h^4)

    plt.xlabel("$h$")
    plt.ylabel("relative error")
    plt.legend()
    plt.show()
    
    ### Convergence tables
    # E(h) = C * h^p  =>  log(E(h)) = log(C) + p * log(h)
    # p = [log(E_i/E_(i+1)) / (log(h_i /h_(i+1)) ]
    rate_t = np.log(error_t[:-1] / error_t[1:]) / np.log(h[:-1] / h[1:]) 
    rate_s = np.log(error_s[:-1] / error_s[1:]) / np.log(h[:-1] / h[1:])

    print("")
    print(f"   h     |   E_T    |     p_T    |   E_S    |    p_S    ")
    print(f"{h[0]:.2e} | {error_t[0]:.2e} |     --     | {error_s[0]:.2e} | -- ")
    for i in range(1,N):
        print(f"{h[i]:.2e} | {error_t[i]:.2e} | {rate_t[i-1]:.4e} | {error_s[i]:.2e} | {rate_s[i-1]:.4e}")
    print("")
