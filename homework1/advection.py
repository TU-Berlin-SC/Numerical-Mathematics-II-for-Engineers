"""
advections.py
- Contains full advection_template code unchanged
- Runs tests for 1D linear advection using FTCS and Upwind schemes
- Investigates CFL effects
"""

# ------------------------
# Begin advection_template
# ------------------------
"""
Numerical Mathematics for Engineers II WS 25/26
Homework 01 Exercise 1.3
1D Linear advection
"""

import math
import matplotlib.pyplot as plt
import numpy as np
import os

# Define the test problem:
# u_t + a u_x = 0
# (d)
def get_testproblem():
    testproblem = {}

    # t_0: initial time
    testproblem["t0"] = 0.0

    # xL,xR : domain [xL,xR]
    testproblem["xL"] = 0.0
    testproblem["xR"] = 1.0

    # advection speed
    testproblem["a"] = 1.0

    # u0: initial data
    testproblem["u0"] = lambda x: np.sin(2 * math.pi * x)

    # uexact: exact solution
    testproblem["uexact"] = lambda x, t: testproblem["u0"]((x-testproblem["a"]*t) % testproblem["xR"])

    # tend: final time
    testproblem["tend"] = 1.0

    return testproblem


# Set parameters for solving the problem
def define_default_parameters():
    parameters = {}

    # nrefine: how many refinements do we do?
    parameters["Nrefine"] = 0

    # N: number of grid points (on coarsest grid)
    parameters["N"] = 40

    # max_steps: maximal number of time steps
    parameters["max_steps"] = 10000

    # CFL: CFL number used
    parameters["CFL"] = 0.8

    # plot_freq: how often do we plot?
    parameters["plot_freq"] = 1

    return parameters


def graph_solution(U, x, time, uexact, xL, xR, method, initialize, final_time):
    # remove the ghost cells at the boundaries
    x_plot = x[1:-1]
    U_comp = U[1:-1]

    # evaluate true solution
    U_true = uexact(x_plot, time)

    if initialize:
        plt.figure(1)
        plt.ion()

    plt.figure(1)
    plt.plot(x[1:-1], U_comp, 'r.', markersize=4, label='computed solution')
    plt.plot(x[1:-1], U_true, 'k-', label='true solution')
    plt.title(method)
    plt.xlabel('x')
    plt.ylabel('u')
    plt.xlim(xL, xR)

    if final_time == 1:
        plt.figure(2)
        plt.plot(x[1:-1], U_comp, 'r.', markersize=4, label='computed solution')
        plt.plot(x[1:-1], U_true, 'k-', label='true solution')
        plt.title(method)
        plt.xlabel('x')
        plt.ylabel('u')
        plt.xlim(xL, xR)
        plt.legend()
        save_str = 'results/' + method + '.eps'
        plt.savefig(save_str)

    plt.show()
    plt.pause(0.01)
    plt.clf()


#### Question (a) 
def compute_dt(CFL, a, dx):
    dt = CFL * dx / a
    return dt


#### Question (b)
#
def update_ftcs(U, a, dt, dx):
    U_temp = U.copy()
    U_temp[1:-1] = U[1:-1] - (a * dt / (2*dx)) * (U[2:] - U[:-2])
    U[:] = U_temp

def update_upwind(U, a, dt, dx):
    U_temp = U.copy()
    U_temp[1:-1] = U[1:-1] - (a * dt / dx) * (U[1:-1] - U[:-2])
    U[:] = U_temp
    

#### Question (c)
def compute_error(x, time, uexact, U):
    U_exact = uexact(x, time)
    err_max = np.max(np.abs(U - U_exact))
    return err_max
#
#######

# Driver
#
# Output:
# - err_max: error in maximum norm
def my_driver(method, testproblem, parameters, N):

    # Extract problem information
    xL = testproblem["xL"]
    xR = testproblem["xR"]
    a  = testproblem["a"]
    t0 = testproblem["t0"]
    tend = testproblem["tend"]
    u0 = testproblem["u0"]
    uexact = testproblem["uexact"]

    max_steps = parameters["max_steps"]
    CFL = parameters["CFL"]
    plot_freq = parameters["plot_freq"]

    # Grid generation
    dx = (xR - xL) / N
    x = np.linspace(xL - dx, xR, N + 2)

    # initialize U
    U = u0(x)
    time = t0
    done = 0

    # Plot initial data
    if plot_freq != 0:
        graph_solution(U, x, time, uexact, xL, xR, method, True, False)

    # Time stepping
    for j in range(1, max_steps + 1):

        U[0] = U[-2]
        U[-1] = U[1]

        dt = compute_dt(CFL, a, dx)

        if (time + dt) > tend:
            dt = tend - time
            done = 1
        time = time + dt

        if (plot_freq != 0) and (j % plot_freq) == 0:
            print('Taking time step %i: \t update from %f \t to %f' % (j, time - dt, time))

        if method == 'FTCS':
            update_ftcs(U, a, dt, dx)
        elif method == 'upwind':
            update_upwind(U, a, dt, dx)
        else:
            raise ValueError('Stop in my_driver. No appropriate method chosen.')

        if (plot_freq != 0) and (j % plot_freq) == 0:
            graph_solution(U, x, time, uexact, xL, xR, method, False, False)

        if done == 1:
            print('Have reached time tend; stop now')
            break

    if j >= max_steps:
        print('Stopped after %i steps.' % max_steps)
        print('Did not suffice to reach the end time %f.' % tend)

    if plot_freq != 0:
        graph_solution(U, x, time, uexact, xL, xR, method, False, True)

    U[-1] = U[1]
    err_max = compute_error(x[1:], time, uexact, U[1:])
    print('Error in maximum norm:\t %3.2e\n' % err_max)

    return err_max

# ------------------------
# End advection_template
# ------------------------


# ------------------------
# Main: Run all tests
# ------------------------
def main():
    if not os.path.exists("results"):
        os.makedirs("results")

    testproblem = get_testproblem()

    # 사용자에게 solver 선택
    print("Choose scheme ('FTCS' or 'upwind'):")
    scheme = input().strip()

    parameters = define_default_parameters()

    # --- Part D: FTCS and Upwind comparison ---
    print("=== FTCS Scheme ===")
    parameters["Nrefine"] = 2  # N = 40, 80, 160
    my_driver("FTCS", testproblem, parameters, parameters["N"])

    print("=== Upwind Scheme ===")
    my_driver("upwind", testproblem, parameters, parameters["N"])

    # --- Part E: CFL influence for Upwind ---
    CFL_values = [0.8, 1.0, 1.2]
    N = 40
    for CFL in CFL_values:
        print(f"\n--- Upwind Scheme with CFL = {CFL} ---")
        parameters = define_default_parameters()
        parameters["CFL"] = CFL
        parameters["Nrefine"] = 0
        my_driver("upwind", testproblem, parameters, N)


if __name__ == "__main__":
    main()
