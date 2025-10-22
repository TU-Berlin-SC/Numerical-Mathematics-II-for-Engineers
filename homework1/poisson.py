"""
* Group 9
* Griselda Lleshi, Hyejun Park, Jacqueline Mc Kenna Reardon
Poisson.py
- Runs tests for const and sin problems, N=[20,40,80,...,1280]
"""

# ------------------------
# Begin Poisson_template
# ------------------------
"""
Numerical Mathematics for Engineers II WS 25/26
Homework 01 Exercise 1.1
1D Poisson equation 
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy
import time
import os
import time # for the t_solution error resolution
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# Define the test problem:
# only -u''(x) = f(x)

def get_testproblem(testfunction):
    testproblem = {}

    # xL,xR : domain [xL,xR]
    testproblem["xL"] = 0.0
    testproblem["xR"] = 1.0

    # choice
    testproblem["choice"] = testfunction
    
    # right-hand-side function and exact solution
    if (testfunction == "const"):
        testproblem["f"] = lambda x: x * 0 + 1 # f(x) = 1
        testproblem["uexact"] = lambda x: 0.5 * x * (1 - x)  # u(x) = 1/2x(1 - x)
    elif (testfunction == "sin"):
        testproblem["f"] = lambda x: np.sin(np.pi * x) # f(x) = sin pi x
        testproblem["uexact"] = lambda x : np.sin(np.pi * x) / (np.pi ** 2)  # u(x) = sin pi x / pi^2
    else:
        raise Exception(
            'Stop in testproblem. Choice of test problem does not exist')

    return testproblem



# Set parameters for solving the problem
def define_default_parameters():
    parameters = {}

    # nrefine: how many refinements do we do?
    parameters["Nrefine"] = 4

    # N: number of inner grid points (on coarsest grid)
    parameters["N"] = 20

    # plot_freq: how often do we plot?
    parameters["plot_freq"] = 0

    return parameters


# plot computed and exact solution
def graph(U_comp, x_plot, uexact, xL, xR):
    # evaluate true solution
    U_true = uexact(x_plot)

    plt.figure(1)
    plt.ion()

    plt.figure(1)
    plt.plot(x_plot, U_comp, 'r.', markersize=4, label='computed solution')
    plt.plot(x_plot, U_true, 'k-', label='true solution')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.xlim(xL, xR)
    plt.legend()
    save_str = 'results/solution.eps'
    plt.savefig(save_str)

    plt.show()
    plt.pause(2.0)
    plt.clf()


def get_rhs_diag(f,N,x,dx):
    # create right-hand-side vector
    # x[0] = xL, x[N+1] = xR, 내부점 x[1:N]
    rhs = dx**2 * f(x[1:-1])
    
    # compute diagonals of matrix and store each of them in a vector,
    # which we will use later on to set up the matrix
    l = -np.ones(N-1)
    d = 2*np.ones(N)
    r = -np.ones(N-1)
    
    return rhs, l, d, r

def my_driver(solver_type, testproblem, parameters, N):
    # Extract problem information
    xL = testproblem["xL"]
    xR = testproblem["xR"]
    f = testproblem["f"]
    uexact = testproblem["uexact"]
    plot_freq = parameters["plot_freq"]

    # Grid generation
    # calculate dx: cell length in space
    dx = (xR - xL) / (N + 1)

    # create mesh with boundary nodes
    x = np.linspace(xL, xR, N + 2)

    solution_left = uexact(xL)
    solution_right = uexact(xR)

    rhs, l, d, r = get_rhs_diag(f,N,x,dx)

    # solve with full matrix
    if solver_type == "full":
       

        # create full matrix
        matrix = np.zeros((N, N)) # init with 0 NxN array
        # zip : pairs elements from two lists together
        # [l,d,r] : left, center, right diagonal values
        # [-1,0,1] : offsets(position) for diagonals -> 0: main diagonal, -1: lower diagonal, 1: upper diagonal
        # np.diag(etries,offset) : add entries vectors to offset diagonals of matrix
        for entries, offset in zip([l, d, r], [-1, 0, 1]):
            matrix += np.diag(entries, offset) # add each diagonal to the matrix -> tridiagonal matrix

        # start time rec.
        t0 = time.time()

        # note: the result of the last three lines can also be accomplished
        # matrix = functools.reduce(
        #     lambda a, b: a + np.diag(b[0], b[1]), zip([l, d, r], [-1, 0, 1]), np.zeros((N, N)))

        # solve
        solution = np.linalg.solve(matrix, rhs)

        # end time rec.
        t_solution = time.time() - t0
       

    elif solver_type == "sparse":
        # create sparse tridiagonal matrix
        offsets = [-1, 0, 1]
        data = [l, d, r]
        sparse_matrix = diags(data, offsets, format='csr')
        
        # start timer
        t0 = time.time()
        
        # solve sparse system
        solution = spsolve(sparse_matrix, rhs)
        
        # end timer
        t_solution = time.time() - t0


    # crate full solution (including boundary nodes)
    full_sol = np.zeros(N + 2)
    full_sol[0] = solution_left
    full_sol[1:-1] = solution
    full_sol[-1] = solution_right

    # plot
    if plot_freq != 0:
        graph(full_sol, x, uexact, xL, xR)

    # compute error
    true_sol = uexact(x)
    err = full_sol - true_sol
    err_max = np.max(np.abs(err))

    print('Error in Max norm:\t %3.2e\n' % err_max)
    print('Solved in:\t\t %3.2e seconds\n' % t_solution)

    return err_max, t_solution

# ------------------------
# End Poisson_template
# ------------------------


# ------------------------
# Main: Run all tests
# ------------------------
if __name__ == "__main__":
    if not os.path.exists("results"):
        os.makedirs("results")

    parameters = define_default_parameters()

    testproblems = ["const", "sin"]
    N_values = [20, 40, 80, 160, 320, 640, 1280]

    # 사용자에게 solver 선택
    print("Choose solver type ('full' or 'sparse'):")
    solver_type = input().strip()

    for problem_name in testproblems:
        print(f"\nSolving test problem: {problem_name}")
        testproblem = get_testproblem(problem_name)
        err_max_list = []

        for N in N_values:
            err_max, _ = my_driver(solver_type, testproblem, parameters, N)
            err_max_list.append(err_max)

        print("\nN\tMax error")
        for N, err in zip(N_values, err_max_list):
            print(f"{N}\t{err:.2e}")

        # 수렴율 계산
        err_max_array = np.array(err_max_list)
        rate = np.log(err_max_array[:-1]/err_max_array[1:])/np.log(2)
        print("Convergence rates:", rate)
