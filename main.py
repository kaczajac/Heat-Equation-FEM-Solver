import matplotlib.pyplot as plt
from helpers import *
from math import sqrt


### Gauss quadrature parameters


gauss_points = [ (-1 / sqrt(3)), (1 / sqrt(3)) ]
gauss_weight = 1


### Predefined functions k(x) and f(x)


def k(x):
    if 0 <= x <= 1: return 1
    else: return 2*x

def f(x):
    return 100 * x


### Plotting function


def plot_result(nodes, u, N):
    plt.figure(num = "Heat Equation | FEM Solution")
    plt.title(f"1D Heat Equation Solution For {N} Elements")
    plt.plot(nodes, u)
    plt.xlabel("x - position")
    plt.ylabel("u(x) - temperature")
    plt.grid()
    plt.show()


### Main computations


def compute_elements(interval):
    B_e = new_matrix(2, 2)
    L_e = new_vector(2)

    for e in gauss_points:
      x = compute_x(interval, e)
      dx_de = compute_x_derivative(interval)
      n = compute_shape_functions(e)
      dn_de = compute_shape_functions_derivatives()

      for i in range(2):
          L_e[i] += -f(x) * n[i] * dx_de * gauss_weight
          for j in range(2):
              B_e[i][j] += -k(x) * dn_de[i] * dn_de[j] * (1 / dx_de) * gauss_weight

    return (B_e, L_e)


def get_matrices(nodes, N):
    B = new_matrix(N + 1, N + 1)
    L = new_vector(N + 1)

    for element in range(N):
        integration_interval = ( nodes[element], nodes[element + 1] )
        B_e, L_e = compute_elements(integration_interval)
        assemble_elements_with_global_matrices(B, L, B_e, L_e, element)
      
    set_boundary_conditions(B, L)
    return (B, L)


def solve_equation(N):
    nodes = get_nodes(N)
    B, L = get_matrices(nodes, N)
    u = solve_gaussian_elimination(B, L)
    plot_result(nodes, u, N)


if __name__ == "__main__":
    N = int(input("Enter the number of elements N: "))
    solve_equation(N)