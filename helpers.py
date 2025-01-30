'''
  Contains helper functions for main.py module
'''


def get_nodes(N):
    ''' 
    Returns:
    - an array of element nodes 
    '''
    return [i * 2 / N for i in range(N + 1)]


def compute_x(interval, e):
    ''' 
    Returns:
    - the value of x(e) function for given e
    '''
    a, b = interval
    return ((a + b) / 2) + ((b - a) / 2) * e


def compute_x_derivative(interval):
    ''' 
    Returns:
    - the value of x(e) function's derivative for given e
    '''
    a, b = interval
    return (b - a) / 2


def compute_shape_functions(e):
    ''' 
    Returns:
    - an array of two shape functions Ni(e) and Nj(e) needed for gaussian quadrature
    '''
    return [0.5 * (1 - e), 0.5 * (1 + e)]


def compute_shape_functions_derivatives():
    ''' 
    Returns:
    - an array of shape functions' derivatives [Ni'(e), Nj'(e)]
    '''
    return [-0.5, 0.5]


def assemble_elements_with_global_matrices(B, L, B_e, L_e, el):
    ''' 
    Adds element matrix B_e and vector L_e to global matrices B and L respectively 
    
    Parameters:
      - B_e: element matrix
      - L_e: element vector
      - B: global matrix
      - L: global vector
      - el: element number
    '''
    for i in range(2):
        L[el + i] += L_e[i]
        for j in range(2):
            B[el + i][el + j] += B_e[i][j]


def set_boundary_conditions(B, L):
    ''' 
    Sets Dirichlet's and Cauchy's boundary conditions on B and L according to the problem statement
    
    Parameters:
      - B: left-hand side matrix of the given differential equation in a weak formulation 
      - L: right-hand side vector
    '''
    N = len(B)

    B[0][0] += 1
    L[0] += 20

    B[-1] = [0] * (N)
    B[-1][-1] = 1
    L[-1] = 0


def new_matrix(rows, cols):
    '''
    Returns:
    - new matrix of size = (rows x cols) with initial values = 0
    '''
    return [[0 for _ in range(cols)] for _ in range(rows)]


def new_vector(length):
    ''' 
    Returns:
      - new vector of size = length with initial values = 0
    '''
    return [0 for _ in range(length)]


def solve_gaussian_elimination(A, B):
    ''' 
    Solves a system of linear equations A * u = B using Gaussian elimination with 
    partial pivoting.

    Parameters:
      - A: The coefficient matrix (n x n) representing the system of equations.
      - B: The right-hand side vector (length n) representing the constants in the equations.
    
    Returns:
      - u: The solution vector containing the values of the unknowns.
    '''

    n = len(A)
    augmented_matrix = [A[i] + [B[i]] for i in range(n)]
      
    for i in range(n):
        if augmented_matrix[i][i] == 0:
            raise ValueError(f"Matrix A is singular at row {i}!")
        divisor = augmented_matrix[i][i]
        for j in range(i, n + 1):
            augmented_matrix[i][j] /= divisor

        for j in range(i + 1, n):
            factor = augmented_matrix[j][i]
            for k in range(i, n + 1):
                augmented_matrix[j][k] -= factor * augmented_matrix[i][k]

    u = [0] * n
    for i in range(n - 1, -1, -1):
        u[i] = augmented_matrix[i][n]
        for j in range(i + 1, n):
            u[i] -= augmented_matrix[i][j] * u[j]

    return u