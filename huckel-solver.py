import numpy as np
import networkx

def get_evals(matrix):
    """
    Returns eigenvalues of a matrix
    """
    eigenvalues = np.linalg.eigvals(matrix)
    return np.sort(eigenvalues)

def make_linear_poly_matrix(n):
    """
    Creates a linear polyene huckel matrix
    """
    
    if n < 1:
        raise ValueError("n must be > 1 for a linear poly-ene")
    
    # Create n x n matrix with zeros
    matrix = np.zeros((n, n))
    
    # fill in matrix
    for i in range(n):
        if i > 0:
            matrix[i, i-1] = -1 # Resonance integral beta
            matrix[i-1, i] = -1 # symmetric beta on the other side
    
    return matrix

def check_poly_eigen(n, eigenvalues, print_all=True, beta=-1):
    """
    Checks eigenvalues against general solution for linear system
    """
    lambdas = []
    # Get eigenvalues for all values of n
    for k in range (1, n + 1):
        # General solution
        eig = 2*beta*np.cos(k * np.pi/ (n+1))
        lambdas.append(eig)
        
    sorted_lambdas = np.sort(lambdas)
    
    diff = eigenvalues - sorted_lambdas
    
    # return diff, with some tolerance for rounding differences
    sig_diffs = [i for i in diff if round(i, 5) > 0]
    
    if print_all:
        return [diff, sig_diffs]
    else:
        return [sig_diffs]
        
def make_and_eval_poly(n):
    """
    2 in 1 function to make huckel matrix and get eigenvalues for convenience
    """
    # Create n x n matrix with zeros
    matrix = np.zeros((n, n))
    
    # fill in matrix
    for i in range(n):
        if i > 0:
            matrix[i, i-1] = -1 # Resonance integral beta
            matrix[i-1, i] = -1 # symmetric beta on the other side
            
    eigenvalues = np.linalg.eigvals(matrix)
    return np.sort(eigenvalues)


def make_cyc_poly_matrix(n):
    """
    Creates a cyclic polyene huckel matrix
    """
    
    if n < 3:
        raise ValueError("At least n = 3 required to make cyclic poly-ene")
    
    # Create n x n matrix with zeros
    matrix = np.zeros((n, n))
    
    # Fill in the matrix
    for i in range(n):
        if i >0:
            matrix[i, i-1] = -1
            matrix[i-1, i] = -1
    
    # For cyclic polyene, there will be beta at positions (n, 1) and (1, n)
    matrix[n-1, 0] = -1
    matrix[0, n-1] = -1
    
    return matrix

def find_degeneracies(eigen, tolerance=1e-8):
    """
    Takes eigenvalues and determines degeneracies by first sorting them and then 
    checking for similar values at a given tolerance to account for differences
    caused by rounding. 
    """
    
    # sort eigenvalues from smallest to biggest
    sorted_eig = np.sort(eigen)
    degens = []
    current_degen = [sorted_eig[0]]
    
    # Loop through sorted eigenvalues
    for eig in sorted_eig[1:]:
        # If difference between eig and smallest eig in current group, add to current group
        if abs(eig - current_degen[0]) < tolerance:
            current_degen.append(eig)
            
        # If difference is larger, append current group to overall group and reassign to eig
        else:
            degens.append(current_degen)
            current_degen = [eig]
        
    degens.append(current_degen)
    
    return degens
    
    
class PlatonicSolidError(Exception):
    pass

def make_plat_matrix(num_vertices):
    """
    Take the number of vertices, if this number corresponds to a platonic solid,
    the NetworkX package is used to generate a mathematical graph which consists of
    nodes (representing atoms) and edges (representing bonds). This graph is then
    converted into an adjacency matrix, which represents connectivity between
    nodes (i.e. bonds). This is essentially the Hückel matrix. The function returns
    the negative adjacency matrix as is convention for Hamiltonian matrices. The 
    sparse matrix produced is also converted to a dense matrix so that it works
    with NumPy operations.
    """
    if num_vertices == 4:
        graph = networkx.tetrahedral_graph()
    
    elif num_vertices == 6:
        graph = networkx.octahedral_graph()
        
    elif num_vertices == 8:
        graph = networkx.cubical_graph()  
    
    elif num_vertices == 12:
        graph = networkx.icosahedral_graph() 

    elif num_vertices == 20:
        graph = networkx.dodecahedral_graph()     

    else:
        raise PlatonicSolidError(f"This is not a platonic solid")
    
    huck_matrix = networkx.adjacency_matrix(graph)
    return - huck_matrix.todense()

# def buckminster():
    


# print(make_poly_matrix(4))

# print(get_evals(make_poly_matrix(4)))

# print(check_poly_eigen(100, get_evals(make_poly_matrix(100)), True))

# print(check_poly_eigen(6, make_and_eval_poly(6), print_all=True))

# print(make_cyc_poly_matrix(6))

# print(find_degeneracies(get_evals(make_cyc_poly_matrix(10))))


print(make_plat_matrix(20))
