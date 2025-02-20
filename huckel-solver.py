import numpy as np
import networkx as nx
import argparse
from tabulate import tabulate

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
    matrix = np.zeros((n, n), dtype=int)
    
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
    matrix = np.zeros((n, n), dtype=int)

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
    matrix = np.zeros((n, n), dtype=int)
    
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
    the NetworkX package is used to generate a hardcoded mathematical graph which consists of
    nodes (representing atoms) and edges (representing bonds). This graph is then
    converted into an adjacency matrix, which represents connectivity between
    nodes (i.e. bonds). This is essentially the Hückel matrix. The function returns
    the negative adjacency matrix as is convention for Hamiltonian matrices. The 
    sparse matrix produced is also converted to a dense matrix (numpy.ndarray) 
    so that it works with NumPy operations.
    """
    if num_vertices == 4:
        graph = nx.tetrahedral_graph()
    
    elif num_vertices == 6:
        graph = nx.octahedral_graph()
        
    elif num_vertices == 8:
        graph = nx.cubical_graph()  
    
    elif num_vertices == 12:
        graph = nx.icosahedral_graph() 

    elif num_vertices == 20:
        graph = nx.dodecahedral_graph()     

    else:
        raise PlatonicSolidError(f"This is not a platonic solid")
    
    huck_matrix = nx.adjacency_matrix(graph)
    return - huck_matrix.todense()

def make_buckyball_matrix(n):
    """
    Since there is no straightfoward way to generate the Huckel matrix for buckyballs algorithmically, this function simply returns
    the hardcoded matrix for buckyballs. The argument n is just a placeholder for the CLI. Not sure if empty arguments are bad practice.
    It does save me an if-else statement I suppose.
    """
    
    # Adjacency matrix taken from https://houseofgraphs.org/graphs/1389
    adj_matrix_str = """
    011100000000000000000000000000000000000000000000000000000000
    100000000000000000000000000000000000000000000000000000110000
    100000000000000000000000000000000000000000000000000000000101
    100000000000000000000000000000000000000000000000000000001010
    000000010000110000000000000000000000000000000000000000000000
    000000010001001000000000000000000000000000000000000000000000
    000000010010000100000000000000000000000000000000000000000000
    000011100000000000000000000000000000000000000000000000000000
    000000000000101000000000100000000000000000000000000000000000
    000000000000010100000000010000000000000000000000000000000000
    000000100001000000000001000000000000000000000000000000000000
    000001000010000000000010000000000000000000000000000000000000
    000010001000000001000000000000000000000000000000000000000000
    000010000100000010000000000000000000000000000000000000000000
    000001001000000000100000000000000000000000000000000000000000
    000000100100000000010000000000000000000000000000000000000000
    000000000000010001000000000001000000000000000000000000000000
    000000000000100010000000000010000000000000000000000000000000
    000000000000001000000000000100010000000000000000000000000000
    000000000000000100000000001000100000000000000000000000000000
    000000000000000000000000010001000000000000000100000000000000
    000000000000000000000000100010000000000000001000000000000000
    000000000001000000000000000100000000010000000000000000000000
    000000000010000000000000001000000000100000000000000000000000
    000000001000000000000100000000010000000000000000000000000000
    000000000100000000001000000000100000000000000000000000000000
    000000000000000000010001000000000010000000000000000000000000
    000000000000000000100010000000000001000000000000000000000000
    000000000000000001000100000000000000000000100000000000000000
    000000000000000010001000000000000000000000010000000000000000
    000000000000000000010000010000000000000001000000000000000000
    000000000000000000100000100000000000000010000000000000000000
    000000000000000000000000000000000010000001000000000010000000
    000000000000000000000000000000000001000010000000000001000000
    000000000000000000000000001000001000000100000000000000000000
    000000000000000000000000000100000100001000000000000000000000
    000000000000000000000001000000000000010100000000000000000000
    000000000000000000000010000000000000101000000000000000000000
    000000000000000000000000000000000001010000000000000100000000
    000000000000000000000000000000000010100000000000001000000000
    000000000000000000000000000000010100000000000000010000000000
    000000000000000000000000000000101000000000000000100000000000
    000000000000000000000000000010000000000000010001000000000000
    000000000000000000000000000001000000000000100010000000000000
    000000000000000000000100000000000000000000000001010000000000
    000000000000000000001000000000000000000000000010100000000000
    000000000000000000000000000000000000000000010100000000000001
    000000000000000000000000000000000000000000101000000000000010
    000000000000000000000000000000000000000001000100000000000100
    000000000000000000000000000000000000000010001000000000001000
    000000000000000000000000000000000000000100000000000100010000
    000000000000000000000000000000000000001000000000001000100000
    000000000000000000000000000000001000000000000000000000010100
    000000000000000000000000000000000100000000000000000000101000
    010000000000000000000000000000000000000000000000000101000000
    010000000000000000000000000000000000000000000000001010000000
    000100000000000000000000000000000000000000000000010001000000
    001000000000000000000000000000000000000000000000100010000000
    000100000000000000000000000000000000000000000001000000000001
    001000000000000000000000000000000000000000000010000000000010
    """

    # .strip() removes whitespace before and after string
    # .split("\n") separates strings into lists by \n
    # list comprehension loops over all lists generated 
    # map(int, line) converts each element in list into an integer while list() makes a new list which contains integer
    adj_matrix = - np.array([list(map(int, line.strip())) for line in adj_matrix_str.strip().split("\n")])

    return adj_matrix

# Command Line Interface
if __name__ == "__main__":
    molecule_types = {
        "linear-polyene": make_linear_poly_matrix,
        "cyclic-polyene": make_cyc_poly_matrix,
        "platonic-solid": make_plat_matrix,
        "buckminsterfullerene": make_buckyball_matrix,
        "bucky": make_buckyball_matrix,
        "l": make_linear_poly_matrix,
        "c": make_cyc_poly_matrix,
        "p": make_plat_matrix,
        "b": make_buckyball_matrix,
        }
    
    # argument parser
    parser = argparse.ArgumentParser(description="Enter type and number of sites")
    
    # positional arguments
    parser.add_argument("type_pos", nargs="?", type=str, help="Type of Molecule")  
    parser.add_argument("num_sites_pos", nargs="?", type=int, help="Number of Sites")

    # keyword arguments also available
    parser.add_argument("-t", "--type", type=str, help="Type of Molecule (alt)")
    parser.add_argument("-n", "--num_sites", type=int, help="Number of Sites (alt)")

    # parse arguments
    args = parser.parse_args()
    
    molecule_type = args.type if args.type else args.type_pos
    num_sites = args.num_sites if args.num_sites else args.num_sites_pos
    
    molecule_type_lower = molecule_type.lower()
    
    if molecule_type_lower in molecule_types:
        new_mat = molecule_types[molecule_type_lower](num_sites)
        new_eigen = get_evals(new_mat)
        new_degen = find_degeneracies(new_eigen)
        
        rounded_energies = [[round(n, 3) for n in i] for i in new_degen]
        # Flatten nested list into regular list
        flattened_energies = [group[0] for group in rounded_energies]
        num_degen = [len(i) for i in new_degen]
        
        total_orbitals = sum(len(sublist) for sublist in rounded_energies)
        
        headers = ["Energies", "Degeneracy"]
        
        data = [[format(flattened_energies[n], ".3f"), num_degen[n]] for n in range(len(flattened_energies))]
        data.reverse()
        
        # print(rounded_energies)
        # print(flattened_energies)
        # print(num_degen)
        # print(data)
        
        print(tabulate(data, headers, tablefmt="grid", floatfmt=".3f"))
        
        print(f"Number of orbitals: {total_orbitals}")
        print(f"Number of unique orbitals: {len(flattened_energies)}")
    
    else:
        print(f"Error: '{molecule_type}' is not a recognized molecule type.")

        # List valid molecule types
        valid_types = ", ".join(molecule_types.keys())
        print("Valid molecule types:")
        print(valid_types)

        print("\nExample Usage:")
        print("  huckel-solver.py bucky 60")
        print("  huckel-solver.py -t cyclic-polyene -n 6")

        # Exit program with error code
        import sys
        sys.exit(1)  
