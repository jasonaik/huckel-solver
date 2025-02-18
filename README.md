# General Hückel Solver

A command-line tool written in Python3 for calculating and printing the Hückel energies and the degeneracies of the π system of a molecule. It is assumed that all atoms are $sp^2$ hybridized and adjacent atoms are the same distance from each other. Calculations can be performed for the following molecules:

1. Linear polyene with n carbons
2. Cyclic polyene with n carbons
3. $sp^2$ hybridized Platonic solids
4. Buckminsterfullerene

## Dependencies

- numpy
- networkx
- argparse
- tabulate

## Usage

### Using positional arguments

python huckel-solver.py c 6

### Using flags

python huckel-solver.py -t cyclic-polyene -n 6

python huckel-solver.py --type cyclic-polyene --num_sites 6

## Examples

$ python huckel-solver.py c 6

| Energies | Degeneracy |
|----------|-----------|
|  2.000   |     1     |
|  1.000   |     2     |
| -1.000   |     2     |
| -2.000   |     1     |

Number of orbitals: 6
Number of unique orbitals: 4
