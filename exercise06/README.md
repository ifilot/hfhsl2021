# Hartree-Fock course for the Han-sur-Lesse Winterschool of 2021::Exercise 06

[![This work is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported License](https://i.creativecommons.org/l/by-nc/3.0/88x31.png)](http://creativecommons.org/licenses/by-nc/3.0/)

## Problem description
With a working implementation for H2, it is logical that we are going to expand
towards non-trivial molecules. Let us stick with two atoms to make matters
not overly complicated and simulate the CO molecule. In this exercise, you are
tasked to modify our earlier implementation to accomplish this goal. Furthermore,
we are going to focus a bit more on the energies of the molecular orbitals
and the linear combination of atomic orbitals that constitute them.

In the past exercises, we build the CGFs corresponding to the atomic orbitals
by explicitly specifying the coefficients. Such practices will become quickly
very impractical. For example, for the CO molecule, we have to define ten sets
of two coefficients to set-up the basis set. As such, we are going to make
use of a custom routine that does this process automatically for us. There is
no magic behind this procedure: the values are stored in a JSON file and the
method simply extracts these from the JSON file and builds the basis set for
us. The code is shown below:

```python
# build molecule
mol = Molecule('CO')
mol.add_atom('C', 0.0, 0.0, -2.116/2.0)
mol.add_atom('O', 0.0, 0.0,  2.116/2.0)

cgfs, nuclei = mol.build_basis('sto3g')
nelec = np.sum([n[1] for n in nuclei])
N = len(cgfs)
```

## Questions
Take [the result of the previous exercise](../exercise05/ex05_solution.py) 
and modify the script to simulate
CO rather than H2. Take a critical look at the for loops and modify the
number of CGFs they need to iterate over. The result of your calculation
should produce an electronic energy of `-111.223473` HT.

Once you have this working, extract the eigenvalues and -vectors from
the simulation and take a critical look at these. Try to answer the
following questions:

1. How many molecular orbitals are there in total? How does this correspond
   to the number of basis functions?
2. How many of the molecular orbitals are occupied?
3. Which of the orbitals are non-bonding? 
4. Which of the orbitals are degenerate?
5. Which 1s orbital lies lower in energy, the one of C or the one of O?

## Solution
The solution is given in [ex06_solution.py](ex06_solution.py)

The answers to the questions are as follows:
1. Given a basis set of 10 AOs (basis functions), all integral matrices
   will be 10x10 matrices and thus we can expect to obtain 10 solutions.
   In general, the number of solutions always matches the number of
   basis functions.
2. The CO molecule has in total 14 electrons (8 for O, 6 for C) and
   thus the number of occupied molecular orbitals is 14/2 = 7.
3. Non-bonding orbitals can be identified by checking for solution
   vectors where there is only one coefficient with a value very
   close to 1.0 while all other values are very close to 0.0.
   For example, the first and second molecular orbitals correspond
   to this pattern and it can be said that these molecular orbitals
   are synonymous to their corresponding atomic orbitals.
4. Degenerate orbitals can be identified by them having the same
   energy. Molecular orbitals (5,6) and (8,9) are two pairs
   of double-degenerate orbitals. These are known as the 1pi
   bonding and 2pi* anti-bonding orbitals.
5. The first molecular orbital corresponds to the 1s atomic orbital
   of oxygen and thus lies lower in energy than the 1s atomic orbital
   of carbon. This result is to be expected as O has a higher
   nuclear charge than C and thus the electron experiences a stronger
   interaction for O than for C.  