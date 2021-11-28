# Hartree-Fock course for the Han-sur-Lesse Winterschool of 2021::Exercise 03

[![This work is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported License](https://i.creativecommons.org/l/by-nc/3.0/88x31.png)](http://creativecommons.org/licenses/by-nc/3.0/)

## Problem description
In this exercise, we are going to evaluate the integrals for the H2 molecule.
Because the H2 molecule has two CGFs (one for each H atom), we are going
to obtain matrices rather than scalar values for the integrals. Instead of a single value 
for the overlap and kinetic integrals, we will obtain 2x2 matrices. For the 
nuclear attraction integral, we in fact obtain two such 2x2 matrices (one for each H atom).
Finally, there will be 16 two-electron integrals which will contain a total
of 4 unique values.

Since the CGF coefficients for both of the H atoms are identical, we can readily
set-up the second one from the first one. Note that this has to be done via
the deepcopy feature of Python, which works as follows

```python
# construct the STO-3g CGF for H
cgf1 = cgf([0.0, 0.0, 0.0])
cgf1.add_gto(0.15432897000000001, 6.3624213899999997, 0, 0, 0)
cgf1.add_gto(0.53532813999999995, 1.1589229999999999, 0, 0, 0)
cgf1.add_gto(0.31364978999999998, 0.44463454000000002, 0, 0, 0)
pos1 = (0,0,0) # position of first H nucleus

# create the second CGF for a H atom at 1.4 a.u. distance from the
# first one; the H-H bonding axis is aligned over the z-axis
pos2 = (0,0,1.4)        # position of second H nucleus
cgf2 = deepcopy(cgf1)
cgf2.p = pos2           # reset its center
```

To store the matrices, we are going to make use of the `numpy` package
of Python. To build a (2x2) matrix of only zeros, one can use the 
following code.

```python
# build 2x2 placeholders
S = np.zeros((2,2))     # overlap matrix
```

To evaluate the integrals, it is useful if we can access the CGFs using
an index rather than directly by their variables. We do so by putting
the CGFs in a list like so

```python
# put the CGFs in a list so that we can use an iterator
# to access them
cgfs = [cgf1, cgf2]
```

To fill the matrices, one can use a twofold set of nested for loops
as follows
```python
for i in range(0,2):
    for j in range(0,2):
        S[i,j] = integrator.overlap(cgfs[i], cgfs[j])
```

## Question 1
Complete the script and calculate the overlap, kinetic energy and nuclear
attraction matrices.

## Two-electron integrals
For a basis set containing two basis functions, there are 16 distinct
combinations of two-electron integrals possible:
```
<11|11> <11|12> <11|21> <11|22>
<12|11> <12|12> <12|21> <12|22>
<21|11> <21|12> <21|21> <21|22>
<22|11> <22|12> <22|21> <22|22>
```

Yet there are only four unique values found due to the following
symmetry
```
<ab|cd> = <ab|dc> = <cd|ab>
```
and due to the fact that e.g. `<11|11>` is equal to `<22|22>` because of the symmetry of the molecule.

And as such
```
<11|11> = <22|22>
<11|22> = <22|11>
<21|11> = <12|11> = <11|12> = <11|21> = <12|22> = <21|22> = <22|21> = <22|12>
<21|21> = <12|12> = <21|12> = <12|21>
```

## Question 2
Calculate the 16 different two-electron integrals and show that there are only four unique
combinations. Use a fourfold nested for loop for this such as the following

```python
for i in range(0,2):
    for j in range(0,2):
        for k in range(0,2):
            for l in range(0,2):
                print(integrator.repulsion(cgfs[i],cgfs[j],cgfs[k],cgfs[l]))
```

## Solution
The solution is given in [ex03_solution.py](ex03_solution.py)
