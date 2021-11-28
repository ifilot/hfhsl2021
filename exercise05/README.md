# Hartree-Fock course for the Han-sur-Lesse Winterschool of 2021::Exercise 05

[![This work is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported License](https://i.creativecommons.org/l/by-nc/3.0/88x31.png)](http://creativecommons.org/licenses/by-nc/3.0/)

## Problem description
In this exercise, we are going to complete the algorithm and build the full
Hartree-Fock algorithm. This is quite a daunting task if you do this for the
first time and as such, we use the concept of 'learning by example'. You are
given a fully working HF code for the H2 molecule and your task is to identify
all the steps in the procedure.

Before we will ook into the algorithm, let us first introduce a new way to
tackle the two-electron integrals. The fact of the matter is that the two-
electron integrals scale with N^4 where N is the number of CGFs. Some of
these two-electron integrals are duplicates with respect to others and as
such we want to avoid re-calculating these duplicates at all cost. To do so,
we introduce a new function, called `teindex` which assigns a value to
each set of four CGF indices in such a way that if two two-electron integrals
have the same value, they will have the same index. Note that this is purely
based on the CGF indices, i.e. if the molecule has some kind of inherent
symmetry, this function does **not** take that into account. The code below
shows the complete procedure.

```python
# build the two-electron integrals: because we want to avoid
# calculating a similar integral twice, we need some additional
# logics that allocates a unique index (or identifier) to each
# integral holding a unique value
N = len(cgfs)   # number of CGFs

# teint_calc keeps track of which integrals have already been calculated.
# if the integral is not yet evaluated, it holds a value of -1 for this
# integral, else it holds a value of 1
teint_calc = np.multiply(np.ones(integrator.teindex(N,N,N,N)), -1.0)

# all two-electron integral results are stored in a big list using
# a unique identifier as obtained from the .teindex function
teint = np.zeros(integrator.teindex(N,N,N,N))
for i, cgf1 in enumerate(cgfs):
    for j, cgf2 in enumerate(cgfs):
        ij = i*(i+1)/2 + j
        for k, cgf3 in enumerate(cgfs):
            for l, cgf4 in enumerate(cgfs):
                kl = k * (k+1)/2 + l
                if ij <= kl:
                    # determine unique identifier for each integral
                    idx = integrator.teindex(i,j,k,l)
                    if teint_calc[idx] < 0:
                        teint_calc[idx] = 1
                        teint[idx] = integrator.repulsion(cgfs[i], cgfs[j], cgfs[k], cgfs[l])
```

The other small matter to pay attention to is that we have introduced two
new variables to make keeping track of nuclei and electrons a bit more
convenient. These are given in the following code

```python
# define nuclei and number of electrons
nuclei = [
    [pos1, 1.0],
    [pos2, 1.0]
]
nelec = 2
```

## Question
The self-consistent field that is one of the major components of a Hartree-Fock
calculation (integral evaluation is the other), is schematically given in the
image below.

![Schematic depiction of the self-consistent field procedure](../img/scf_cycle.png)

Open the file [ex05_template.py](ex05_template.py) and indicate where in the
code each of the steps as shown in the image above are executed. Note that
`ex05_template.py` is already a fully working script, so feel free to execute
it. It is perhaps also insightful to print a few matrices to the terminal
for each of the iterations to see how everything is progressing.

## Solution
The solution is given in [ex05_solution.py](ex05_solution.py)