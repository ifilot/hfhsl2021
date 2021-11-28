# Hartree-Fock course for the Han-sur-Lesse Winterschool of 2021::Exercise 02

[![This work is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported License](https://i.creativecommons.org/l/by-nc/3.0/88x31.png)](http://creativecommons.org/licenses/by-nc/3.0/)

## Problem description
We expand upon the previous exercise and start to look into a He atom. He contains two electrons
and as such, a two-electron integral occurs corresponding to the repulsion between these two
electrons. In this exercise, you will learn how to evaluate a two-electron integral and
calculate the total electronic energy for the He atom.

Every atom uses slightly different coefficients for their CGFs and GTOs which are tuned to
optimally describe the atom when it is bound to other atoms, i.e., when it resides
in a molecule. For the He atom, the coefficients to build the CGF are as follows.

```python
cgf1 = cgf([0.0, 0.0, 0.0])
cgf1.add_gto(0.15432897000000001, 6.3624213899999997, 0, 0, 0)
cgf1.add_gto(0.53532813999999995, 1.1589229999999999, 0, 0, 0)
cgf1.add_gto(0.31364978999999998, 0.44463454000000002, 0, 0, 0)
```

Helium has a single pair of electrons at opposite spin. Thus, it does have a
two-electron integral corresponding to the repulsion between these two electrons but
it does not have an exchange integral. To evaluate a two-electron integral, you
can use the following method.

```python
te = integrator.repulsion(<CGF1>, <CGF2>, <CGF3>, <CGF4>)
```

With two electrons in the same spatial orbital, the total energy for the He
atom is given by:
```
e_total = 2.0 * (t + v) + te
```

wherein `t` is the kinetic energy, `v` is the nuclear attraction, and `te` is the
electron-electron repulsion energy.

## Questions
Calculate the electron-electron repulsion energy of He and calculate its ground
state energy. Take a good look at the workout of exercise 01 to get you started.
Do not forget that the nuclear charge of He is 2 and not 1 as is the case of H.

## Solution
The solution is given in [ex02_solution.py](ex02_solution.py)
