# Hartree-Fock course for the Han-sur-Lesse Winterschool of 2021::Exercise 01

[![This work is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported License](https://i.creativecommons.org/l/by-nc/3.0/88x31.png)](http://creativecommons.org/licenses/by-nc/3.0/)

## Problem description
This exercise teaches you about contracted Gaussian functions and the various different integrals 
that can be constructed. We start simple by exploring a single H atom which has a single electron.

Contacted Gaussian functions (CGFs) are a combination of Gaussian Type Orbitals (GTOs). One can
build CGFs by specifying a series of coefficients that describe the GTOs. For example, for the
1s atomic orbital in H according to the STO-3g basis set, the following instructions build
the CGF.

```python
cgf1 = cgf([0.0, 0.0, 0.0])
cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)
```

To calculate the overlap, kinetic and nuclear attraction integral of this CGF with itself,
an integrator object needs to be constructed, which is done via the following one-liner

```python
integrator = PyQInt()
```

This integrator class has the following methods to evaluate the integrals
```python
s = integrator.overlap(<CGF1>, <CGF2>)
t = integrator.kinetic(<CGF1>, <CGF2>)
v = integrator.nuclear(<CGF1>, <CGF2>, <POSITION OF NUCLEUS>, <CHARGE OF NUCLEUS>)
```

## Questions
Calculate the overlap, kinetic and nuclear attraction integral of the 1s atomic orbital
as represented by the STO-3g basis set. Use [ex01_template.py](ex01_template.py) to get
started.

## Solution
The solution is given in [ex01_solution.py](ex01_solution.py)