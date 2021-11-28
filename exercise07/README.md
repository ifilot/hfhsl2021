# Hartree-Fock course for the Han-sur-Lesse Winterschool of 2021::Exercise 06

[![This work is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported License](https://i.creativecommons.org/l/by-nc/3.0/88x31.png)](http://creativecommons.org/licenses/by-nc/3.0/)

## Problem description
There are two common ways to visualize molecular orbitals and these correspond
to projecting the wave function onto a plane (contour plot) or by connecting
all points in three-dimensional space to form a surface (isosurface). In this
exercise, we will explore how to build contour plots.

Once the self-consistent field procedure is converged, we can readily
visualize the results. The coefficient matrix contains the proportionality
constants, i,e. the linear coefficients, how much each of the basis functions
contribute to the molecular orbital. As such, by looping over these coefficients
and the basis functions, we can calculate the value of the wave function of
the molecular orbital as function of a three-dimensional coordinate. Evaluation
of the value of the wave function is something that is already encoded in the
`integrator` class of `PyQInt` and thus we can readily use this function. But
before we can use this function, we also need to tell the function at which
grid points we wish to evaluate the wave function.

The code below demonstrates this procedure. First a rectangular equidistant
grid of 100x100 points between [-2,2] on the xz plane is built. Next, the
values for the wave function is evaluated at each of these grid points. Finally,
a heat-map of the values projected onto the plane is produced.

```python
# produce an isosurface of the molecular orbitals
x = np.linspace(-2, 2, 100)
z = np.linspace(-2, 2, 100)
xx, zz = np.meshgrid(x,z)
yy = np.zeros(len(x) * len(z))
grid = np.vstack([xx.flatten(), yy, zz.flatten()]).reshape(3,-1).T

plt.figure(dpi=300)
res = integrator.plot_wavefunction(grid, C[:,0], cgfs).reshape((len(z), len(x)))
plt.imshow(res, origin='lower', extent=[-2,2,-2,2], cmap='PiYG')
plt.colorbar()
plt.xlabel('x [a.u.]')
plt.ylabel('z [a.u.]')
plt.title('Contour plot of the MO #1 of CO')
plt.show()
plt.close()
```

## Questions
Produce a heat-map for each of the molecular orbitals of CO by modifying the
[template file](ex07_template.py). Your results should look 
something like this (though not necessarily in one image).

![Molecular orbitals of CO](../img/mos_co.png)

## Solution
The solution is given in [ex07_solution.py](ex07_solution.py)
