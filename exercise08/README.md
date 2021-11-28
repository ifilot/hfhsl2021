# Hartree-Fock course for the Han-sur-Lesse Winterschool of 2021::Exercise 06

[![This work is licensed under a Creative Commons Attribution-NonCommercial 3.0 Unported License](https://i.creativecommons.org/l/by-nc/3.0/88x31.png)](http://creativecommons.org/licenses/by-nc/3.0/)

## Problem description
In this exercise, we will visualize the molecular orbitals using
contour plots. Our Python script will generate a .abo file which
can be opened in the Managlyph program, allowing you to easily
browse through the molecular orbitals.

Make sure you have downloaded and installed [Managlyph](https://www.managlyph.nl/download).

In this exercise, you will not do a lot of programming as the programming
involved would be rather complex. We mainly want to explain the process
of isosurface creation to you so that you obtain conceptual understanding
of how this works.

Isosurface creation is basically connecting all points in a three-dimensional space
that have the same value, called the isovalue. There are several algorithms
available which can perform the surface creation but one of the most commonly
used and also one of the most efficient ones is the [marching cubes](https://en.wikipedia.org/wiki/Marching_cubes)
algorithm. Programming this algorithm yourself is moderately challenging
and as such, we readily make use of the [PyTessel](https://github.com/ifilot/pytessel) to
do the heavy lifting for us.

Three-dimensional surfaces can be represented by a myriad of files. In a very creative
scenario, you could even print such surfaces on a 3D-printer. In this exercise,
we have chosen the `.abo` format, which is a native (though open-source) format used
by Managlyph. An `.abo` file can contain the atoms of a molecule as well as a number
of molecular orbitals. Furthermore, it is relatively compact and thus ideal if you want
to share it with your colleagues. Before showing you how to build an `.abo` file from
your HF calculation, let me briefly mention here that PyTessel also supports outputting
`.ply` files. These files are very common and supported by a broad set of programs. Creating
`.ply` files is the way to go if you want to use your molecular orbitals in another program
such as [Blender](https://www.blender.org/) for rendering or any of the [slicer](https://en.wikipedia.org/wiki/Slicer_(3D_printing))
programs for 3D-printing. Further information on how to do this is found in the README.md
files of [PyQInt](https://github.com/ifilot/pyqint) and [PyTessel](https://github.com/ifilot/pytessel).

Without further ado, to create a `.abo` file, a handy function is provided as 
shown below. Admittedly, this is a moderately complex function but using the function
should be very simple. The arguments are `outfile` which specifies the location
where you want to store your `.abo` file, `nuclei` is a list of nuclei, `cgfs` 
is a list of contracted Gaussian functions, `coeff` is the coefficient matrix,
`energies` is a list of energies and finally `isovalue` is the isovalue to use
for isosurfaces (0.1 is a good value). Note that `nuclei`, `cgfs`, `coeff` and
`energies` are obtained from the previous HF calculation.

```python
build_abo(outfile, nuclei, cgfs, coeff, energies, isovalue):
```

[An example script](ex08_solution.py) where a HF calculation of CO and the generation of the isosurfaces of
its molecular orbitals is already given to you. Open the script in Spyder and execute it. In the
same folder as where you have found `ex08_solution.py`, you should now find a file called `co.abo`. Open
this file in Managlyph. You should see a result as follows.

![Molecular orbital visualized using Managlyph](../img/managlyph_co.png)

Using the arrows, you can browser between the molecular orbitals. Take
a look at the results and compare them with the contour plots as obtained
in the previous exercise.

## Exercises

1. Compare the isosurfaces with the contour plots. What are the advantages and disadvantages of 
   visualizing molecular orbitals using one versus the other.
2. How are the double-degenerate molecular orbitals related to each other?

## Solutions
1. Isosurfaces are great for showing the three-dimensional shape and do not
   depend on the choice of surface on which you project the values. That being
   said, they have a hard time conveying the internal structure of a molecular orbital.
   For example, if there would be concentric nodal sphere inside the molecular orbital,
   e.g. as in any of the 2s AOs, if would not be visible. Such features would be
   readily visible in contour plots. In other words, there is no silver bullet
   and you are recommended to use both if you want to get a more holistic picture.
2. The double-degenerate molecular orbitals are related to each other by a 90 degree
   rotation around the C-O bonding axis (here, the z-axis). This 90 degree rotation
   almost immediately shows that despite these orbitals lying at the same energy, they
   remain orthogonal with respect to each other.

## Further details
Let me share a bit more information on building isosurfaces. Isosurfaces are the
three-dimensional surfaces that connect all points with the same isovalue. For
the wave functions of the molecular orbitals, we want to visualize both the
positive and negative lobes. Thus instead of building a single isosurface,
we have to build two, one for each sign. One handy feature of the `.abo` format
is that it allows to store multiple 3D objects together with a molecule. Other
programs, including Blender, would require you to import two objects. That
is of course not too big of an obstacle, but you would like to streamline the
process as much as possible, which is why I opted for this format.
