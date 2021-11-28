# -*- coding: utf-8 -*-

from pyqint import PyQInt, cgf

# construct the STO-3g CGF for H
cgf1 = cgf([0.0, 0.0, 0.0])
cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)
pos = (0,0,0) # position of H nucleus

# create integrator object
integrator = PyQInt()

## start coding from here