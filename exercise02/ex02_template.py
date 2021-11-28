# -*- coding: utf-8 -*-

from pyqint import PyQInt, cgf

# construct the STO-3g CGF for He
cgf1 = cgf([0.0, 0.0, 0.0])
cgf1.add_gto(0.15432897000000001, 6.3624213899999997, 0, 0, 0)
cgf1.add_gto(0.53532813999999995, 1.1589229999999999, 0, 0, 0)
cgf1.add_gto(0.44463454000000002, 0.31364978999999998, 0, 0, 0)
pos = (0,0,0) # position of He nucleus

# create integrator object
integrator = PyQInt()

# start programming from here
