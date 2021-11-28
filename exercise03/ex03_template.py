# -*- coding: utf-8 -*-

from pyqint import PyQInt, cgf
from copy import deepcopy
import numpy as np

# construct the STO-3g CGF for H
cgf1 = cgf([0.0, 0.0, 0.0])
cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)
pos1 = (0,0,0) # position of first H nucleus

# create the second CGF for a H atom at 1.4 a.u. distance from the
# first one; the H-H bonding axis is aligned over the z-axis
pos2 = (0,0,1.4)        # position of second H nucleus
cgf2 = deepcopy(cgf1)
cgf2.p = pos2           # reset its center

# build 2x2 placeholders
S = np.zeros((2,2))     # overlap matrix
T = np.zeros((2,2))     # kinetic matrix
V1 = np.zeros((2,2))    # nuclear attraction matrix for H1
V2 = np.zeros((2,2))    # nuclear attraction matrix for H2

# build integrator object
integrator = PyQInt()

# put the CGFs in a list so that we can use an iterator
# to access them
cgfs = [cgf1, cgf2]
