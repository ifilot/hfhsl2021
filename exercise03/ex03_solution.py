# -*- coding: utf-8 -*-

# 
# This file is part of the HFHSL2021 distribution (https://github.com/ifilot/hfhsl2021).
# Copyright (c) 2021 Ivo Filot <i.a.w.filot@tue.nl>
# 
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

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

# calculate the overlap, kinetic energy and nuclear attraction matrices
for i in range(0,2):
    for j in range(0,2):
        S[i,j] = integrator.overlap(cgfs[i], cgfs[j])
        T[i,j] = integrator.kinetic(cgfs[i], cgfs[j])
        V1[i,j] = integrator.nuclear(cgfs[i], cgfs[j], pos1, 1.0)
        V2[i,j] = integrator.nuclear(cgfs[i], cgfs[j], pos2, 1.0)

print('Overlap matrix:\n', S)
print('Kinetic energy matrix:\n', T)
print('Nuclear attraction matrix #1:\n', V1)
print('Nuclear attraction matrix #2:\n', V2)

# calculate two-electron integrals and automatically produce labels
# to identify them
tes = []
labels = []
for i in range(0,2):
    for j in range(0,2):
        for k in range(0,2):
            for l in range(0,2):
                tes.append(integrator.repulsion(cgfs[i],cgfs[j],cgfs[k],cgfs[l]))
                labels.append('(%i,%i,%i,%i)' % (i+1, j+1, k+1, l+1))

# print all two-electron integrals
print('Two electron integrals list:')
for label, te in zip(labels, tes):
    print(label, te)
    
# Note: the code below is a bit more advanced and is only meant to
# show you which two-electron integrals are equal
tedict = {}
for label, te in zip(labels, tes):
    value = '%0.4f' % te
    if value not in tedict:
        tedict[value] = []
    tedict[value].append(label)

print('The following integrals have the same values:')
for key in tedict.keys():
    print(tedict[key], key)