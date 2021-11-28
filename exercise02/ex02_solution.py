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

# construct the STO-3g CGF for He
cgf1 = cgf([0.0, 0.0, 0.0])
cgf1.add_gto(0.15432897000000001, 6.3624213899999997, 0, 0, 0)
cgf1.add_gto(0.53532813999999995, 1.1589229999999999, 0, 0, 0)
cgf1.add_gto(0.44463454000000002, 0.31364978999999998, 0, 0, 0)
pos = (0,0,0) # position of He nucleus

# create integrator object
integrator = PyQInt()

# calculate integrals
s = integrator.overlap(cgf1, cgf1)
t = integrator.kinetic(cgf1, cgf1)
v = integrator.nuclear(cgf1, cgf1, pos, 2.0)
te = integrator.repulsion(cgf1, cgf1, cgf1, cgf1)

e_total = 2.0 * (t + v) + te

# print results to terminal
print('Overlap: %f' % s)
print('Kinetic: %f' % t)
print('Coulombic interaction: %f' % v)
print('Electron-electron repulsion: %f' % te)
print('Total energy: %f' % e_total)
