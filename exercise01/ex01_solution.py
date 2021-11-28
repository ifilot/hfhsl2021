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

# construct the STO-3g CGF for H
cgf1 = cgf([0.0, 0.0, 0.0])
cgf1.add_gto(0.154329, 3.425251, 0, 0, 0)
cgf1.add_gto(0.535328, 0.623914, 0, 0, 0)
cgf1.add_gto(0.444635, 0.168855, 0, 0, 0)
pos = (0,0,0) # position of H nucleus

# create integrator object
integrator = PyQInt()

# calculate integrals
s = integrator.overlap(cgf1, cgf1)
t = integrator.kinetic(cgf1, cgf1)
v = integrator.nuclear(cgf1, cgf1, pos, 1.0)

# print results to terminal
print('Overlap: %f' % s)
print('Kinetic: %f' % t)
print('Coulombic interaction: %f' % v)