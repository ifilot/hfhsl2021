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

from pyqint import PyQInt, Molecule
import numpy as np
import matplotlib.pyplot as plt

############################################
#
# STEP 1: Define nuclei and basis functions
#
############################################

# build molecule
mol = Molecule('CO')
mol.add_atom('C', 0.0, 0.0, -2.116/2.0)
mol.add_atom('O', 0.0, 0.0,  2.116/2.0)

cgfs, nuclei = mol.build_basis('sto3g')
nelec = np.sum([n[1] for n in nuclei])
N = len(cgfs)

# build 2x2 placeholders
S = np.zeros((N,N))     # overlap matrix
T = np.zeros((N,N))     # kinetic matrix
V1 = np.zeros((N,N))    # nuclear attraction matrix for H1
V2 = np.zeros((N,N))    # nuclear attraction matrix for H2

# build integrator object
integrator = PyQInt()

############################################
#
# STEP 2: Calculate S,T,V,H,TEINT integrals
#
############################################

# calculate the overlap, kinetic energy and nuclear attraction matrices
for i in range(0,N):
    for j in range(0,N):
        S[i,j] = integrator.overlap(cgfs[i], cgfs[j])
        T[i,j] = integrator.kinetic(cgfs[i], cgfs[j])
        V1[i,j] = integrator.nuclear(cgfs[i], cgfs[j], nuclei[0][0], nuclei[0][1])
        V2[i,j] = integrator.nuclear(cgfs[i], cgfs[j], nuclei[1][0], nuclei[1][1])

# add the two nuclear attraction matrices
V = V1 + V2 # note, this is elementwise addition

# calculate two-electron integrals
teint_calc = np.multiply(np.ones(integrator.teindex(N,N,N,N)), -1.0)
teint = np.zeros(integrator.teindex(N,N,N,N))
for i, cgf1 in enumerate(cgfs):
    for j, cgf2 in enumerate(cgfs):
        ij = i*(i+1)/2 + j
        for k, cgf3 in enumerate(cgfs):
            for l, cgf4 in enumerate(cgfs):
                kl = k * (k+1)/2 + l
                if ij <= kl:
                    # determine unique identifier for each integral
                    idx = integrator.teindex(i,j,k,l)
                    if teint_calc[idx] < 0:
                        teint_calc[idx] = 1
                        teint[idx] = integrator.repulsion(cgfs[i], cgfs[j], cgfs[k], cgfs[l])

############################################
#
# STEP 3: Calculate transformation matrix
#
############################################

# diagonalize S
s, U = np.linalg.eigh(S)

# construct transformation matrix X
X = U.dot(np.diag(1.0/np.sqrt(s)))

#################################################
#
# STEP 4: Obtain initial guess for density matrix
#
#################################################

# create empty P matrix as initial guess
P = np.zeros(S.shape)

# start iterative procedure; it is always good practice to set an
# upper bound to the number of cycles (here: 100)
energies = []
for niter in range(0,100):
    #################################################
    #
    # STEP 5: Calculate G,H,F,F' from P
    #
    #################################################
    
    G = np.zeros(S.shape)
    for i in range(S.shape[0]):
        for j in range(S.shape[0]):
            for k in range(S.shape[0]):
                for l in range(S.shape[0]):
                    idx_rep = integrator.teindex(i,j,l,k)
                    idx_exc = integrator.teindex(i,k,l,j)
                    G[i,j] += P[k,l] * (teint[idx_rep] - 0.5 * teint[idx_exc])
    
    # build Fock matrix
    F = T + V + G
    
    # transform Fock matrix
    Fprime = X.transpose().dot(F).dot(X)
    
    #################################################
    #
    # STEP 6: Diagonalize F' to obtain C' and e
    #
    #################################################
    
    # diagonalize F
    e, Cprime = np.linalg.eigh(Fprime)
    
    #################################################
    #
    # STEP 7: Calculate C from C'
    #
    #################################################
    
    # back-transform
    C = X.dot(Cprime)
    
    # calculate energy E
    energy = 0.0
    M = T + V + F
    for i in range(S.shape[0]):
        for j in range(S.shape[0]):
            energy += 0.5 * P[j,i] * M[i,j]
    
    # calculate repulsion of the nuclei
    for i in range(0, len(nuclei)):
        for j in range(i+1, len(nuclei)):
            r = np.linalg.norm(np.array(nuclei[i][0]) - np.array(nuclei[j][0]))
            energy += nuclei[i][1] * nuclei[j][1] / r
    
    #################################################
    #
    # STEP 8: Calculate P from C
    #
    #################################################
    
    # calculate a new P
    P = np.zeros(S.shape)
    for i in range(S.shape[0]):
        for j in range(S.shape[0]):
            for k in range(0,int(nelec/2)):
                P[i,j] += 2.0 * C[i,k] * C[j,k]

    #################################################
    #
    # PRINT INFORMATION OF THIS ITERATION AND CHECK
    # FOR convergence
    #
    #################################################

    # print info for this iteration
    print("Iteration: %i Energy: %f" % (niter, energy))
    
    # calculate energy difference between this and the previous
    # iteration; terminate the loop when energy difference is less
    # than threshold;
    # note that convergence is here based purely on the energies,
    # alternatively, it can be based on the values of the density
    # matrix
    if niter > 1:
        ediff = np.abs(energy - energies[-1])
        if ediff < 1e-5:
            print("Stopping SCF cycle, convergence reached.")
            break
    
    # store energy for next iteration
    energies.append(energy)
    
# produce an isosurface of the molecular orbitals
x = np.linspace(-2, 2, 100)
z = np.linspace(-2, 2, 100)
xx, zz = np.meshgrid(x,z)
yy = np.zeros(len(x) * len(z))
grid = np.vstack([xx.flatten(), yy, zz.flatten()]).reshape(3,-1).T

for i in range(0, len(e)):
    plt.figure(dpi=300)
    res = integrator.plot_wavefunction(grid, C[:,i], cgfs).reshape((len(z), len(x)))
    plt.imshow(res, origin='lower', extent=[-2,2,-2,2], cmap='PiYG')
    plt.colorbar()
    plt.xlabel('x [a.u.]')
    plt.ylabel('z [a.u.]')
    plt.title('Contour plot of the MO #%i of CO' % (i+1))
    plt.show()
    plt.close()