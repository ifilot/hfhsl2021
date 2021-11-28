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
from mendeleev import element
from pytessel import PyTessel
import os

def main():
    nuclei, cgfs, coeff, energies = calculate_ch4()
    build_abo('ch4.abo', nuclei, cgfs, coeff, energies, 0.1)

def build_abo(outfile, nuclei, cgfs, coeff, energies, isovalue):
    """
    Build managlyph atom/bonds/orbitals file from
    previous HF calculation
    """
    # build integrator
    integrator = PyQInt()
    
    # set transparency
    alpha = 0.97
    
    # specify colors for occupied and virtual orbitals
    colors = [
        np.array([0.592, 0.796, 0.369, alpha], dtype=np.float32),
        np.array([0.831, 0.322, 0.604, alpha], dtype=np.float32),
        np.array([1.000, 0.612, 0.000, alpha], dtype=np.float32),
        np.array([0.400, 0.831, 0.706, alpha], dtype=np.float32)
    ]
    
    # build pytessel object
    pytessel = PyTessel()

    # build output file
    f = open(outfile, 'wb')

    # write number of frames
    nr_frames = len(cgfs) + 1
    f.write(nr_frames.to_bytes(2, byteorder='little'))

    #
    # First write the bare geometry of the molecule
    #

    # write frame_idx
    f.write(int(1).to_bytes(2, byteorder='little'))

    descriptor = 'Geometry'

    f.write(len(descriptor).to_bytes(2, byteorder='little'))
    f.write(bytearray(descriptor, encoding='utf8'))

    # write nr_atoms
    f.write(len(nuclei).to_bytes(2, byteorder='little'))
    for atom in nuclei:
        f.write(element(atom[1]).atomic_number.to_bytes(1, byteorder='little'))
        f.write(np.array(atom[0], dtype=np.float32).tobytes())

    # write number of models
    f.write(int(0).to_bytes(1, byteorder='little'))
    f.write(int(0).to_bytes(1, byteorder='little'))

    # calculate number of electrons
    nelec = np.sum([atom[1] for atom in nuclei])

    #
    # Write the geometry including the orbitals
    #
    for i in range(1, nr_frames):
        # write frame_idx
        f.write((i+1).to_bytes(2, byteorder='little'))

        descriptor = 'Molecular orbital %i\nEnergy: %.4f eV' % (i,energies[i-1])

        f.write(len(descriptor).to_bytes(2, byteorder='little'))
        f.write(bytearray(descriptor, encoding='utf8'))

        # write nr_atoms
        f.write(len(nuclei).to_bytes(2, byteorder='little'))
        for atom in nuclei:
            f.write(element(atom[1]).atomic_number.to_bytes(1, byteorder='little'))
            f.write(np.array(atom[0], dtype=np.float32).tobytes())

        print('Writing MO #%02i' % i)

        # write number of models
        f.write(int(2).to_bytes(2, byteorder='little'))
        for j in range(0, 2):
            # build the pos and negative isosurfaces from the cubefiles
            sz = 100
            grid = integrator.build_rectgrid3d(-5, 5, sz)
            scalarfield = np.reshape(integrator.plot_wavefunction(grid, coeff[:,i-1], cgfs), (sz, sz, sz))
            unitcell = np.diag(np.ones(3) * 10.0)
            vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), isovalue if j==1 else -isovalue)
            vertices_normals = np.hstack([vertices, normals])
            
            # write model idx
            f.write(j.to_bytes(2, byteorder='little'))
            
            # write model color
            if i < nelec / 2:
                color = np.array(colors[j])
            else:
                color = np.array(colors[j+2])
            f.write(color.tobytes())
            
            # write number of vertices
            f.write(vertices_normals.shape[0].to_bytes(4, byteorder='little'))
            
            # write vertices
            f.write(vertices_normals.tobytes())
            
            # write number of indices
            f.write(int(len(indices)/3).to_bytes(4, byteorder='little'))
            
            # write indices
            f.write(indices.tobytes())
            
            if j == 0:
                print('    Writing positive lobe: %i vertices and %i facets' % (vertices_normals.shape[0], indices.shape[0] / 3))
            else:
                print('    Writing negative lobe: %i vertices and %i facets' % (vertices_normals.shape[0], indices.shape[0] / 3))

    f.close()

    # report filesize
    print("Creating file: %s" % outfile)
    print("Size: %f MB" % (os.stat(outfile).st_size / (1024*1024)))

def calculate_ch4():
    ############################################
    #
    # STEP 1: Define nuclei and basis functions
    #
    ############################################
    
    # build molecule
    dist = np.power(1.09, 1.0/3.0)
    mol = Molecule('CH4')
    mol.add_atom('C', 0.00000000, 0.00000000, 0.0),
    mol.add_atom('H',  dist,  dist,  dist)
    mol.add_atom('H', -dist, -dist,  dist)
    mol.add_atom('H', -dist,  dist, -dist)
    mol.add_atom('H',  dist, -dist, -dist)
    
    cgfs, nuclei = mol.build_basis('sto3g')
    nelec = np.sum([n[1] for n in nuclei])
    N = len(cgfs)
    
    # build 2x2 placeholders
    S = np.zeros((N,N))     # overlap matrix
    T = np.zeros((N,N))     # kinetic matrix
    V = np.zeros((N,N))     # nuclear attraction matrix for H1
    
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
            for k in range(0, len(nuclei)):
                V[i,j] += integrator.nuclear(cgfs[i], cgfs[j], nuclei[k][0], nuclei[k][1])
    
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
    
    return nuclei, cgfs, C, e

if __name__ == '__main__':
    main()