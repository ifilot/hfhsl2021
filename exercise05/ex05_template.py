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

# define nuclei and number of electrons
nuclei = [
    [pos1, 1.0],
    [pos2, 1.0]
]
nelec = 2

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

# add the two nuclear attraction matrices
V = V1 + V2 # note, this is elementwise addition

# build the two-electron integrals: because we want to avoid
# calculating a similar integral twice, we need some additional
# logics that allocates a unique index (or identifier) to each
# integral holding a unique value
N = len(cgfs)   # number of CGFs

# teint_calc keeps track of which integrals have already been calculated.
# if the integral is not yet evaluated, it holds a value of -1 for this
# integral, else it holds a value of 1
teint_calc = np.multiply(np.ones(integrator.teindex(N,N,N,N)), -1.0)

# all two-electron integral results are stored in a big list using
# a unique identifier as obtained from the .teindex function
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

# diagonalize S
s, U = np.linalg.eigh(S)

# construct transformation matrix X
X = U.dot(np.diag(1.0/np.sqrt(s)))

# create empty P matrix as initial guess
P = np.zeros(S.shape)

# start iterative procedure; it is always good practice to set an
# upper bound to the number of cycles (here: 100)
energies = []
for niter in range(0,100):
    G = np.zeros(S.shape)
    for i in range(S.shape[0]):
        for j in range(S.shape[0]):
            for k in range(S.shape[0]):
                for l in range(S.shape[0]):
                    # here, we first establish the index for each required
                    # two-electron integral by which we can obtain its value
                    # from the list that was built earlier
                    idx_rep = integrator.teindex(i,j,l,k)
                    idx_exc = integrator.teindex(i,k,l,j)
                    G[i,j] += P[k,l] * (teint[idx_rep] - 0.5 * teint[idx_exc])
    
    # build Fock matrix
    F = T + V + G
    
    # transform Fock matrix
    Fprime = X.transpose().dot(F).dot(X)
    
    # diagonalize F
    e, Cprime = np.linalg.eigh(Fprime)
    
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
    
    # calculate a new P
    P = np.zeros(S.shape)
    for i in range(S.shape[0]):
        for j in range(S.shape[0]):
            for k in range(0,int(nelec/2)):
                P[i,j] += 2.0 * C[i,k] * C[j,k]

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