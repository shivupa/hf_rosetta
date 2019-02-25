#!/usr/bin/env python
# imports
import numpy as np
import scipy.linalg as spla

# load data
convergence_DM = np.loadtxt("../../data/convergence_DM.txt")
convergence_E = np.loadtxt("../../data/convergence_E.txt")
S = np.loadtxt("../../data/S.txt")
T = np.loadtxt("../../data/T.txt")
V = np.loadtxt("../../data/V.txt")
eri = np.loadtxt("../../data/eri.txt")
E_nuc = np.loadtxt("../../data/E_nuc.txt")
iteration_max = (int)(np.loadtxt("../../data/iteration_max.txt"))
num_ao = (int)(np.loadtxt("../../data/num_ao.txt"))
num_elec_alpha = (int)(np.loadtxt("../../data/num_elec_alpha.txt"))
num_elec_beta = (int)(np.loadtxt("../../data/num_elec_beta.txt"))
iteration_max = (int)(np.loadtxt("../../data/iteration_max.txt"))
# Code


def idx2(i, j):
    if i >= j:
        return int(i*(i+1)/2+j)
    else:
        return int(j*(j+1)/2+i)


def idx4(i, j, k, l):
    return idx2(idx2(i, j), idx2(k, l))


D = np.zeros((num_ao, num_ao))
# loop variables
iteration_num = 0
E_total = 0
E_elec = 0.0
iteration_E_diff = 0.0
iteration_rmsc_dm = 0.0
converged = False
exceeded_iterations = False

H = T + V

while (not converged and not exceeded_iterations):
    # store last iteration and increment counters
    iteration_num += 1
    E_elec_last = E_elec
    D_last = np.copy(D)
    # form G matrix
    G = np.zeros((num_ao, num_ao))
    for i in range(num_ao):
        for j in range(num_ao):
            for k in range(num_ao):
                for l in range(num_ao):
                    G[i, j] += D[k, l] * ((2.0*(eri[idx4(i, j, k, l)])) -
                                          (eri[idx4(i, k, j, l)]))
    # build fock matrix
    F = H + G
    # solve the generalized eigenvalue problem
    E_orbitals, C = spla.eigh(F, S)
    # compute new density matrix
    D = np.zeros((num_ao, num_ao))
    for i in range(num_ao):
        for j in range(num_ao):
            for k in range(num_elec_alpha):
                D[i, j] += C[i, k] * C[j, k]
    # calculate electronic energy
    E_elec = np.sum(np.multiply(D, (H + F)))
    # calculate energy change of iteration
    iteration_E_diff = np.abs(E_elec - E_elec_last)
    # rms change of density matrix
    iteration_rmsc_dm = np.sqrt(np.sum((D - D_last)**2))
    if(np.abs(iteration_E_diff) < convergence_E and iteration_rmsc_dm < convergence_DM):
        converged = True
    if(iteration_num == iteration_max):
        exceeded_iterations = True

E_total = E_elec + E_nuc

print("{:^20.15f}".format(E_total))
