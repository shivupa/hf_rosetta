import pyscf
from pyscf import gto, scf
import numpy as np

with open("../eqh2o.xyz",'r') as f:
    atomxyz = ""
    for count,line in enumerate(f):
        if count > 1:
            atomxyz += "{};".format(line.strip("\n"))
        
mol = pyscf.gto.M(
    atom = atomxyz,
    basis = 'cc-pvdz',
    unit="Ang",
    verbose=0,
    symmetry=False,
    spin=0,
    charge=0
)

# Nuclear repulsion energy
E_nuc = mol.energy_nuc()
E_nuc = np.array([E_nuc])
np.savetxt("E_nuc.txt",E_nuc)
# calculate overlap integrals
S = mol.intor('cint1e_ovlp_sph')
np.savetxt("S.txt",S)
# calculate kinetic energy integrals
T = mol.intor('cint1e_kin_sph')
np.savetxt("T.txt",T)
# calculate nuclear attraction integrals
V = mol.intor('cint1e_nuc_sph')
np.savetxt("V.txt",V)
# calculate two electron integrals
eri = mol.intor('cint2e_sph', aosym='s8')
np.savetxt("eri.txt",eri)
