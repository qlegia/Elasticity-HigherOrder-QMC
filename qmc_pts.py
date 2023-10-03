import sys
from pyQMC.lattice import PolynomialLattice

# number of dimensions
s = 32 

m = 4
# name of file containing generating vector
qmc_dir = "./pyQMC-master1"
filename = qmc_dir+"/staircase2d_spod_a2_C0.1_SPOD_2dstaircase_t0/SPOD_2dstaircase_t0_m"+str(m)+".json"
lat = PolynomialLattice(filename, s)
for p in lat:
    print (p)


