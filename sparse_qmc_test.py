# Adapted from
# https://github.com/jorgensd/dolfinx-tutorial/commit/e3669f2fe760748cbe79ff26d6bec3dd5f7251c6#
# to run the file, 
#  conda activate fenicsx-env
#  python3 sparse_qmc L k q 
#  python3 sparse_qmc 
import numpy as np

import sys 

# implement the sparse grid sum (6.6)
# N_1^{(j)} = 2^{L-k},  L-k=9,...,1
# N_2^{(k)} = 2^k,      k=1,...,9
#
# compute
#   SN = 2^L Q_{N^{L-k,k}} = sum_{j=1}^{2^{L-k}} sum_{i=1}^{2^k} calF(p(j), q(i))
# in partion qq for qq = 1,2,3,4 for L=10,
# let cnt be the counter run over 2 loops of p(j) and q(i) 
#  partition 1 = { 0 < cnt <=256}  
#  partition 2 = { 256 < cnt <= 512}  
#  partition 3 = { 512 < cnt <= 512+256=768}  
#  partition 4 = { 768 < cnt <= 1024}  
#  
# each k will be run in parallel, see sparse_qmc.pbs 
LL = int(sys.argv[1]) # first parameter from the command line
kk = int(sys.argv[2]) # second parameter from the command line
qq = int(sys.argv[3]) # third parameter from the command line

from pyQMC.lattice import PolynomialLattice
# QMC points
dim_s = 16   # number of truncated dimensions
print('dims=',dim_s)
# name of file containing generating vector
qmc_dir = "./pyQMC-master1"
N1 = 2**(LL-kk)  # LL=10 so LL-kk=9,8,7,6,5,4,3,2,1
N2 = 2**kk #                    k=1,2,3,4,5,6,7,8,9
print('N1= ',N1, ' N2= ', N2)
filename1 = qmc_dir+"/staircase2d_spod_a2_C0.1_SPOD_2dstaircase_t0/SPOD_2dstaircase_t0_m"+str(LL-kk)+".json"
# lat is a list containing all the QMC quadrature points
lat1 = PolynomialLattice(filename1, dim_s)
filename2 = qmc_dir+"/staircase2d_spod_a2_C0.1_SPOD_2dstaircase_t0/SPOD_2dstaircase_t0_m"+str(kk)+".json"
lat2 = PolynomialLattice(filename2, dim_s)

k1 = np.arange(1,dim_s)  #{1,2,3,4,5,...,s}
k2 = np.arange(1,2*dim_s,2)  # {1,3,5,7,9,...,2s+1}
step = 256



from math import sin, pi

cnt = 0  # counter for both p and q
SN = 0.0
cmin = 2**LL
cmax = 0
for qmc_p in lat1:
      #
      for qmc_q in lat2:
         cnt = cnt + 1
         if ((cnt > (qq-1)*step)&(cnt <= step*qq)):
             if (cnt < cmin):
                  cmin = cnt 
             if (cnt > cmax):
                  cmax = cnt
print('SN= ',SN,' cmin=', cmin, ' cmax=',cmax)
