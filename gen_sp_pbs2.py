# scripts for katana
# usage
#   python gen_sp_pbs2.py L p_min p_max
import sys
L  = int(sys.argv[1])
p_min  = int(sys.argv[2])
p_max  = int(sys.argv[3])
print("#!/bin/bash ")
print("#PBS -N E_p%d_%d " % (p_min,p_max))  
print("#PBS -l select=1:ncpus=2:mem=10gb ")
print("#PBS -l walltime=24:00:00 ")
print("#PBS -m ae ")
print("#PBS -M qlegia@unsw.edu.au ")
print("#PBS -J 1-%d \n" % L)
print("cd /home/z9701564/ElasticityQMC ")
print("echo 'I am now working on job ${PBS_ARRAY_INDEX}' ")
print("conda init bash")
print("source activate fenicsx-env")
for p in range(p_min,p_max+1):
  print("python3 sparse_qmc2.py %d ${PBS_ARRAY_INDEX} %d" % (L,p))

