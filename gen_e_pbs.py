# generate pbs scripts for katana
# usage
#   python gen_e_pbs.py p_min p_max
import sys
p_min  = int(sys.argv[1])
p_max  = int(sys.argv[2])
print("#!/bin/bash ")
print("#PBS -N E_p%d_%d " % (p_min,p_max))  
print("#PBS -l select=1:ncpus=2:mem=60gb ")
print("#PBS -l walltime=24:00:00 ")
print("#PBS -m ae ")
print("#PBS -M qlegia@unsw.edu.au ")
print("#PBS -J %d-%d \n"%(p_min,p_max))
print("cd /home/z9701564/ElasticityQMC ")
print("echo 'I am now working on job ${PBS_ARRAY_INDEX}' ")
print("conda init bash")
print("source activate fenicsx-env")
print("python3 par_elas_qmc.py ${PBS_ARRAY_INDEX} 1")
print("python3 par_elas_qmc.py ${PBS_ARRAY_INDEX} 2")
print("python3 par_elas_qmc.py ${PBS_ARRAY_INDEX} 3")
print("python3 par_elas_qmc.py ${PBS_ARRAY_INDEX} 4")
print("python3 par_elas_qmc.py ${PBS_ARRAY_INDEX} 5")
