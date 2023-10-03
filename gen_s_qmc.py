# usage
#  python gen_s_qmc.py > qsub_all_sp.sh
#  source qsub_all_sp.sh
L = 15
maxk = 32 
for k in range(0,maxk):
    print("qsub sL%d_qmc%d.pbs" % (L,k))
    print("qsub s2L%d_qmc%d.pbs" % (L,k))
