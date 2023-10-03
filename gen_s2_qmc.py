# usage
#  python gen_s2_qmc.py > qsub_all_sp2.sh
#  source qsub_all_sp2.sh
for k in range(0,16):
    print("qsub s2_qmc%d.pbs" % k)
