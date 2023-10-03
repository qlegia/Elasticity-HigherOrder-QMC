# usage
#  python gen_sh2.py > qsub_all.sh
#  source qsub_all.sh
for k in range(0,103):
    print("qsub f_qmc%d.pbs" % k)
