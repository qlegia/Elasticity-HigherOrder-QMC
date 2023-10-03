# usage
#  python gen_sh.py > qsub_all.sh
#  source qsub_all.sh
for k in range(3,51):
    print("qsub qmc%d.pbs" % k)
