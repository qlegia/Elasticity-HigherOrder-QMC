# generate all pbs scripts using gen_e_pbs2
# usage:
#  python gen_all_pbs2.py > script_all_pbs.sh
#  source script_all_pbs.sh
for k in range(0,102):
    p0 = k*20+1
    p1 = (k+1)*20
    print("python gen_e_pbs2.py %d %d > f_qmc%d.pbs" % (p0,p1,k))
print("python gen_e_pbs2.py 2041 2048 > f_qmc102.pbs")
