# generate all pbs scripts using gen_e_pbs
# usage:
#  python gen_all_pbs.py > script_all_pbs.sh
#  source script_all_pbs.sh
for k in range(3,50):
    p0 = k*20+1
    p1 = (k+1)*20
    print("python gen_e_pbs.py %d %d > qmc%d.pbs" % (p0,p1,k))
print("python gen_e_pbs.py 1001 1024 > qmc50.pbs")
