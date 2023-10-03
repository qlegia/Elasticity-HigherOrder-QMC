# generate all pbs scripts using gen_sp_pbs2
# usage:
#  python gen_all_sp2.py > all_sp_pbs2.sh
#  source all_sp_pbs2.sh
L = 15
maxk = 32 
for k in range(0,maxk):
    p0 = k*4+1
    p1 = (k+1)*4
    print("python gen_sp_pbs2.py %d %d %d > s2L%d_qmc%d.pbs" % (L,p0,p1,L,k))

