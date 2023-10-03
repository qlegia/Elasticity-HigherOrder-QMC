# generate all pbs scripts using gen_sp_pbs
# usage:
#  python gen_all_sp.py > all_sp_pbs.sh
#  source all_sp_pbs.sh
# (L,maxk) = (13, 8); (14, 16); (15, 32)
#  eg. for k in range(0,16): 
#     b/c for L=15 , 2^15/256/8 = 16 (2^L/step/ 8 partition per pbs)   
L = 15
maxk = 32 
for k in range(0,maxk):
    p0 = k*4+1
    p1 = (k+1)*4
    print("python gen_sp_pbs.py %d %d %d > sL%d_qmc%d.pbs" % (L,p0,p1,L,k))

