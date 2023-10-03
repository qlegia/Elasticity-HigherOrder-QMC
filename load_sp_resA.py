# load results
import numpy as np
fdir = './sp_resL11/'
fdir2 = './sp_resL11_2/'
SS = 0.0
cnt = 0
L = 11
Q = 0.0
for k in range(1,11):
    N1 = 2**(L-k)  # LL=10 so LL-kk=10,9,8,7,6,5,4,3,2,1
    N2 = 2**k #                    k=1,2,3,4,5,6,7,8,9,10 
    N2b = 2**(k-1)
    SS = SSb = 0.0
    for q in range(1,9):
      cnt = cnt +1
      fname = fdir+'spS_N1_'+str(N1)+'_N2_'+str(N2)+'_q_'+str(q)+'.npz' 
      myf = np.load(fname)
      SN = myf['mySN']
      SS = SS + SN
      if (k>1):
        fname2 = fdir2+'spS_N1_'+str(N1)+'_N2_'+str(N2b)+'_q_'+str(q)+'.npz' 
        myf2 = np.load(fname2)
        SNb = myf2['mySN']
        SSb = SSb + SNb
    Qk = SS/ (2**L) - SSb/ (2**(L-1))
    print('k= ',k,' Qk=',Qk)
    Q = Q+Qk
print('Q = ', Q)
# Q = -0.12598821219998504 
