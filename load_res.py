# load results
import numpy as np
N = 1024
N = 2048
if (N == 1024):
  fdir = './resN1024/'
  maxq = 4  # 4 partition, each has 256 qmc points
elif (N == 2048):
  fdir = './resN2048/'
  maxq = 8 
SS = 0.0
cnt = 0
for p in range(1,N+1):
   for q in range(1,maxq+1):
      cnt = cnt +1
      fname = fdir+'S_N'+str(N)+'_p_'+str(p)+'_q_'+str(q)+'.npz'
      myf = np.load(fname)
      SN = myf['mySN']
      SS = SS + SN

print('cnt= ', cnt, SS)
print('QN=',SS/(256*cnt))
# QN= -0.1567797096295199 (for N=1024)
# QN= -0.15678003292188647 (for N=2048)
