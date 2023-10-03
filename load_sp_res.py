# load results
# usage: python load_sp_res.py L
import numpy as np
import sys
L = 12
L = int(sys.argv[1])
if (L==9):
   fdir = './sp_resL9/'
   fdir2 = './sp_resL9_2/'
   qmax = 2
elif (L==10):
   fdir = './sp_resL10/'
   fdir2 = './sp_resL10_2/'
   qmax = 4
elif (L==11):
   fdir = './sp_resL11/'
   fdir2 = './sp_resL11_2/'
   qmax = 8
elif (L==12):
   fdir = './sp_resL12/'
   fdir2 = './sp_resL12_2/'
   qmax = 16 
elif (L==13):
   fdir = './sp_resL13/'
   fdir2 = './sp_resL13_2/'
   qmax = 32 
elif (L==14):
   fdir = './sp_resL14/'
   fdir2 = './sp_resL14_2/'
   qmax = 64 
elif (L==15):
   fdir = './sp_resL15/'
   fdir2 = './sp_resL15_2/'
   qmax = 128

#Qfull= -0.1567797096295199  # for full grid N=1024
Qfull = -0.15678003292188647 # (for N=2048)
SS = 0.0
cnt = 0
Q = 0.0
print('L = ',L)
for k in range(0,L+1):
    N1 = 2**(L-k)  # LL=10 so LL-kk=9,8,7,6,5,4,3,2,1
    N2 = 2**k #                    k=1,2,3,4,5,6,7,8,9 
    N2b = 2**(k-1)
    SS = SSb = 0.0
    for q in range(1,qmax+1):
      cnt = cnt +1
      fname = fdir+'est1_N1_'+str(N1)+'_N2_'+str(N2)+'_q_'+str(q)+'.npz' 
      myf = np.load(fname)
      SN = myf['myS1']
      SS = SS + SN
      if (k>0):
        fname2 = fdir2+'est2_N1_'+str(N1)+'_N2_'+str(N2b)+'_q_'+str(q)+'.npz' 
        myf2 = np.load(fname2)
        SNb = myf2['mySN']
        SSb = SSb + SNb
    Qk = SS/ (2**L) - SSb/ (2**(L-1))
    print('k= ',k,' Qk=',Qk, 'Qk_1=', SS/(2**L))
    Q = Q+Qk

print('Q = ', Q, 'err =', np.abs(Q-Qfull))
# for the table in the paper
err= np.abs(Q-Qfull)
upper = 2**(2-L)*L
print(f'{L:2d} {err:9.4e} {upper:9.4e}')
