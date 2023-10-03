# load results
import numpy as np
import sys
L = 9 
L = int(sys.argv[1])
if (L==9):
   fdir = './sp_resL9/'
   fdir2 = './sp_resL9_2/'
elif (L==10):
   fdir = './sp_resL10/'
   fdir2 = './sp_resL10_2/'
elif (L==11):
   fdir = './sp_resL11/'
   fdir2 = './sp_resL11_2/'
elif (L==12):
   fdir = './sp_resL12/'
   fdir2 = './sp_resL12_2/'
elif (L==13):
   fdir = './sp_resL13/'
   fdir2 = './sp_resL13_2/'
elif (L==14):
   fdir = './sp_resL14/'
   fdir2 = './sp_resL14_2/'
elif (L==15):
   fdir = './sp_resL15/'
   fdir2 = './sp_resL15_2/'

for k in range(0,L+1):
    N1 = 2**(L-k)  # LL=10 so LL-kk=9,8,7,6,5,4,3,2,1
    N2 = 2**k #                    k=1,2,3,4,5,6,7,8,9 
    N2b = 2**(k-1)
    fname = 'est1_N1_'+str(N1)+'_N2_'+str(N2)+'_q_*.npz' 
    print('mv ',fname,' ',fdir)
    if (k>0):
        fname2 = 'est2_N1_'+str(N1)+'_N2_'+str(N2b)+'_q_*.npz' 
        print('mv ',fname2,' ',fdir2)
