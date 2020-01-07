from __future__ import division
import numpy as np
import pylab as pl
from matplotlib import pyplot as plt

#Moving Average Module
def runningMeanFast(x, N):
    return np.convolve(x.reshape(-1,), np.ones((N,))/N)[(N-1):]

def calcBinaryMixEntropy(Lpd1_Lpd1,Lpd2_Lpd2,Lpd1_Lpd2):
    SL = Lpd1_Lpd1 + Lpd2_Lpd2
    DL = Lpd1_Lpd2
    
    X_SL = np.true_divide(SL, SL + DL) #Similar Lipid
    X_DL = np.true_divide(DL, SL + DL) #Dissimilar Lipid
    
    binary_mix_entropy = -((X_SL * np.log2(X_SL)) +( X_DL * np.log2(X_DL)))

    return binary_mix_entropy
    
WRKDIR       = "."

lpd11_u      = np.loadtxt(r'%s/upper_Lpd1_Lpd1.R1.dat'%WRKDIR)
lpd11_l      = np.loadtxt(r'%s/lower_Lpd1_Lpd1.R1.dat'%WRKDIR)

lpd22_u      = np.loadtxt(r'%s/upper_Lpd2_Lpd2.R1.dat'%WRKDIR)
lpd22_l      = np.loadtxt(r'%s/lower_Lpd2_Lpd2.R1.dat'%WRKDIR)

lpd12_u      = np.loadtxt(r'%s/upper_Lpd1_Lpd2.R1.dat'%WRKDIR)
lpd12_l      = np.loadtxt(r'%s/lower_Lpd1_Lpd2.R1.dat'%WRKDIR)

binary_mix_entropy_u = calcBinaryMixEntropy(lpd11_u,lpd22_u,lpd12_u)
binary_mix_entropy_l = calcBinaryMixEntropy(lpd11_l,lpd22_l,lpd12_l)
            

arr_length = np.shape(binary_mix_entropy_u)
time       = np.arange(0,30.0,30.0/arr_length[0])

plt.plot(time,binary_mix_entropy_u,'-',color='red',alpha=0.75, linewidth=1,label='Upper-leaflet') 
plt.plot(time,binary_mix_entropy_l,'-',color='blue',alpha=0.75, linewidth=1,label='Lower-leaflet') 

plt.ylim(0.0,1.0)
plt.xlim(0,30)
plt.ylabel('$<S_{mix}>$')
plt.xlabel('time (${\mu}$s)')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Binary-mixing-entropy.png')
#plt.show()