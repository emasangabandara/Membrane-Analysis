from __future__ import division
import numpy as np
import pylab as pl
from matplotlib import pyplot as plt
import seaborn as sns

sns.set(style="ticks",context="poster", font_scale=1.2, rc={"lines.linewidth": 4})

sim_time = 10.0 #Simulation time in microseconds

#Moving Average Module for Sommething the Timeseries
def runningMeanFast(x, N):
    return np.convolve(x.reshape(-1,), np.ones((N,))/N)[(N-1):]

#Custom plotting module for time series    
def plot_timeseries(ts_rep,color,label,linestyle) :
    mean_ts                   = np.mean(ts_rep,axis=0)  
    std_ts                    = np.std(ts_rep,axis=0)

    #Creating dummy arrey for the time axis
    arr_length                  = np.shape(mean_ts)
    time                        = np.arange(0,sim_time,sim_time/arr_length[0])

    plt.plot(time,mean_ts,linestyle=linestyle,color=color,alpha=1.0, linewidth=2,label=label) 
    plt.fill_between(time ,mean_ts-std_ts,mean_ts+std_ts,alpha=0.2,facecolor=color,edgecolor='none')
 
WORKDIR  = '.'
lpd11_u      = np.loadtxt(r'%s/upper_Lpd1_Lpd1.dat'%WRKDIR)
lpd11_l      = np.loadtxt(r'%s/lower_Lpd1_Lpd1.dat'%WRKDIR)
    
lpd11_rep = runningMeanFast(lpd11_u,20)
lpd11_rep = runningMeanFast(lpd11_l,20)

lpd12_u      = np.loadtxt(r'%s/upper_Lpd1_Lpd2.dat'%WRKDIR)
lpd12_l      = np.loadtxt(r'%s/lower_Lpd1_Lpd2.dat'%WRKDIR)
    
lpd12_rep = runningMeanFast(lpd12_u,20)
lpd12_rep = runningMeanFast(lpd12_l,20)

lpd13_u      = np.loadtxt(r'%s/upper_Lpd1_Lpd3.dat'%WRKDIR)
lpd13_l      = np.loadtxt(r'%s/lower_Lpd1_Lpd3.dat'%WRKDIR)
    
lpd13_rep = runningMeanFast(lpd13_u,20)
lpd13_rep = runningMeanFast(lpd13_l,20)

lpd22_u      = np.loadtxt(r'%s/upper_Lpd2_Lpd2.dat'%WRKDIR)
lpd22_l      = np.loadtxt(r'%s/lower_Lpd2_Lpd2.dat'%WRKDIR)
    
lpd22_rep = runningMeanFast(lpd22_u,20)
lpd22_rep = runningMeanFast(lpd22_l,20)        

lpd23_u      = np.loadtxt(r'%s/upper_Lpd2_Lpd3.dat'%WRKDIR)
lpd23_l      = np.loadtxt(r'%s/lower_Lpd2_Lpd3.dat'%WRKDIR)
    
lpd23_rep = runningMeanFast(lpd23_u,20)
lpd23_rep = runningMeanFast(lpd23_l,20)

lpd33_u      = np.loadtxt(r'%s/upper_Lpd3_Lpd3.dat'%WRKDIR)
lpd33_l      = np.loadtxt(r'%s/lower_Lpd3_Lpd3.dat'%WRKDIR)
    
lpd33_rep = runningMeanFast(lpd33_u,20)
lpd33_rep = runningMeanFast(lpd33_l,20)


prtlpd1_u      = np.loadtxt(r'%s/upper_Prot_Lpd1.dat'%WRKDIR)
prtlpd1_l      = np.loadtxt(r'%s/lower_Prot_Lpd1.dat'%WRKDIR)
    
prtlpd1_rep = runningMeanFast(prtlpd1_u,20)
prtlpd1_rep = runningMeanFast(prtlpd1_l,20)

prtlpd2_u      = np.loadtxt(r'%s/upper_Prot_Lpd2.dat'%WRKDIR)
prtlpd2_l      = np.loadtxt(r'%s/lower_Prot_Lpd2.dat'%WRKDIR)
    
prtlpd2_rep = runningMeanFast(prtlpd2_u,20)
prtlpd2_rep = runningMeanFast(prtlpd2_l,20)

prtlpd3_u      = np.loadtxt(r'%s/upper_Prot_Lpd3.dat'%WRKDIR)
prtlpd3_l      = np.loadtxt(r'%s/lower_Prot_Lpd3.dat'%WRKDIR)
    
prtlpd3_rep = runningMeanFast(prtlpd3_u,20)
prtlpd3_rep = runningMeanFast(prtlpd3_l,20)

prtprt_u      = np.loadtxt(r'%s/upper_Prot_Prot.dat'%WRKDIR)
prtprt_l      = np.loadtxt(r'%s/lower_Prot_Prot.dat'%WRKDIR)
    
prtprt_rep = runningMeanFast(prtprt_u,20)
prtprt_rep = runningMeanFast(prtprt_l,20) 
     
tot_lpd_contacts = np.sum((lpd11_rep,lpd22_rep,lpd33_rep,lpd12_rep,lpd13_rep,lpd23_rep), axis=0)
tot_prot_contacts = np.sum((prtlpd1_rep,prtlpd2_rep,prtlpd3_rep), axis=0)
plt.subplot(211)  
plot_timeseries(100*np.asarray(lpd11_rep)/tot_lpd_contacts,'blue','DPPC-DPPC','-')
plot_timeseries(100*np.asarray(lpd22_rep)/tot_lpd_contacts,'red','DIPC-DIPC','-')
plot_timeseries(100*np.asarray(lpd33_rep)/tot_lpd_contacts,'green','CHOL-CHOL','-')
plot_timeseries(100*np.asarray(lpd12_rep)/tot_lpd_contacts,'magenta','DPPC-DIPC','-')
plot_timeseries(100*np.asarray(lpd13_rep)/tot_lpd_contacts,'cyan','DPPC-CHOL','-')
plot_timeseries(100*np.asarray(lpd23_rep)/tot_lpd_contacts,'brown','DIPC-CHOL','-') 
plt.ylim(0,35)
plt.xlim(0,16.8)
plt.ylabel('# of contacts')
plt.xlabel('time (${/mu}$s)')
plt.legend(loc='best')

plt.subplot(212)      
plot_timeseries(100*np.asarray(prtlpd1_rep)/tot_prot_contacts,'blue','Prot-DPPC','-')
plot_timeseries(100*np.asarray(prtlpd2_rep)/tot_prot_contacts,'red','Prot-DIPC','-')
plot_timeseries(100*np.asarray(prtlpd3_rep)/tot_prot_contacts,'green','Prot-CHOL','-')
#plot_timeseries(prtprt_rep,'grey','Prot-Prot','-')   

plt.ylim(0,100)
plt.xlim(0,16.8)
plt.ylabel('% contacts')
plt.xlabel('time (${/mu}$s)')
plt.legend(loc='best')

plt.tight_layout()
plt.savefig(r'plot_lipd_prot_percent_contact_ts_%s.svg')   
    






