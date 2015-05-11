import numpy as np
from matplotlib import pyplot as plt 
from scipy import signal
import h5py
from time import time
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import xml.etree.cElementTree as ET
from statsmodels.tsa.stattools import grangercausalitytests
data_path='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Real Data/Neuroscience'
cluster_path='/agbs/cluster/naji/Linear Filters/real_data'
manuscript_path='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Manuscripts/'

mice_name='vvp01'
session_name='2006-4-9_18-43-47'
#extracting the lfp sampling rate:
tree=ET.ElementTree(file=data_path+'/'+mice_name+'/'+session_name+'/'+session_name+'.xml')
atype=tree.findall('fieldPotentials')
for btype in atype[0].findall('lfpSamplingRate'):
    sampling_rate=btype.text

with h5py.File(data_path+'/'+mice_name+'_'+session_name+'_CA1.h5','r') as f:
    Y = f['/data'].value
with h5py.File(data_path+'/'+mice_name+'_'+session_name+'_CA3.h5','r') as f:
    X = f['/data'].value

for i in range(32):
    plt.plot(np.arange(100*1250,112*1250)/1250.,X[100*1250:112*1250,i]+i*10000)
plt.xlabel(r'${\rm Time(s)}$')
plt.ylabel(r'${\rm LFP}$')    
plt.gca().axes.yaxis.set_ticklabels([])
plt.savefig(manuscript_path+'Linear_Filters/Figures/LFP[104,134]_CA3_vvp01_2006-4-9_18-43-87.eps',transparent=True)
plt.figure()
for i in range(32):
    plt.plot(np.arange(100*1250,112*1250)/1250.,Y[100*1250:112*1250,i]+i*10000)
plt.xlabel(r'${\rm Time(s)}$')
plt.ylabel(r'${\rm LFP}$')    
plt.gca().axes.yaxis.set_ticklabels([])
plt.savefig(manuscript_path+'Linear_Filters/Figures/LFP[104,134]_CA1_vvp01_2006-4-9_18-43-87.eps',transparent=True)

#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0., prop={'size':20})    
plt.show()
#for i in range(32):
#    ax[0].plot(np.arange(105*1250,108*1250)/1250.,X[105*1250:108*1250,i]+i*10000)
#    ax[1].plot(np.arange(105*1250,108*1250)/1250.,Y[105*1250:108*1250,i]+i*10000)
#ax[0].axvspan(106,107,facecolor='r', alpha=0.5)
#ax[1].axvspan(106,107,facecolor='r', alpha=0.5)
#font = {'family' : 'normal',
#        'weight' : 'normal',
#        'size'   : 20}
#plt.rc('font', **font)

#ax[0].set_xlabel(r'${\rm Time(s)}$')
#ax[0].set_ylabel(r'${\rm CA3\ \ LFP}$')
#ax[1].set_ylabel(r'${\rm CA1\ \ LFP}$')
#plt.rc('text',usetex=True)
#plt.savefig(manuscript_path+'Linear_Filters/Figures/LFP[100,112]_CA1_and_CA3_vvp01_2006-4-9_18-43-87.png',dpi=200,transparent=True)
#fig.set_size_inches(30.40,20.64)
#ax[0].axes.yaxis.set_ticklabels([])






