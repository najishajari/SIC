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

mice_name='vvp01'
session_name='2006-4-9_17-29-30'
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
    plt.plot(np.arange(180*1250)/1250.,X[:180*1250,i]+i*10000)
plt.xlabel(r'${\rm Time(s)}$')
plt.ylabel(r'${\rm LFP}$')    
plt.gca().axes.yaxis.set_ticklabels([])
plt.figure()
for i in range(32):
    plt.plot(np.arange(180*1250)/1250.,Y[:180*1250,i]+i*10000)
plt.xlabel(r'${\rm Time(s)}$')
plt.ylabel(r'${\rm LFP}$')    
plt.gca().axes.yaxis.set_ticklabels([])

#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0., prop={'size':20})    
plt.show()