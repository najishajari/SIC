from __future__ import print_function
import numpy as np
#from matplotlib import pyplot as plt 
from scipy import signal
import h5py
from time import time
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import xml.etree.cElementTree as ET
from statsmodels.tsa.stattools import grangercausalitytests

import sys,os
#a line to add Codes path to be able to use mod_welch module!
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from mod_welch import welch

def delta_estimator_1(h,S_x):
    return (np.mean(np.multiply(h,S_x))-(np.mean(h)*np.mean(S_x)))/(np.mean(np.multiply(h,S_x)))
def delta_estimator_2(h,S_x):
    return (np.mean(np.multiply(h,S_x))-(np.mean(h)*np.mean(S_x)))/(np.sum(h)*np.max(S_x))
def delta_estimator_3(h,S_x):
    int_h_S_x=np.mean(np.multiply(h,S_x))
    int_h=np.mean(h)
    int_S_x=np.mean(S_x)
    return np.log(int_h_S_x)-np.log(np.mean(int_h))-np.log(int_S_x)
    
def win_sig(x,nperseg):
    """A function just to cut a multidimensional time series into pieces of specific length (nperseg) """
    #index set manipulation for generating the splitted version of the signals faster                        
    sig_len=x.shape[0]-nperseg+1
    idx_temp=np.indices((sig_len,nperseg))
    idx_temp=idx_temp[0]+idx_temp[1]
    
    return x[idx_temp].reshape((-1,nperseg))
def AR_fit(x,y,order):
    inp_len=x.shape[0]
    x_mat=win_sig(x,order)
    y_mat=y[order-1:]
    w=np.linalg.lstsq(x_mat,y_mat)[0]
    return w
    
def SIC_method(X,Y,order=100):   
    #low-passing to take LFP only
    
    h_for=AR_fit(X,Y,order)
    y_new=signal.fftconvolve(h_for,X)
    h_back=AR_fit(Y,X,order)
    x_new=signal.fftconvolve(h_back,Y)

    #Sx=welch(x_new,nperseg=1000)[1]
    #Sy=welch(y_new,nperseg=1000)[1]

#    Sy=welch(Y,nperseg=1000)[1]
#    Sx=welch(X,nperseg=1000)[1]
#
#    X_Y=delta_estimator_3(Sy/Sx,Sx)
#    Y_X=delta_estimator_3(Sx/Sy,Sy)
            
    #mask1=Sy!=0
    #mask2=Sx[mask1]!=0
    #X_Y=eval('delta_estimator_'+str(method_no))(Sy[mask1][mask2][1:-1]/Sx[mask1][mask2][1:-1],Sx[mask1][mask2][1:-1])
    #Y_X=eval('delta_estimator_'+str(method_no))(Sx[mask1][mask2][1:-1]/Sy[mask1][mask2][1:-1],Sy[mask1][mask2][1:-1])
    #X_Y=eval('delta_estimator_'+str(method_no))(Sy[mask1][mask2][1:-1]/Sx[mask1][mask2][1:-1],Sx[mask1][mask2][1:-1])
    #Y_X=eval('delta_estimator_'+str(method_no))(Sx[mask1][mask2][1:-1]/Sy[mask1][mask2][1:-1],Sy[mask1][mask2][1:-1])

    X_Y=np.var(y_new)/float(np.sum(h_for**2)*np.var(X))
    Y_X=np.var(x_new)/float(np.sum(h_back**2)*np.var(Y))
    
    return X_Y,Y_X




parser = ArgumentParser(description='--t: Number of process', formatter_class=RawTextHelpFormatter)
parser.add_argument('--t',required=True)
args = vars(parser.parse_args())
task_number=int(eval(args['t']))-1

data_path='/is/ei/naji/Real Data/Neuroscience'
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
    f.close()
with h5py.File(data_path+'/'+mice_name+'_'+session_name+'_CA3.h5','r') as f:
    X = f['/data'].value
    f.close()
print (Y.shape)
print (X.shape)
res_SIC=[]
sample_len=int(sampling_rate)*2
f=open(cluster_path+'/'+mice_name+'/'+session_name+'/output/'+str(task_number/32)+','+str(np.mod(task_number,32))+'.txt','w')    
#f_GC=open(cluster_path+'/'+mice_name+'/'+session_name+'/GC_output/'+str(task_number/32)+','+str(np.mod(task_number,32))+'.txt','w')    
for i in np.arange(300):
    row=task_number/32
    column=np.mod(task_number,32)
    #calculating granger causality for time intervals
    X_window=X[i*sample_len:(i+1)*sample_len,row]
    Y_window=Y[i*sample_len:(i+1)*sample_len,column]

    #Granger causality checking step:
    #GC_XY=grangercausalitytests(np.vstack((Y_window,X_window)).T,maxlag=10,verbose=False)[10][0]['ssr_ftest']
    #GC_YX=grangercausalitytests(np.vstack((X_window,Y_window)).T,maxlag=10,verbose=False)[10][0]['ssr_ftest']
    #print (GC_XY[1],';',GC_YX[1],file=f_GC)
    #print (np.asarray(grangercausalitytests(np.vstack((Y_window,X_window)).T,maxlag=10,verbose=False))[:][0])

    res_SIC=SIC_method(X_window,Y_window)
    print (res_SIC[0],';',res_SIC[1],file=f)
    del X_window
    del Y_window
f.close()
#f_GC.close()