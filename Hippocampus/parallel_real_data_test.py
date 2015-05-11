from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt 
from scipy import signal
import h5py
from time import time
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
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


def SIC_method(X,Y,method_no,order=200):   
    #low-passing to take LFP only

    h_for=AR_fit(X,Y,order)
    y_new=signal.fftconvolve(h_for,X)
    h_back=AR_fit(Y,X,order)
    x_new=signal.fftconvolve(h_back,Y)
    
    Sx=welch(x_new,nperseg=500)[1]
    Sy=welch(y_new,nperseg=500)[1]

    #Sy=welch(Y,nperseg=500)[1]
    #Sx=welch(X,nperseg=500)[1]
    
    mask1=Sy!=0
    mask2=Sx[mask1]!=0
    #plt.plot(Sx)
    #plt.show()
    X_Y=eval('delta_estimator_'+str(method_no))(Sy[mask1][mask2][1:-1]/Sx[mask1][mask2][1:-1],Sx[mask1][mask2][1:-1])
    Y_X=eval('delta_estimator_'+str(method_no))(Sx[mask1][mask2][1:-1]/Sy[mask1][mask2][1:-1],Sy[mask1][mask2][1:-1])
    #return abs(X_Y)<abs(Y_X)
    return X_Y,Y_X

parser = ArgumentParser(description='--t: Number of process', formatter_class=RawTextHelpFormatter)
parser.add_argument('--t',required=True)
args = vars(parser.parse_args())
task_number=int(eval(args['t']))-1
x_len=2000

data_path='/is/ei/naji/Real Data/Neuroscience'
cluster_path='/agbs/cluster/naji/Linear Filters/real_data'

mice_name='ec016'

with h5py.File(data_path+'/'+mice_name+'_DG.h5','r') as f:
    X = f['/data'].value
with h5py.File(data_path+'/'+mice_name+'_CA3.h5','r') as f:
    Y = f['/data'].value
print (X.shape)
print (Y.shape)


res_1=[]
res_2=[]
res_3=[]
t=time()
CA1_size_reduction=0
CA3_size_reduction=0
idx=[]
divident=task_number/X.shape[1]
rem=task_number-divident*X.shape[1]

res_1=SIC_method(X[:x_len,rem],Y[:x_len,divident],1)
res_2=SIC_method(X[:x_len,rem],Y[:x_len,divident],2)
res_3=SIC_method(X[:x_len,rem],Y[:x_len,divident],3)
f=open(cluster_path+'/'+mice_name+'/output/'+str(divident)+','+str(rem)+'.txt','w')    
print (res_1[0],';',res_1[1],';',res_2[0],';',res_2[1],';',res_3[0],';',res_3[1],file=f)
f.close()
print (time()-t)
##
