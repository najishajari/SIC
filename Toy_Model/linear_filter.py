"""script to generate linear dynamical system"""
from __future__ import print_function
import sys,os
#a line to add Codes path to be able to use mod_welch module!
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from mod_welch import welch

import numpy as np
from scipy import signal
#from scipy.stats import truncnorm
from time import time
from matplotlib.font_manager import FontProperties
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from scipy import fftpack
def cov(a,b):
    return np.mean(np.multiply((a-np.mean(a)),(b-np.mean(b))))
def delta_estimator_1(h,S_x):
    return (np.mean(np.multiply(h,S_x))-(np.mean(h)*np.mean(S_x)))/(np.mean(np.multiply(h,S_x)))
def delta_estimator_2(h,S_x):
    return (np.mean(np.multiply(h,S_x))-(np.mean(h)*np.mean(S_x)))/(np.sum(h)*np.max(S_x))

def SDR_estimator(h,S_x):
    S_x=np.ma.masked_invalid(S_x)    
    h=np.ma.masked_invalid(h)    
    int_h_S_x=np.mean(np.multiply(h,S_x))
    int_h=np.mean(h)
    int_S_x=np.mean(S_x)
    return float(int_h_S_x)/(float(np.mean(int_h))*float(int_S_x))

def delta_estimator_4(h,S_x):
    """SDR_estimator() with log differnece rather than ratio"""
    
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


def linear_filter_SIC(FO=5, BO=5,ksim_size=1000, x_size=10000, AR_amp_upperbound=0.1, noise_amp=0.1, in_noise=False, out_noise=False, denoise=False, win_length=500):
    """ The main function for generating a toy model consisting of two linear filters. The input of the first filter is white noise and the output
    of it is feeded to the second filter as X_t (the cause). The output of the second filter in turn is considered to be Y_t. The coefficients of the filters
    are also i.i.d gaussian variables. The function finally returns the SDRs in both directions for a number of trials specified by ksim_size.
    The linear filters are defined according scipy IIR filters:  
                                
                                scipy.signal.lfilter(b, a, x, axis=-1, zi=None)     (**)
                                
    (See "http://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html" for more info on IIR filters)
        
    input:
        FO:                     int: Forward order of IIR filter(length of b in (**))
        BO:                     int: Backward order of IIR filter(length of a in (**))
        ksim_size:              int: the number of trials of these toy model
        x_size:                 int: the length of the input white noise to the first filter
        AR_amp_upperbound:      float: a coefficient to control the largness of the feedback coefficients (to maitain BIBO stability)
        noise_amp:              float: the amplitude of noise (standard deviation of the white additive noise)
        in_noise:               boolean: whether there is input additive noise
        out_noise:              boolean: whether there is output additive noise
        denoise:                boolean: whether we denoise the observations by fitting an autoregressive filter in both directions (presummably the real model!)
    
    output:
        Dir                    numpy ndarray.float64: stored values of SDR for both directions for different trials and also the variances of X_t and Y_t
    """
    #initializing Dir matrix
    Dir=np.empty((ksim_size,2,2))
    
    #looping for different trials
    for ksim in xrange(ksim_size):
        #setting the initial coefficients of filter assuming that the first element of coefficients of filters in both directions is 1 (for symmetry)
        A_=np.random.randn(FO)*AR_amp_upperbound
        A=np.hstack((np.ones(1),-A_))
        B_=np.random.randn(FO)*AR_amp_upperbound
        B=np.hstack((np.ones(1),-B_))
        
        E_=np.random.randn(BO-1)*AR_amp_upperbound
        F_=np.random.randn(BO-1)*AR_amp_upperbound
        E=np.hstack((np.ones(1),-E_))
        F=np.hstack((np.ones(1),-F_))

        #re-setting the coefficients in case either of the two filters are not BIBO stable        
        while np.any(np.abs(np.roots(E))>=0.99) or np.any(np.abs(np.roots(A))>=0.99):
            A_=np.random.randn(FO)*AR_amp_upperbound
            A=np.hstack((np.ones(1),-A_))
            E_=np.random.randn(BO-1)*AR_amp_upperbound
            E=np.hstack((np.ones(1),-E_))
        while np.any(np.abs(np.roots(F))>0.99) or np.any(np.abs(np.roots(B))>=0.99):
            B_=np.random.randn(FO)*AR_amp_upperbound
            B=np.hstack((np.ones(1),-B_))
            F_=np.random.randn(BO-1)*AR_amp_upperbound
            F=np.hstack((np.ones(1),-F_))
        
        # printing the final selected coefficients
        print ("A:",A)
        print ("B:",B)
        
        # generating the cause, X_t (here x) and the effect, Y_t (here y)
        x=signal.lfilter(A,E,np.random.randn(x_size))
        y=signal.lfilter(B,F,x)
        
        # storing total variances of cause and effect time series
        Dir[ksim,0,0]=np.var(x)
        Dir[ksim,1,1]=np.var(y)
        
        # adding noise
        if out_noise:
            y+=np.random.randn(x_size)*noise_amp
        if in_noise:
            x+=np.random.randn(x_size)*noise_amp
        
        # denoising if necessarily
        if denoise:
            order=20 #the AR order
            h_for=AR_fit(x,y,order)
            y_new=signal.fftconvolve(h_for,x)
            h_back=AR_fit(y,x,order)
            x_new=signal.fftconvolve(h_back,y)

            Dir[ksim,0,1]=np.var(y_new)/float(np.sum(h_for**2)*np.var(x))
            Dir[ksim,1,0]=np.var(x_new)/float(np.sum(h_back**2)*np.var(y))

        else:
            Sy=welch(y,nperseg=win_length)[1]
            Sx=welch(x,nperseg=win_length)[1]            
            Dir[ksim,0,1]=SDR_estimator(Sy/Sx,Sx)    
            Dir[ksim,1,0]=SDR_estimator(Sx/Sy,Sy)
            
    return Dir
def linear_filter_decreasing_impulse(FO=5,BO=5,ksim_size=1000,x_size=10000,AR_amp_upperbound=0.1,noise_amp=0.1,in_noise=False,out_noise=False,denoise=False):
    """linear_filter_SIC() rewritten where the coefficients of IR functions are picked in a way that their absolute value decreases by the increase of time. For the definitions
    of arguments or the description of the logic of the code, see the body of linear_filter_SIC()."""
    Dir=np.empty((ksim_size,2,2))
    for ksim in xrange(ksim_size):
        A_=np.random.randn(FO)*AR_amp_upperbound
        A_idx=abs(A_).argsort()[::-1]
        A=A_[A_idx]
        B_=np.random.randn(FO)*AR_amp_upperbound
        B_idx=abs(B_).argsort()[::-1]
        B=B_[B_idx]
        
        E_=np.random.randn(BO-1)*AR_amp_upperbound
        F_=np.random.randn(BO-1)*AR_amp_upperbound
        E_idx=abs(E_).argsort()[::-1]
        F_idx=abs(F_).argsort()[::-1]
        E=np.hstack((np.ones(1),E_[E_idx]))
        F=np.hstack((np.ones(1),F_[F_idx]))
        while np.any(np.abs(np.roots(E))>=0.98) or np.any(np.abs(np.roots(A))>=0.98):
            A_=np.random.randn(FO)*AR_amp_upperbound
            A_idx=abs(A_).argsort()[::-1]
            A=A_[A_idx]
            E_=np.random.randn(BO-1)*AR_amp_upperbound
            E_idx=abs(E_).argsort()[::-1]
            E=np.hstack((np.ones(1),E_[E_idx]))
        while np.any(np.abs(np.roots(F))>0.98) or np.any(np.abs(np.roots(B))>=0.98):
            B_=np.random.randn(FO)*AR_amp_upperbound
            B_idx=abs(B_).argsort()[::-1]
            B=B_[B_idx]
            F_=np.random.randn(BO-1)*AR_amp_upperbound
            F_idx=abs(F_).argsort()[::-1]
            F=np.hstack((np.ones(1),F_[F_idx]))
        print ("A:",A)
        print ("B:",B)
                
        x=signal.lfilter(A,E,np.random.randn(x_size))
        y=signal.lfilter(B,F,x)
        Dir[ksim,0,0]=np.var(x)
        Dir[ksim,1,1]=np.var(y)
        
        if out_noise:
            y+=np.random.randn(x_size)*noise_amp
        if in_noise:
            x+=np.random.randn(x_size)*noise_amp
        
        if denoise:
            order=100
            
            h_for=AR_fit(x,y,order)
            y_new=signal.fftconvolve(h_for,x)
            h_back=AR_fit(y,x,order)
            x_new=signal.fftconvolve(h_back,y)

            Dir[ksim,0,1]=np.var(y_new)/float(np.sum(h_for**2)*np.var(x))
            Dir[ksim,1,0]=np.var(x_new)/float(np.sum(h_back**2)*np.var(y))

        else:
            Sy=welch(y,nperseg=1000)[1]
            Sx=welch(x,nperseg=1000)[1]
                        
            Dir[ksim,0,1]=SDR_estimator(Sy/Sx,Sx)
    
            Dir[ksim,1,0]=SDR_estimator(Sx/Sy,Sy)
    return Dir
    
def linear_filter_CoM(FO=10,BO=10,ksim_size=1000,x_size=10000,AR_amp_upperbound=0.2,noise_amp=3.1):
    """
    !!! HAS NOT BEEN UPDATED !!!
    A function to test the concentration of measure for the proved theorem in the paper
    """
    Dir=np.empty((ksim_size,2,2))
    for ksim in xrange(ksim_size):
        A=np.random.randn(FO)
        B=np.random.randn(FO)
        E_=np.random.randn(BO-1)*AR_amp_upperbound
        F_=np.random.randn(BO-1)*AR_amp_upperbound
        E=np.hstack((np.ones(1),E_))
        F=np.hstack((np.ones(1),F_))
        while np.any(np.abs(np.roots(E))>=1):
            E_=np.random.randn(BO-1)*AR_amp_upperbound
            E=np.hstack((np.ones(1),E_))
            print (FO)
        while np.any(np.abs(np.roots(F))>=1):
            F_=np.random.randn(BO-1)*AR_amp_upperbound
            F=np.hstack((np.ones(1),F_))
                        
        x=signal.lfilter(A,E,np.random.randn(x_size))
        y=signal.lfilter(B,F,x)+np.random.randn(x_size)*noise_amp        
        x=x+np.random.randn(x_size)*noise_amp
        Sy=signal.welch(y,nperseg=256)[1]
        Sx=signal.welch(x,nperseg=256)[1]
        
        Dir[ksim,0,0]=np.sum(Sx)
        Dir[ksim,0,1]=cov(Sy[1:-1]/Sx[1:-1],Sx[1:-1])/(np.max(Sy[1:-1]/Sx[1:-1])*np.max(Sx[1:-1]))

        Dir[ksim,1,0]=cov(Sx[1:-1]/Sy[1:-1],Sy[1:-1])/(np.max(Sx[1:-1]/Sy[1:-1])*np.max(Sy[1:-1]))
        Dir[ksim,1,1]=np.sum(Sy)
    return Dir

def SIC_performance_dim(trials=12,ksim_size=1000,xsize=10000,order_increment=5,param=None):
    """ NOT IN USE!!!
    
    Running the filter on cluster and returning the performance
        arguments:
            -param =            0 if feedforward order is 0
                                1 if feedback order is 0
                                2 if none are zero
    """
    Dirs=np.empty((trials,ksim_size,2,2))
    for i in xrange(trials):
        if param==0:
            Dirs[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=1,BO=order_increment*(i+1)-1,AR_amp_upperbound=0.1,noise_amp=.0)
            output_name='order=[2,20],x_size=10000,ksim _size=1000,noise=null,fw_dim=0'
        elif param==1:
            Dirs[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=order_increment*(i+1),BO=order_increment*(i+1)-1,AR_amp_upperbound=0,noise_amp=.0)
            output_name='order=[2,20],x_size=10000,ksim_size=1000,noise=null,bw_dim=0'
        else:
            Dirs[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=order_increment*(i+1),BO=order_increment*(i+1)-1,AR_amp_upperbound=0.1,noise_amp=.0)
            output_name='order=[2,20],x_size=10000,ksim_size=1000,noise=null,fw_dim=bw_dim'
    f=open('/is/ei/naji/Dropbox/Winter Semester 2014/Master Thesis/Programming/Multivariate IGCI/Linear Filters/Dimension Change/'+output_name,'w')
    print (Dirs[:,0,0,1].shape)
    print (np.mean(Dirs[:,:,0,1]>Dirs[:,:,1,0],axis=1),file=f)
    f.close()





