from mod_welch import welch
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt
def SIC(X,Y):
    Sx=welch(X,nperseg=500)[1]
    Sy=welch(Y,nperseg=500)[1]
    Delta_X_Y=np.sum(Sy)/(np.sum(Sx)*np.sum(Sy/Sx))
    Delta_Y_X=np.sum(Sx)/(np.sum(Sy)*np.sum(Sx/Sy))
    if Delta_X_Y>Delta_Y_X:
        print "X_t causes Y_t"
    else:
        print "Y_t causes X_t"

# Illustration of example in the paper
FO=10
x_size=10000
AR_amp_upperbound=0.1
A=np.random.randn(FO)*AR_amp_upperbound
#input as random noise
X=np.random.randn(x_size)
#setting up an IIR filter with white noise as inpuit and coeffcieints :A
Y=signal.lfilter(A,[1],X)
SIC(X,Y)


