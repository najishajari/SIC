# -*- coding: utf-8 -*-
from __future__ import print_function
import wave
import sys
import numpy as np
from scipy.io import wavfile
from matplotlib import pyplot as plt
from scipy.fftpack import fft, ifft,fftfreq
from scipy.signal import fftconvolve,convolve
from scipy.ndimage import convolve1d
from os import listdir
from utility import pcm2float,float2pcm
import prettyplotlib as ppl
from prettyplotlib import brewer2mpl
from scipy import signal 
from multidim_welch import ndim_welch
from time import time 
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

def delta_estimator_1(h,S_x):
    return (np.mean(np.multiply(h,S_x))-(np.mean(h)*np.mean(S_x)))/(np.mean(np.multiply(h,S_x)))
def delta_estimator_2(h,S_x):
    S_x=np.ma.masked_invalid(S_x)    
    return (np.mean(np.multiply(h,S_x))-(np.mean(h)*np.mean(S_x)))/(np.sum(h)*np.max(S_x))
def delta_estimator_3(h,S_x):
    #S_x=np.ma.masked_invalid(S_x)    
    #h=np.ma.masked_invalid(h)    
    int_h_S_x=np.mean(np.multiply(h,S_x))
    int_h=np.mean(h)
    int_S_x=np.mean(S_x)
    return np.log(int_h_S_x)-np.log(np.mean(int_h))-np.log(int_S_x)


main_path='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/'
#main_path='/Users/Naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/'

#input_dir=main_path+'Sounds/Winter-MPH/Original/'
#output_dir=main_path+'Sounds/Winter-MPH/Trials/2/'

song_name='Winter'

input_dir=main_path+'Sounds/'+song_name+'/Original/2/'
output_dir=main_path+'Sounds/'+song_name+'/Experiments/Hall/3/'


input_file_names=(listdir(input_dir+'Segments/'))
output_file_names=(listdir(output_dir+'Segments/'))
def estim_diff(percent=256):
    sound_counter=0
    res=np.empty(len(input_file_names))
    for i in range(res.shape[0]):
        input_rate,input_sig=wavfile.read(input_dir+'Segments/'+input_file_names[i])
        output_rate,output_sig=wavfile.read(output_dir+'Segments/'+output_file_names[i])
    
        input_sig=pcm2float(input_sig,'float32')
        output_sig=pcm2float(output_sig,'float32')
        
        min_size=np.min((input_sig[:,0].shape[0],output_sig[:,0].shape[0]))
        #print min_size,min_size*percent
        #S_inp=np.absolute(fft(input_sig[:min_size,0]-np.mean(input_sig[:min_size,0])))
        #S_out=np.absolute(fft(output_sig[:min_size,0]-np.mean(output_sig[:min_size,0])))
    
        t=time()
        nperseg=int(min_size*percent)-np.mod(int(min_size*percent),10)
        real_perc=float(float(nperseg)/int(min_size*percent))
        S_inp=signal.welch(input_sig[:min_size,0],nperseg=nperseg)[1]    
        S_out=signal.welch(output_sig[:min_size,0],nperseg=nperseg)[1]    
        #S_inp=ndim_welch(input_sig[:min_size,0][None,...],nperseg=int(min_size*percent))[1]    
        #S_out=ndim_welch(output_sig[:min_size,0][None,...],nperseg=int(min_size*percent))[1]    
        #print time()-t
        #print S_inp_1,S_inp_2
        res[sound_counter]=delta_estimator_3(S_out/S_inp,S_inp)-delta_estimator_3(S_inp/S_out,S_out)
        #out=float2pcm(output_sig,'int16')
        sound_counter+=1
    return real_perc,int(min_size*percent),res
#matplot_figure=plt.figure()
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}
plt.rc('font', **font)
plt.rc('text', usetex=True)
#plt.imshow(res,aspect='auto',extent=[0,len(IR_file_names),0,len(sound_file_names)],interpolation='none')
#plt.title(r'${\rm Lacrimosa-Mozart}$')
#plt.colorbar()
#plt.plot(res)
#plt.errorbar(np.arange(len(res)),res,color='b',ecolor='b',fmt='--o',label=r'${\rm Max Planck Hall\}$')
#plt.savefig(output_dir+song_name+'.eps',transparent=True)
#py.plot_mpl(matplot_figure)
#print repr(res)
#plt.plot(range(len(res)),len(res)*[0.],'--')
#perc_range=1./(1.+np.exp(-np.arange(-5,3,.15)))
def perf_eval(p):
    fin_res=[]
    real_perc_range=[]
    real_perc,size,res=estim_diff(percent=p/2000.)
    real_perc_range.append(real_perc)
    f=open('/agbs/cluster/naji/Linear Filters/Echo/out/32/'+str(p)+'.txt','w')
    print (p/3000.,file=f)
    print (np.mean(res>0),file=f)

parser = ArgumentParser(description='--t: Number of process', formatter_class=RawTextHelpFormatter)
parser.add_argument('--t',required=True)
args = vars(parser.parse_args())
task_number=int(eval(args['t']))
perf_eval(task_number)
