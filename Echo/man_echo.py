# -*- coding: utf-8 -*-
from __future__ import print_function
import wave

import sys,os
#a line to add Codes path to be able to use mod_welch module!
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '..'))
from mod_welch import welch


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
#from scipy.signal import welch
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
output_dir=main_path+'Sounds/'+song_name+'/Experiments/Room/1/'


#input_file_names=(listdir(input_dir+'Segments/'))
#output_file_names=(listdir(output_dir+'Segments/'))
def estim_diff(input_sig, input_seg_len, output_sig, output_seg_len, nperseg, num_of_seg, nperseg_step):
    """The difference between estimators in both directions. 
    inputs:
        seg_num:        number of segments to divide the music piece into. In practice
                        we always divide the music to seg_num+1 segments ignoring to 
                        last segment.
                                
        nperseg:        initial Welch method window size
    
    """
    res=np.empty(num_of_seg)
    
    for i in range(num_of_seg):
        # converting int16 wav file to a float signal and cutting the signal to 
        # pieces of size seg_len. Also selecting the mono signal

        temp_input_sig=pcm2float(input_sig[input_seg_len*i:input_seg_len*(i+1)],'float32')[:,0]
        temp_output_sig=pcm2float(output_sig[output_seg_len*i:output_seg_len*(i+1)],'float32')[:,0]
        t=time()
        S_inp=welch(temp_input_sig,nperseg=nperseg*nperseg_step)[1]    
        S_out=welch(temp_output_sig,nperseg=nperseg*nperseg_step)[1]    

        res[i]=delta_estimator_3(S_out/S_inp,S_inp)-delta_estimator_3(S_inp/S_out,S_out)
        #out=float2pcm(output_sig,'int16')
    return res
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

#list of number of pieces that a music gets chopped into
num_of_segs=[8,16,32,64,128]

#size of increment for nperseg in each step
nperseg_step=500

def perf_eval(param):
    # wrtitten in 3000 basis, finding the nperseg value from param
    nperseg=param % 3000
    # wrtitten in 3000 basis, finding the number of music segments
    num_of_seg_idx=(param-nperseg)/3000
    num_of_seg=num_of_segs[num_of_seg_idx]

    input_rate,input_sig=wavfile.read(input_dir+song_name+'.wav')
    output_rate,output_sig=wavfile.read(output_dir+song_name+'.wav')
    
    #the +1 in denominator is because we exclude the last piece of music to only
    # consider music pieces of the same size.
    input_seg_len=input_sig.shape[0]/(num_of_seg+1)
    output_seg_len=output_sig.shape[0]/(num_of_seg+1)
    
    if input_rate!=output_rate:
        print ("Rate Mistmatch!")
        sys.exit(0)
    
    #print (nperseg,nperseg_step,input_seg_len,output_seg_len)
    if np.min((input_seg_len,output_seg_len))*0.7 < nperseg * nperseg_step:
        print ("Nothing to do!")
        sys.exit(0)
        
    res=estim_diff(input_sig, input_seg_len, output_sig, output_seg_len, nperseg, num_of_seg, nperseg_step)

    
    f=open('/agbs/cluster/naji/Linear Filters/Echo/out/Winter/Room/'+str(num_of_seg)+'/'+str(nperseg)+'.txt','w')
    print (nperseg,file=f)
    print (np.mean(res>0),file=f)

parser = ArgumentParser(description='--t: Number of process', formatter_class=RawTextHelpFormatter)
parser.add_argument('--t',required=True)
args = vars(parser.parse_args())
task_number=int(eval(args['t']))
t=time()
#for task_number in np.arange(1000,1500):
#    if not os.path.isfile('/agbs/cluster/naji/Linear Filters/Echo/out/Winter/Room/8/'+str(task_number)+'.txt'):
#        print ("For "+str(task_number)+" we are in")
perf_eval(task_number)
print (time()-t)
