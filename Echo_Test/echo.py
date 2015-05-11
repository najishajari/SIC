# -*- coding: utf-8 -*-
#import pyaudio
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
from scipy import signal
from multidim_welch import ndim_welch

def delta_estimator_1(h,S_x):
    return (np.mean(np.multiply(h,S_x))-(np.mean(h)*np.mean(S_x)))/(np.mean(np.multiply(h,S_x)))
def delta_estimator_2(h,S_x):
    S_x=np.ma.masked_invalid(S_x)    
    return (np.mean(np.multiply(h,S_x))-(np.mean(h)*np.mean(S_x)))/(np.sum(h)*np.max(S_x))
def delta_estimator_3(h,S_x):
    S_x=np.ma.masked_invalid(S_x)    
    h=np.ma.masked_invalid(h)    
    int_h_S_x=np.mean(np.multiply(h,S_x))
    int_h=np.mean(h)
    int_S_x=np.mean(S_x)
    return np.log(int_h_S_x)-np.log(np.mean(int_h))-np.log(int_S_x)
def delta_estimator_4(h,S_x):
    S_x=np.ma.masked_invalid(S_x)    
    h=np.ma.masked_invalid(h)    
    int_h_S_x=np.mean(np.multiply(h,S_x))
    int_h=np.mean(h)
    int_S_x=np.mean(S_x)
    return float(int_h_S_x)/(float(np.mean(int_h))*float(int_S_x))
    
def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

file_path='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/'
#file_path='/Users/Naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/'
def SIC_echo(sound_piece):
    sound_path=file_path+'Sounds/'+sound_piece+'/Synthetic Echo'
    sound_file_names=(listdir(sound_path+'/Segments/'))
    #input_dir='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/'
    IR_path=file_path+'Places/Elevden Hall/'
    IR_file_names=listdir(IR_path)
    
    sound_counter=0
    res=np.empty((len(sound_file_names),len(IR_file_names)))
    for sound_file_name in sound_file_names:
        IR_counter=0
        for IR_file_name in IR_file_names:
            #if in_counter==2:
            #    break
            #print sound_path+sound_file_name
            print sound_path+'/Segments/'+sound_file_name,IR_file_name
            input_rate,input_sig=wavfile.read(sound_path+'/Segments/'+sound_file_name)
            input_sig=pcm2float(input_sig,'float32')
    
            #test: artificial signals
            #input_rate=44100
            #input_sig=np.random.randint(low=1,high=16,size=(1000,2))   
            
            IR_rate,IR_sig=wavfile.read(IR_path+IR_file_name)
            IR_sig=pcm2float(IR_sig,'float32')
            
            if input_rate!=IR_rate:
                print "Size mismatch"
                sys.exit(-1)
            else:
                rate=input_rate
            print sound_file_name
            con_len=-1
            out_0=fftconvolve(input_sig[:,0],IR_sig[:,0])
            out_0=out_0/np.max(np.abs(out_0))
            out_1=fftconvolve(input_sig[:,1],IR_sig[:,1])
            out_1=out_1/np.max(np.abs(out_1))
            output_sig=np.vstack((out_0,out_1)).T
            #S_inp=np.abs(fft(input_sig[:,0]))**2
            #S_out=np.abs(fft(out[:input_sig[:,0].shape[0],0]))**2
            min_size=np.min((input_sig[:,0].shape[0],output_sig[:,0].shape[0]))
    
            S_inp=ndim_welch(input_sig[:min_size,0][None,...],nperseg=min_size)[1]    
            S_out=ndim_welch(output_sig[:min_size,0][None,...],nperseg=min_size)[1]    
                    
            res[sound_counter,IR_counter]=delta_estimator_4(S_out/S_inp,S_inp)-delta_estimator_4(S_inp/S_out,S_out)
            out=float2pcm(output_sig,'int16')
            out[:input_sig[:,0].shape[0]]=out[:input_sig[:,0].shape[0]]+input_sig
            wavfile.write(sound_path+'/echoed outputs/'+sound_file_name+'_'+IR_file_name+'_echoed.wav',rate,out)
            IR_counter+=1
        sound_counter+=1
    return IR_file_names, sound_file_names, res
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}
plt.rc('font', **font)
plt.rc('text', usetex=True)
tableau20 = np.asarray([(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] )/255.
IR_names,lacrimosa_segments,Lacrimosa_res=SIC_echo('Lacrimosa')
IR_names,winter_segments,Winter_res=SIC_echo('Winter')
fig,ax=plt.subplots(2,1)
for i in range(len(IR_names)):
    ax[0].plot(np.arange(1,len(lacrimosa_segments)+1),Lacrimosa_res[:,i],label=r'${\rm '+IR_names[i]+r'}$',color=tableau20[i*3])
    ax[1].plot(np.arange(1,len(winter_segments)+1),Winter_res[:,i],label=r'${\rm '+IR_names[i]+r'}$',color=tableau20[i*3])
    ax[0].set_ylabel(r'$\Delta^{\infty}_{X_t\to Y_t}-\Delta^{\infty}_{Y_t\to X_t}$')
    ax[0].set_title(r'${\rm Lacrimosa\ movement}$')
    ax[0].set_xlabel(r'${\rm Audio\ segment}$')
    ax[1].set_ylabel(r'$\Delta^{\infty}_{X_t\to Y_t}-\Delta^{\infty}_{Y_t\to X_t}$')
    ax[1].set_xlabel(r'${\rm Audio\ segment}$')
    ax[1].set_title(r'${\rm Winter\ movement}$')


ax[0].legend(bbox_to_anchor=(0.0, 1.05, 1., .102), loc=3, ncol=2, mode="expand")
plt.gcf().set_size_inches(13,15)
plt.savefig(file_path+'Winter,Lacrimosa.eps',transparent=True,bbox_inches='tight')
plt.show()
