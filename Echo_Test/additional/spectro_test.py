#!/usr/bin/env python
from matplotlib import pyplot as plt
from scipy.io import wavfile
import numpy as np
from scipy.signal import decimate,chirp
from scipy.fftpack import fft, ifft,fftfreq
input_dir='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/'
#s=wavfile.read(input_dir+'Linchirp.wav')
#x=s[1][:200000]
rate=10000
t=np.linspace(0,5,rate*5)
x=chirp(t, f0=2000, t1=5, f1=5000, method='linear')
print x.shape
#x = s1 + s2 + nse # the signal
NFFT = 1000       # the length of the windowing segments
Fs = 1. # the sampling frequency

# Pxx is the segments x freqs array of instantaneous power, freqs is
# the frequency vector, bins are the centers of the time bins in which
# the power is computed, and im is the matplotlib.image.AxesImage
# instance
wavfile.write(input_dir+'chirp_test.wav',rate,x)
ax1 = plt.subplot(211) 
plt.plot(t, x)
plt.subplot(212, sharex=ax1)
Pxx, freqs, bins, im = plt.specgram(x, NFFT=NFFT, Fs=Fs, noverlap=900,cmap=plt.cm.gist_heat)
plt.figure()
#print fftfreq(x.shape[0])
plt.plot(fftfreq(x.shape[0]),abs(fft(x)))
plt.show()