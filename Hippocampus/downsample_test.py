from scipy import signal
import numpy as np
x=np.random.randn(1000)
y=signal.decimate(x,100)
print y
