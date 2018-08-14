#!/usr/bin/python
import pylab as plt
import numpy as np
import os.path

width=3.25 #inches
font=11
fig=plt.figure(figsize=(width,width*0.7))
fig.patch.set_facecolor('white')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=font)
plt.rc('legend', fontsize=font-2)

dir="mu0.1p2.5"
file=dir+"/tcf"
if not os.path.isfile(file):
    raise ValueError(file+ ' file doesn''t exist')
data=np.genfromtxt(file)
plt.plot(data[:,0],data[:,1])
plt.plot(data[:,0],data[:,2])

plt.legend(["Real","Imag"],handletextpad=0.1)
plt.xlabel("Time, $t$ (ps)")
plt.ylabel("Time Correlation Function")
#plt.xlim(2900,3900)
#plt.ylim(0,40)

plt.savefig('tcf.pdf', bbox_inches='tight')

#plt.tight_layout()
#plt.show()

plt.close(fig)
