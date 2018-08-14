#!/usr/bin/python
import pylab as plt
import numpy as np
import os.path

def myplot(filename,myline='-'):
    if not os.path.isfile(filename):
        raise ValueError(filename+ ' file doesn''t exist')
    data=np.genfromtxt(filename)
    line, =plt.plot(data[:,0],data[:,1],linestyle=myline)
    return line

width=3.25 #inches
font=11
fig=plt.figure(figsize=(width,width*0.7))
fig.patch.set_facecolor('white')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=font)
plt.rc('legend', fontsize=font-2)

#plot uncoupled for only one mu, b/c it's indep
myplot("mu0.1p2.5/Wunc_hist").set_color('k')

#plot AS data
myplot('../AuerSkinnerData/Pu.txt',':').set_color('k')

plt.legend(["Mine","AS2008"],handletextpad=0.1,loc="upper left")
plt.xlabel("Frequency, $\omega$ (cm$^{-1}$)")
plt.ylabel("P(\omega)")
plt.xlim(2900,3900)
plt.ylim(0,0.003)

plt.savefig('comparePu.pdf', bbox_inches='tight')

#plt.tight_layout()
#plt.show()

plt.close()
