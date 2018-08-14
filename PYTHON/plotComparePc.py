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

dirlist=["mu0.1p2.5","mu0.2p2.5","mu1p2.5"]
leglist=["0.16","0.26","1.0"]
leglist=["$\mu_g^\prime = "+s+"$" for s in leglist]
leglist.append('AS2008')

#plot coupled hist for all mu
for dir in dirlist:
    myplot(dir+"/Wcoup_hist")

#plot AS data
myplot('../AuerSkinnerData/Pc.txt',':').set_color('black')

plt.legend(leglist,handletextpad=0.1,loc="upper left")
plt.xlabel("Frequency, $\omega$ (cm$^{-1}$)")
plt.ylabel("P(\omega)")
plt.xlim(2900,3900)
plt.ylim(0,0.003)

plt.savefig('comparePc.pdf', bbox_inches='tight')

#plt.tight_layout()
#plt.show()

plt.close(fig)
