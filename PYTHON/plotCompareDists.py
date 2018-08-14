#!/usr/bin/python
import matplotlib.pyplot as plt
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

dir="fixQ"
myplot(dir+"/Wunc_hist")
myplot(dir+"/Wcoup_hist")
plt.gca().set_prop_cycle(None) #reset color cycle

#plot AS data
dir="fixQ1.9"
myplot(dir+"/Wunc_hist",':')
myplot(dir+"/Wcoup_hist",':')

plt.legend(["$P_u$","$P_c$","$\cos\phi/m_O = -1.9$"],handletextpad=0.1,loc="upper left")
plt.xlabel("Frequency, $\omega$ (cm$^{-1}$)")
plt.ylabel("P(\omega)")
plt.xlim(2900,3900)
plt.ylim(0,0.003)

plt.savefig('compareDists.pdf', bbox_inches='tight')

#plt.tight_layout()
#plt.show()

plt.close()
