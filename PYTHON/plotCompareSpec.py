#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os.path

def myplot(filename,myline='-'):
    if not os.path.isfile(filename):
        raise ValueError(filename+ ' file doesn''t exist')
    data=np.genfromtxt(filename)
    max=np.amax(data[:,1])
    data[:,1]/=max
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
myplot("nve/spec")
myplot("nve_defaultAvgF/spec")
legList=["3415 cm$^{-1}$","3449 cm$^{-1}$"]

#plot AS data
#legList=["Mine","AS2008"]
#plt.gca().set_prop_cycle(None) #reset color cycle
#myplot('../AuerSkinnerData/IR_theory.txt',':')

plt.legend(legList,handletextpad=0.1,loc="upper left")
plt.xlabel("Frequency, $\omega$ (cm$^{-1}$)")
plt.ylabel("Intensity $I(\omega)$ (arb.)")
plt.xlim(2000,4500)
plt.ylim(0.0,1.05)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    right=False,      # ticks along the bottom edge are off
    left=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off

plt.savefig('compareSpec_avgF.pdf', bbox_inches='tight')

#plt.tight_layout()
#plt.show()

plt.close()
