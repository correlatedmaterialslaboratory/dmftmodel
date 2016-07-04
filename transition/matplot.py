from scipy import *
import matplotlib.pyplot as plt

savepng="MIT_mu1.0_U2.05"

fig, ax1 = plt.subplots()
fig.suptitle(savepng, fontsize=14, fontweight='bold')


[x,y2]=loadtxt('docc.dat').T
ax1.plot(x, y2, 'r.-')
ax1.set_ylabel('docc', color='r')

#plt.show()
F = plt.gcf()
#F.savefig(savepng+".png",dpi = 250,format='png',bbox_inches='tight')
F.savefig(savepng+".pdf",dpi = 300,format='pdf',bbox_inches='tight')
