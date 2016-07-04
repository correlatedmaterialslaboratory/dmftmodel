import sys
from scipy import *
import os
from pylab import *
#nf_L = []
#which = []

Nitt=int(sys.argv[1])

Nsite=20
ilimit = Nitt-10
nf_is_it = zeros((Nsite,Nitt),float)
if os.path.exists('data'):
	myfile = open('data','r')
	data = myfile.readlines()
	for iline, line in enumerate(data):
		word=line.split()
		if word:
			if word[0]=="index=":
				itt = int(word[1])
				isite = int(word[2])
				word2 = data[iline+2].split()
				nf = 0.5-float(word2[1])
				if itt<Nitt:
					nf_is_it[isite,itt] = nf

				#nf_L.append(doping)
				#which.append(itt)

nf_ave=zeros(Nsite,float)
for ii,i in enumerate(range(-Nsite/2,Nsite/2)):
	#nf_ave[ii]=sum(nf_is_it[ii])/Nitt
	nf_ave[ii]=sum(nf_is_it[ii][ilimit:])/(Nitt-ilimit)
	print i,nf_ave[ii]


#if len(sys.argv)>1:
#	imax=99
#	imin=80
#	step=3
#	for i in range(imin,imax)[::step]:
#		plot(range(-Nsite/2,Nsite/2),nf_is_it[:,i],'.-',label=str(i))
#	legend()
#	show()
