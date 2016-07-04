from scipy import *
import os
import sys
from pylab import *

Nstart=int(sys.argv[1])
Nlimit=int(sys.argv[2])

Nsite=20
do_is_it = zeros((Nsite,Nlimit),float)
if os.path.exists('data'):
	myfile = open('data','r')
	data = myfile.readlines()
	U=float(data[0].split()[2][:-1])
	print "#U=",U
	for iline, line in enumerate(data):
		word=line.split()
		if word:
			if word[0]=="index=":
				itt = int(word[1])
				isite = int(word[2])
				word2 = data[iline+1].split()
				do = float(word2[1])/U
				if itt<Nlimit:
					do_is_it[isite,itt] = do

				#do_L.append(doping)
				#which.append(itt)

do_ave=zeros(Nsite,float)
for ii,i in enumerate(range(-Nsite/2,Nsite/2)):
	#do_ave[ii]=sum(do_is_it[ii])/Nitt
	do_ave[ii]=sum(do_is_it[ii][Nstart:])/(Nlimit-Nstart)
	print i,do_ave[ii]


#if len(sys.argv)>1:
#	imax=99
#	imin=80
#	step=3
#	for i in range(imin,imax)[::step]:
#		plot(range(-Nsite/2,Nsite/2),do_is_it[:,i],'.-',label=str(i))
#	legend()
#	show()
