import sys
from scipy import *
import os
which = range(300)
En_list = []
doping_list = []
do_list = []
S = []
for iR in which:
	if iR%10==0:
		#if iR%100==0:
		r = iR/100.
		s = "%1.1f" %r
		S.append(s)
		#else:
		r = iR/100.
		s = "%1.2f" %r
		S.append(s)
	else:
		r = iR/100.
		s = "%1.2f" %r
		S.append(s)


#print ('2.00' in S)
#sys.exit()
#print S

nlow=80 #lowest itt from which you take average 
print "nlow = ",nlow

for s in S:
	#if os.path.exists('954'+str(i)+'naivephya/data'):
	#	myfile = open('954'+str(i)+'naivephya/data','r')
	nitt = 0
	En_res = 0.0
	nf_res = 0.0
	if os.path.exists(s+'/data.out'):
		print s
		myfile = open(s+'/data.out','r')
		data = myfile.readlines()
		for iline,line in enumerate(data):
			word = line.split()
			if word and word[0]=="itt=":
				itt = int(word[1])
				if itt > nlow:
					last_line = data[iline+5].split()
					nf = 0.5-float(data[iline+4].split()[1])
					En = float(data[iline+3].split()[1])
					EU = float(data[iline+2].split()[1])
					U = float(last_line[1]) #U
					do=EU/U #double occupancy
					mu = float(last_line[3])/(0.5*U) #mu
					nitt += 1
		if nitt > 0:
			#print "data = ", s,	"	number of data picked = ",nitt
			En_list.append([U,En])
			do_list.append([U,do])
			doping_list.append([U,nf])

		else:
			continue

if En_list:
	En_list = array(En_list)
	do_list = array(do_list)
	doping_list = array(doping_list)
	#print ER
	savetxt('do_U.dat',do_list,delimiter='	', fmt ="%.3f %0.6f")
	savetxt('en_mu.dat',En_list,delimiter='	', fmt ="%.3f %0.6f")
	savetxt('doping_mu.dat',doping_list,delimiter='	', fmt ="%.3f %0.6f")
else:
	print "No data found, EXIT!"
