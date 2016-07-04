"""
3/25/2016
"""

import sys
import os
from scipy import *
from scipy import integrate
from scipy import interpolate
from pylab import *

#from eval_Ekin import Cmp_Ekin_Sigma#(emesh,dos,freq,Sigma,mu,beta)


def ftn_to_array(omega,fn):
	"""
	Get (N,3) structure array from a complex fn
	with 0 oms and 1,2 as real and imag respectively
	"""
	return array([omega,real(fn),imag(fn)]).transpose()


def Cnumber_Sigs(Sigs):
	Sigs_re = Sigs[1]
	Sigs_im = Sigs[2]
	return Sigs_re+1j*Sigs_im


def fermi(e,beta):
  return 1./(exp(beta*e)+1.)


def Cmp_Gloc_Delta(Sigs):
	oms1 = Sigs[0] #supposed to be log mesh
	if len(oms1)!=len(oms_ind):
		print "frequency number does not match, EXIT!"
		sys.exit()
	
	#array structure to complex number
	Sigma = Cnumber_Sigs(Sigs)

	Delta=zeros((len(Sigs[0]),3),dtype=float)
	Delta[:,0] = Sigs[0]
	Glocal=zeros((len(Sigs[0]),3),dtype=float)
	Glocal[:,0] = Sigs[0]
			
	#Retrieve old Delta from Delta.dat (regular mesh)
	Dfile = loadtxt(fileDelta)
	Delta_old = zeros( (len(oms_ind), shape(Dfile)[1]), dtype=float )
	for i in range(len(oms_ind)):
		Delta_old[i,:] = Dfile[oms_ind[i],:]

	oms_limit =  2*pi/beta*Dlimit
	
	for iw,w in enumerate(oms1):
		z = w*1j+mu-Sigma[iw]
		G = 1./(z-emesh)
		Gloc = integrate.simps(G*dos,emesh)

		#Replacing tail with its Coeff/iw for w > oms_limit.
		if w < oms_limit:
			Delt = z-1./Gloc
		else:
			Delt = -1j/w*Coeff_Delta
		#mix with old data
		Delt = Delt*(1.-mixr) + (Delta_old[iw,1]+1j*Delta_old[iw,2])*mixr
		
		Glocal[iw,1:] = array([Gloc.real,Gloc.imag])
		Delta [iw,1:] = array([Delt.real,Delt.imag])
	
	return Delta,Glocal

if __name__ == '__main__':

##############Options & Global variables#####################
	U = 1.975
	mu = 0.95*U/2.
	
	data = loadtxt('DOS').T
	emesh = data[0]
	dos = data[1]

	
	rSig=loadtxt('sig/data1/real_sig_me.out')
	oms=rSig[:,0]
	sig_oo=loadtxt('../Sig.out')[-1,1]
	print "mu= %2.5f, sig_oo=, %2.5f" %(mu,sig_oo)
	
	#array structure to complex number
	Sigma = rSig[:,1]+1j*rSig[:,2]+sig_oo

	Glocal=zeros((len(oms),2),dtype=float)
	Glocal[:,0] = oms
			
	for iw,w in enumerate(oms):
		z = w+mu-Sigma[iw]
		G = 1./(z-emesh)
		Gloc = integrate.simps(G*dos,emesh)
		Glocal[iw,1] = -1*(Gloc.imag)/pi

	savetxt('Aw.dat',Glocal)

