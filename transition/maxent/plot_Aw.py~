"""
2/1/2016
@When Sig@real_freq is given
calculate G_ii(w)
"""

import sys
import os
from scipy import *
from scipy import integrate
from scipy import interpolate
from pylab import *

def ftn_to_array(omega,fn):
	"""
	Get (Noms,3) structure array from a complex fn
	with 0 oms and 1,2 as real and imag respectively
	"""
	return array([omega,real(fn),imag(fn)]).transpose()


def Cnumber_Sigs(Sigs):
	Sigs_re = Sigs[1]
	Sigs_im = Sigs[2]
	return Sigs_re+1j*Sigs_im

def Cnumber_Sigs_Array(Sigs_site):
	"""
	Sigs_site = (Nsite,3,Nw_small)
	"""
	Sig_iw_i = zeros((len(oms_f_log),Nsite),complex)
	for i in range(Nsite):
		Sigs_re = Sigs_site[i,1]
		Sigs_im = Sigs_site[i,2]
		Sig_iw_i[:,i] =  Sigs_re+1j*Sigs_im
	return Sig_iw_i


def Cmp_Spectra():
	"""
	@inputs:
		Sigs_mit = (Nsite,3,Nw_small)
	@global inputs
		mu = given chemical potential
		Sig_met = (3,Nw_small) = [oms, Sig_met.re, Sig_met.im]
		Sig_ins = (3,Nw_small) = [oms, Sig_ins.re, Sig_ins.im]
		Nsite = number of sites
		emesh = 2D energy mesh (-2/3,2/3) = (-4t,4t)
		dos = 2D DOS(emesh)
		oms1 = real frequency (from real_sig_met)
	@local variables
		Sig_met(ins)_v = cnumber of Sig_mit(ins)
		t_mat = transfer matrix (off-diagonal Hamiltonian)
		Gew_mit = G(e,iw) = 1/(jw+mu-e-Sig-F)
		z_met(ins) = jw+mu-e-Sig_met(ins)_v
		F_met = 0.5*(z_met-sqrt(z_met-1./9))
		D = diagoinal part of Hamilitonian 
		G_inv = jw+mu-e-Sig-F
		
	@OUTPUT: Aw = (noms,Nsite)
	@Energy Unit = 6t
	"""
	
	def F_z(z):
		"""
		1/F_m(i) = iw+mu-e-Sig_m(i)-t**2*F_m(i)
		->F_m = 1/t*(a-sqrt(a**2-1))
			where a = (iw+mu-e-Sig)/(2t)
		Then 
		t**2*F_m = 3t*(z-sqrt(z**2-1/9)) 
						 = 6t*0.5*(z-sqrt(z**2-1/9))
		where z = (iw+mu-e-Sig)/(6t)
		"""
		if real(z)>=0:
			return 0.5*(z-sqrt(z**2-1./9))
		else:
			return 0.5*(z+sqrt(z**2-1./9))
	Sigs_site =zeros((noms,Nsite),complex)

	for i in range(Nsite):
		sigdata=loadtxt('sigs/sig.'+str(i)+'/'+sigfile).T
		sig_oo=loadtxt('sigs/Sig.out.'+str(i))[-1,1]
		oms=sigdata[0]
		print "mu=%2.5f, sig_oo=%2.5f" %(mu,sig_oo)
		Sigs_site[:,i]=sigdata[1]+sigdata[2]*1j+sig_oo
		
	Sig_met_v = Cnumber_Sigs(Sig_met)+Sig_m_oo
	Sig_ins_v = Cnumber_Sigs(Sig_ins)+Sig_i_oo
	
	t_mat = zeros((Nsite,Nsite),float) #hopping matrix
	for i in range(Nsite-1):
		t_mat[i,i+1]= -1./6 #energy unit is 6t(half-bandwidth)
		t_mat[i+1,i]= -1./6
	
	Aw=zeros((noms,Nsite),float)
	for iw,w in enumerate(oms1):
		Gew_mit = zeros((len(emesh),Nsite,Nsite),complex) #G(e,iw) / emesh=2D / Nsite=# of sites
		for ie,e in enumerate(emesh):
			z_met = w+mu-e-Sig_met_v[iw]
			z_ins = w+mu-e-Sig_ins_v[iw]
			D = diag(w+mu-e-Sigs_site[iw])
			G_inv = D - t_mat 
			G_inv[0,0] += -F_z(z_met)
			G_inv[Nsite-1,Nsite-1] += -F_z(z_ins)
			Gew_mit[ie] = linalg.inv(G_inv) #<- the most time consuming part
		for i in range(Nsite):
			G_E = Gew_mit[:,i,i]
			Gwloc_mit = integrate.simps(G_E*dos,emesh)
			Aw[iw,i]=-1./pi*imag(Gwloc_mit)
		print iw,Aw[iw,1]

	return Aw

if __name__ == '__main__':

##############Options & Global variables#####################
	"""
	@All iw-dependent variables will be first evaluated with
	small mesh (logged tail) of of iw_n then splined to the regular(large) mesh.
	"""
	
	mu=0.85*0.5*2.13 #mu=a*U/2
	
	site_list=range(20)
	Nsite=len(site_list)
	
	data = loadtxt('DOS_2D').T
	emesh = data[0]
	dos = data[1]
	
	sigfile='real_sig_me.out'
	Sig_met = loadtxt('sigs/sig_m/'+sigfile).T
	Sig_m_oo = loadtxt('sigs/Sig_m.out')[-1,1]
	Sig_ins = loadtxt('sigs/sig_i/'+sigfile).T
	Sig_i_oo = loadtxt('sigs/Sig_i.out')[-1,1]
	oms1=Sig_met[0]
	noms=len(oms1)

	sdir="spectral_data/"
	if not os.path.exists(sdir):
		os.mkdir(sdir)
	Aw = Cmp_Spectra()
	for i in range(Nsite):
		plot(oms1,Aw[:,i],label=str(i))
		savetxt(sdir+"spectral."+str(i),array([oms1,Aw[:,i]]).T)

	legend()
	show()


