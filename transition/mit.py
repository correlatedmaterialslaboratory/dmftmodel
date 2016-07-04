"""
7/31/2015
@Starting to develop the main code to calculate 
metal-intermediate-insulator
system to see the critical length.
"""

import sys
import os
from scipy import *
from scipy import integrate
from scipy import interpolate
from pylab import *

from ftns_Fourier import Create_om_mesh,Spline_Real
from include_DMFT import IMP

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


def Cmp_Delta_site(Sigs_mit,itt=0):
	"""
	@inputs:
		Sigs_mit = (Nsite,3,Nw_small)
	@global inputs
		mu = given chemical potential
		Sig_met = (3,Nw_small) = [oms, Sig_met.re, Sig_met.im]
		Sig_ins = (3,Nw_small) = [oms, Sig_ins.re, Sig_ins.im]
		Nsite = number of sites
		oms_f_log = fermion frequency with log tail
		oms_ind = log tail index for regular mesh
		emesh = 2D energy mesh (-2/3,2/3) = (-4t,4t)
		dos = 2D DOS(emesh)
		mixr = mixing ratio with old file
	@local variables
		oms1 = oms_f_log (from Sigs_mit[0])
		Sig_met(ins)_v = cnumber of Sig_mit(ins)
		t_mat = transfer matrix (off-diagonal Hamiltonian)
		Gew_mit = G(e,iw) = 1/(jw+mu-e-Sig-F)
		z_met(ins) = jw+mu-e-Sig_met(ins)_v
		F_met = 0.5*(z_met-sqrt(z_met-1./9))
		D = diagoinal part of Hamilitonian 
		G_inv = jw+mu-e-Sig-F
		
	@OUTPUT: Delta_site = (Nsite,Nw_small,3)
	@energy unit = 6t
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

	oms1 = Sigs_mit[0,0] #supposed to be log mesh for its tail.
	if len(oms1)!=len(oms_ind):
		print "frequency number does not match, EXIT!"
		sys.exit()

	Sig_mit_v = Cnumber_Sigs_Array(Sigs_mit)
	#print shape(Sig_mit_v)

	Sig_met_v = Cnumber_Sigs(Sig_met)
	if False: #test
		Fz = F_z(1j*oms1+mu-Sig_met_v)
		savetxt('F_z_met.dat', array([oms1,Fz.real,Fz.imag]).T)
		sys.exit()
	Sig_ins_v = Cnumber_Sigs(Sig_ins)
	
	t_mat = zeros((Nsite,Nsite),float) #hopping matrix
	for i in range(Nsite-1):
		t_mat[i,i+1]= -1./6 #energy unit is 6t(half-bandwidth)
		t_mat[i+1,i]= -1./6

	Delta_site = zeros((Nsite,Nw_small,3),float)
	
	#f = open('Delta_ii_2D.dat','w')
	for iw,w in enumerate(oms1):
		Gew_mit = zeros((len(emesh),Nsite,Nsite),complex) #G(e,iw) / emesh=2D / Nsite=# of sites
		for ie,e in enumerate(emesh):
			z_met = 1j*w+mu-e-Sig_met_v[iw]
			z_ins = 1j*w+mu-e-Sig_ins_v[iw]
			D = diag(1j*w+mu-e-Sig_mit_v[iw])
			G_inv = D - t_mat 
			#if ie==0:
			#	print w,abs(G_inv[0,0]),abs(G_inv[Nsite-1,Nsite-1]),abs(F_z(z_met)),abs(F_z(z_ins))
			G_inv[0,0] += -F_z(z_met)
			G_inv[Nsite-1,Nsite-1] += -F_z(z_ins)
			Gew_mit[ie] = linalg.inv(G_inv) #<- the most time consuming part
		for i in range(Nsite):
			G_E = Gew_mit[:,i,i]
			Gwloc_mit = integrate.simps(G_E*dos,emesh)
			Delta_v = 1j*w+mu-Sig_mit_v[iw,i]-1/Gwloc_mit
			Delta_site[i,iw,0] = w
			Delta_site[i,iw,1] = Delta_v.real
			Delta_site[i,iw,2] = Delta_v.imag
			#if i == 0:
			#	print >>f, w, Delta_.real,Delta_.imag, iw
	#f.close()


	#mixing with old Delta file if exists
	fileDelta = 'Delta_site.dat'
	if (os.path.exists(fileDelta)):
		Dfile = loadtxt(fileDelta).T
		Delta_site_old = zeros((Nsite,Nw_small,3),float)
		for isite in range(Nsite):
			Delta_i = zeros((Nw_small,3),float)
			Delta_i[:,0]=Dfile[0]
			Delta_i[:,1]=Dfile[2*isite+1]
			Delta_i[:,2]=Dfile[2*isite+2]
			Delta_site_old[isite] = Delta_i

		Delta_site = Delta_site*(1-mixr)+mixr*Delta_site_old

	Dfile = zeros((Nw_small,2*Nsite+1),float)
	for isite in range(Nsite):
		Dfile[:,0] = oms1
		Dfile[:,2*isite+1] = Delta_site[isite,:,1]
		Dfile[:,2*isite+2] = Delta_site[isite,:,2]
	savetxt("Delta_site.dat",Dfile)
	if (itt+1)%10==0:
		savetxt("Delta_site.dat."+str(itt),Dfile)

	return Delta_site

if __name__ == '__main__':

##############Options & Global variables#####################
	"""
	@All iw-dependent variables will be first evaluated with
	small mesh (logged tail) of of iw_n then splined to the regular(large) mesh.
	"""
	
	op_U = False #if false: EXIT!
	op_beta = False	#if false: default beta = 100.
	op_mu = False #if false: default mu = U/2
	op_Nitt = False #if false: default Nitt = 15
	op_Mstep = False #if false: default Mstep = 10e6
	op_mixr = False #if false: default mixr = 0.5
	op_Nsite = False #if false: default mixr = 0.5

	options = sys.argv

	for iop, op in enumerate(options):
		if op=='U': 
			U = float(options[iop+1])
			op_U = True
		if op=='beta': 
			beta = float(options[iop+1])
			op_beta = True
		if op=='mu': 
			mu = float(options[iop+1])*U/2.
			op_mu = True
		if op=='Nitt':
			Nitt = int(options[iop+1])
			op_Nitt = True
		if op=='mixr':
			mixr = float(options[iop+1])
			op_mixr = True
		if op=='Mstep':
			Mstep = float(options[iop+1])
			op_Mstep = True
		if op=='Nsite': 
			Nsite = int(options[iop+1])
			op_Nsite = True


	if not op_U:
		print "U is not given, EXIT!!"
		sys.exit()
	if not op_beta:
		print "No beta value given, take default value B=100."
		beta = 100.
	
	if not op_mu:
		print "No chemical potential given: half-filled."
		mu = U/2.

	if not op_Nitt:
		print "No initial Nitt: default value"
		Nitt = 15

	if not op_mixr:
		print "No mixr: default value"
		mixr = 0.5
	
	if not op_Mstep:
		print "No Mstep: default value 2e6"
		Mstep = 2e6
	
	if not op_Nsite:
		print "No Nsite: default value Nsite=20"
		Nsite = 20
	
	
	Nomega = 5000
	nom = 100 #number of sampling in ctqmc
	nom_tail = 50 #number of tail
	oms_f, oms_f_log, oms_b, oms_b_log, oms_equal_f,oms_ind \
			= Create_om_mesh(beta, Nomega, nom, nom_tail) 
	#oms_f(b) = regular mesh
	#oms_f(b)_log = small mesh with log tail
	#oms_ind = index of small mesh in terms of regular mesh
	Nw_small = len(oms_f_log)
	Nw = len(oms_f)
	#load density of states
	data = loadtxt('DOS_2D').T
	emesh = data[0]
	dos = data[1]

	#print options
	print "U = %2.3f, mu = %2.3f, B = %2.3f, Nitt = %d, Mstep = %2.1f ,mixr = %2.2f" %(U,mu,beta,Nitt,Mstep,mixr)
	print '# mesh_omega = ', len(oms_f)
	print '# small_mesh = ', len(oms_f_log)

############# initial-DMFT ###############################


	Sig_met = loadtxt('Sig_m.out').T
	Sig_ins = loadtxt('Sig_i.out').T

	fileDelta = 'Delta_site.dat'
	Sigs_site = zeros((Nsite,3,Nw_small),float)
	if (os.path.exists(fileDelta)):
		Dfile = loadtxt(fileDelta).T
		print "Old state exists!"
		print
		
		Delta_site = zeros((Nsite,Nw_small,3),float)
		for isite in range(Nsite):
			Delta_i = zeros((Nw_small,3),float)
			Delta_i[:,0]=Dfile[0]
			Delta_i[:,1]=Dfile[2*isite+1]
			Delta_i[:,2]=Dfile[2*isite+2]
			Delta_site[isite] = Delta_i

	else:
		print "No old file."
		print 

		for isite in range(Nsite):
			if isite < Nsite/2:
				Sigs_site[isite] = Sig_met #half with metallic sols
			else:
				Sigs_site[isite] = Sig_ins #half with insulating sols

		Delta_site = Cmp_Delta_site(Sigs_site)
		#savetxt('Delta_1.dat',Delta_site[0])
		#savetxt('Delta_2.dat',Delta_site[1])
		##check point 8/10/2015##
	#sys.exit()

###############self-consistent loop################
	ctqmc = IMP(nom,Nsite,beta)
	p_ener = 0.0
	for itt in range(Nitt):
		for isite in range(Nsite):
			Delta = Delta_site[isite]
		
			(Sigs,Sigs_regular, n_imp, TrSigmaG) = \
					ctqmc.Run(U, mu, oms_f, oms_ind, isite,Delta,Mstep)
			ctqmc.copyfiles(isite,itt)
			Sigs_site[isite] = Sigs

			EU = TrSigmaG # = 1/2*Tr(Sigma*G)
			print 
			print "index= ",itt,isite
			print "EU= %2.7f" %EU
			print "nf= %2.7f" %(n_imp/2)
			print
			#print "U= %2.6f mu= %2.6f T= %2.6f Etotal= %2.6f" %(U,mu,1/beta,Etotal)
		
		#New Delta
		Delta_site = Cmp_Delta_site(Sigs_site,itt)

