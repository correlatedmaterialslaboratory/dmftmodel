"""
7/16/2015
@mu_c calculation with fixed U at very low temp like b=400/800
@then nom=200/400
"""

import sys, os
from scipy import *
from scipy import integrate
from scipy import interpolate
from pylab import *

from ftns_Fourier import Create_om_mesh,Spline_Real
from include_DMFT import IMP

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

def Cmp_Ekin_Sigma(Sig_reg):
	"""
	Sig_reg: (3,Nomega) structure

	Ek= 2*Int(D(E)*E*(1/beta)*sum 1/(iw+mu-E-Sigma[iw]))
		
		= 2*Int(D(E)*E*(1/beta)*sum 1/(iw+mu-E-Sigma[iw]))
		 -2*Int(D(E)*E*(1/beta)*sum 1/(iw+mu-E-Sigma_oo))
		 +2*Int(D(E)*E* fermi(E-mu+Sigma_oo))
	"""
	freq = Sig_reg[0] #Take its frequency
	Sigma = Cnumber_Sigs(Sig_reg) #data str to complex number
	Sigma_oo = real(Sigma[-1]) #Get  Sigma(oo)
	#print "Sigma_oo = ",Sigma_oo
	
	DE = zeros(len(emesh),float)
	corr = zeros(len(emesh),float)
	for iE,E in enumerate(emesh):
		An = 1./(1j*freq+mu-E-Sigma)
		Bn = 1./(1j*freq+mu-E-Sigma_oo)
		#print shape(Bn)
		Mat_sum = 2*real(sum(An))/beta
		corr[iE] = -2*real(sum(Bn))/beta +fermi(E+Sigma_oo-mu,beta)
		
		DE[iE] = Mat_sum+corr[iE]
	
	Ek = 2*integrate.simps(dos*emesh*DE,emesh)
	
	return Ek


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
	"""
	@All iw-dependent variables will be first evaluated with
	log mesh(small) of iw_n then splined to the regular mesh(large).

	@
	"""
	
	Cubic = True #if false: Bethe lattice calculation
	op_U = False #if false: EXIT!
	op_beta = False	#if false: default beta = 100.
	op_mu = False #if false: default mu = U/2
	op_Nitt = False #if false: default Nitt = 15
	op_Mstep = False #if false: default Mstep = 10e6
	op_mixr = False #if false: default mixr = 0.5

	options = sys.argv

	for iop, op in enumerate(options):
		if op=='noCubic': Cubic = False
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
		print "No Mstep: default value 10e6"
		Mstep = 1e6
	
	Nomega = int(beta*20)
	nom = int(beta)/2 #number of sampling in ctqmc
	#Nomega = 2000 #total Matsubara number
	#nom = 50 #total Matsubara number
	nom_tail = 50 #number of tail
	oms_f, oms_f_log, oms_b, oms_b_log, oms_equal_f,oms_ind \
			= Create_om_mesh(beta, Nomega, nom, nom_tail) #regular mesh, log mesh of iw_n
	
	if Cubic:
		print "3D cubic lattice + DMFT"
		data = loadtxt('DOS').T
		emesh = data[0]
		dos = data[1]
		##symmetrize to negative energy
		#Ne = len(emesh0)
		#emesh = zeros(2*Ne-1,float)
		#dos = zeros(2*Ne-1,float)
		#for i in range(Ne):
		#	emesh[i] = -emesh0[-(i+1)]
		#	emesh[-(i+1)] = emesh0[-(i+1)]
		#	dos[i] = dos0[-(i+1)]
		#	dos[-(i+1)] = dos0[-(i+1)]
		#savetxt("DOS_",array([emesh,dos]).T)
		#plot(emesh,dos)
		#show()
		#sys.exit()

	else:
		print "Bethe lattice + DMFT"
		emesh = linspace(-1,1.,2000)
		dos = 2./pi*sqrt(1.-emesh**2)
		
	#print options
	print "U = %2.3f, mu = %2.3f, B = %2.3f, Nitt = %d, Mstep = %2.1f ,mixr = %2.2f" %(U,mu,beta,Nitt,Mstep,mixr)
	print '# mesh_omega = ', len(oms_f)


	#For Delta tail, 
	Coeff_Delta = integrate.simps(emesh**2*dos,emesh)
	print "Coeff for Delta tail = ",Coeff_Delta
	Dlimit = 300
	oms_cut = 2*pi/beta*Dlimit
	print "For Delta calculation, oms is cut by %2.2f and replaced by Coeff/iw_n." %(oms_cut)
	#sys.exit()

############# initial-DMFT ###############################
	
	ctqmc = IMP(nom,beta)
	fileDelta = 'Delta.dat'
	if (os.path.exists(fileDelta)):
		print "Old state exists!"
		print 
		#New Delta
		Dfile = loadtxt(fileDelta)
		Delta = zeros( (len(oms_ind), shape(Dfile)[1]), dtype=float )
		for i in range(len(oms_ind)):
			Delta[i,:] = Dfile[oms_ind[i],:]

	else:
		print "No old file."
		print 
		Delt = 0.5 *1j*(oms_f_log-sqrt(oms_f_log**2+1)) #1/4*G
		Delta = ftn_to_array(oms_f_log,Delt)
	
###############self-consistent loop################
	p_ener = 0.0
	for itt in range(Nitt):
		savetxt("Delta.dat."+str(itt),Delta)
		if itt!=0:
			savetxt("Gloc.dat."+str(itt),Glocal)
		(Sigs,Sigs_regular, n_imp, TrSigmaG) = \
				ctqmc.Run(U, mu, oms_f, oms_ind, Delta,Mstep)
		ctqmc.copyfiles(itt)

		EU = TrSigmaG # = 1/2*Tr(Sigma*G)
		Ekin = Cmp_Ekin_Sigma(Sigs_regular)
		Etotal = Ekin+EU
		print 
		print "itt= ",itt
		print "Ekin= %2.4f" %Ekin
		print "EU= %2.4f" %EU
		print "Etotal= %2.4f" %Etotal
		print "nf= %2.4f" %(n_imp/2)
		print "U= %2.4f mu= %2.4f T= %2.4f Etotal= %2.4f" %(U,mu,1/beta,Etotal)
		
		#New Delta
		Delta,Glocal = Cmp_Gloc_Delta(Sigs)

		if itt==Nitt-1:
			print 
			print "Final state"
			print "U= %2.4f mu= %2.4f T= %2.4f Etotal= %2.4f" %(U,mu,1/beta,Etotal)
