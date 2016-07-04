from scipy import *
from scipy import linalg
from scipy import interpolate
import sys
import time


def log_mesh(om,nom, ntail_,del_nom = 1):
    """Creates logarithmic mesh on Matsubara axis
       Takes first istart points from mesh om and
       the rest of om mesh is replaced by ntail poinst
       redistribued logarithmically.
       Input:
           om      -- original long mesh
           istart  -- first istart points unchanged
           ntail   -- tail replaced by ntail points only
       Output:
           som             -- smaller mesh created from big om mesh
           sSig[nom,nc]    -- Sig on small mesh
       Also computed but not returned:
           ind_om  -- index array which conatins index to
                      kept Matsubara points
    """
    
    istart = min(nom, len(om))
    ntail = min(ntail_, len(om)-istart)

    ind_om=[]
    alpha = log((len(om)-1.)/istart)/(ntail-1.)
    for i in range(istart)[::del_nom]:
        ind_om.append(i)
    for i in range(ntail):
        t = int(istart*exp(alpha*i)+0.5)
        if (t != ind_om[-1]):
            ind_om.append(t)

    ind_oms_equal = [[0]]
    for it in range(1,len(ind_om)-1):
        istart = int(0.5*(ind_om[it-1]+ind_om[it])+0.51)
        iend = int(0.5*(ind_om[it]+ind_om[it+1])-0.01)
        equal = [i for i in range(istart,iend+1)]
        ind_oms_equal.append(equal)
    istart = int(0.5*(ind_om[-2]+ind_om[-1])+0.51)
    equal = [i for i in range(istart,ind_om[-1]+1)]
    ind_oms_equal.append(equal)

    oms_equal=[]
    for ind in ind_oms_equal:
        oms_equal.append( array([om[i] for i in ind]) )
    
    return (ind_om,oms_equal)




def Spline_Real(Arr,mesh_log,mesh_reg):
	#Arr = 1D real array
	#Arr = Arr[:10]
	#mesh_log = mesh_log[:10]
	#print Arr, mesh_log

	Re = interpolate.UnivariateSpline(mesh_log,Arr,s=0)
	Arr_large = Re(mesh_reg)
	#Nmesh = len(mesh_log)
	#print mesh_log[Nmesh/2]
	#print Re(0.05)
	#sys.exit()
	return Arr_large


def Spline_Complex(Arr,mesh_log,mesh_reg):
	#Arr = 1D complex array
	Re = interpolate.UnivariateSpline(mesh_log,real(Arr),s=0)
	Im = interpolate.UnivariateSpline(mesh_log,imag(Arr),s=0)
	Arr_large = Re(mesh_reg)+1j*Im(mesh_reg)
	return Arr_large


def Create_om_mesh(beta, Nomega, nom, nom_tail):
	oms_f = [] #fermion frequency
	for i in range(Nomega):
		oms_f.append((2*i+1)*pi/beta)
	oms_f = array(oms_f)

	oms_b = [] #boson frequency
	for i in range(Nomega):
		s = (2*i)*pi/beta
		oms_b.append(s)
	oms_b = array(oms_b)

	oms_ind, oms_equal_f = log_mesh(oms_f,nom,nom_tail)

	oms_f_log = []
	for i in oms_ind:
		oms_f_log.append((2*i+1)*pi/beta)
	oms_f_log = array(oms_f_log)

	oms_b_log = []
	for i in oms_ind:
		s = (2*i)*pi/beta
		oms_b_log.append(s)
	oms_b_log = array(oms_b_log)
	return oms_f, oms_f_log, oms_b, oms_b_log, oms_equal_f, oms_ind


if __name__ =="__name__":
	beta = 30.
	Ntau = 2000
	tau_safe = 10.
	ntau = 50
	tau_mesh, tau_large = Create_tau_mesh2(beta,Ntau,ntau,tau_safe)
