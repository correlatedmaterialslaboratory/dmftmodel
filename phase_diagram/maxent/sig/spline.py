from scipy import *
from scipy import linalg
from scipy import interpolate
import sys
import time

in_sig=str(sys.argv[1])
out_sig=str(sys.argv[2])

def Spline_Real(Arr,mesh_log,mesh_reg):
	Re = interpolate.UnivariateSpline(mesh_log,Arr,s=0)
	Arr_large = Re(mesh_reg)
	return Arr_large

beta=100
Nomega=5000

data=loadtxt(in_sig).T
s_oo=data[1][-1]
oms_s = data[0]
oms=pi/beta*(2*arange(Nomega)+1)
sig_r = Spline_Real(data[1]-s_oo,oms_s,oms)
sig_i = Spline_Real(data[2],oms_s,oms)
savetxt(out_sig,array([oms,sig_r,sig_i]).T)

