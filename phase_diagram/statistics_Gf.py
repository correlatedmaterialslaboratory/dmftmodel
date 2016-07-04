
from scipy import *
from scipy import interpolate
import sys
import os
from pylab import *


def Spline_Real(Arr,mesh_log,mesh_reg):
	Re = interpolate.UnivariateSpline(mesh_log,Arr,s=0)
	Arr_large = Re(mesh_reg)
	return Arr_large

def Take_im0(GfFilename):
	#take imG(iw=0)
	x0 = 0.
	data = loadtxt(GfFilename).transpose()
	
	Xx = data[0]
	Yy = data[2]
	Yspline = interpolate.UnivariateSpline(Xx,Yy,s=0)
	y0 = Yspline(x0)
	return y0
	
itt_min=0
itt_max=200
which = range(itt_min,itt_max)

which_= []
Sig = []
Gf = []
for itt in which:
	GfFile = 'Gf.out.'+str(itt)
	if os.path.exists(GfFile):
		ImGw0 = Take_im0(GfFile)
		print "%d, %4.6f" %(itt,ImGw0)
		Gf.append(-ImGw0/pi)
		which_.append(itt)


if len(sys.argv)>1:
	plot(which_,Gf)
	show()
	
