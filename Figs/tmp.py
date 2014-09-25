import numpy as np
import scipy.special as SS
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

fcst=open("../Constants.jl","r")
constants={}
for l in fcst:
	line=l.split()
	if line and "const" in line[0]:
		try:
			constants[line[1]]=eval(line[3].replace(';',''))
		except Exception as e:
			print e

def tausr_base(PC_rp,PC_f,PC_q,PF_sig ,PF_n,GRD_nx,GRD_dx,PC_v):
	v=PC_v
	a=PC_f#+PF_sig*.7#*sqrt(2*log(PF_n*sqrt(2*pi*PF_sig**2))) /2
	b=GRD_nx*GRD_dx/2.
	t2=(1.-PC_rp)/PC_rp #avg time of straight flight
	t1=PC_rp/(1.-PC_rp) #avg time of wait 
	t2opt= (a/v * sqrt( log(b/a) -1./2.))
	popt=1/(t2opt+1)
	print 1.-popt
	return v,a,b,t1,t2
	
def tausr1(PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v):  #Static mode
	v,a,b,t1,t2=tausr_base(PC_rp,PC_f,PC_q,PF_sig ,PF_n ,GRD_nx,GRD_dx,PC_v)
	k=10000000./PC_q #rate of catching
	xy=sqrt( 2*k*t1/(1+k*t1) )/v/t2
	x=a*xy
	y=b*xy
	result1=xy**1/a*(1+k*t1)*(b**2-a**2)**2*SS.iv(0,x)/SS.iv(1,x)
	result2=  8*b**2 + xy**2 * (1+k*t1)*(4*b**4*log(b/a) + (b**2-a**2)*(a**2-3*b**2+8/xy**2))
	#result2=1
	result3= result1+ result2/4.
	result4= (t1+t2)/(2*k*t1*b**2) * result3
	#return (t1+t2)/(2*k*t1*y**2)*result2/4.
	return result4

def tausr2(*args):  #Diffusive mode
	v,a,b,t1,t2=tausr_base(*args)
	D=1.
	Db=v**2*t2/2
	alp=sqrt(1./(D*t1) + 1./(Db*t2) )
	elem1=SS.iv(1,b*alp) * SS.kv(1,a*alp)-SS.iv(1,a*alp) * SS.kv(1,b*alp)
	elem2=SS.iv(1,b*alp) * SS.kv(0,a*alp)+SS.iv(0,a*alp) * SS.kv(1,b*alp)
	Lelem1=SS.iv(0,a/sqrt(Db*t2) )* elem1
	Lelem2= alp * sqrt(Db*t2) *SS.iv(1,a/sqrt(Db*t2) )* elem2
	Lp=Lelem1 + Lelem2
	Lm= Lelem1-Lelem2
	M= SS.iv(0,a/sqrt(Db*t2) )* elem2  -4 * a**2* sqrt(Db*t2) /(alp* (b**2-a**2)**2) *SS.iv(1,a/sqrt(Db*t2) )* elem1
	result1=a*alp*(b**2/a**2-1) * M/2/Lp - Lm / Lp
	result2= (3-4*log(b/a))*b**4 - 4*a**2 *b**2 + a**4
	result3= result1 - a**2 * D *t1/ (8*Db*t2) * result2 / (b**2-a**2)
	result4=(t1+t2) * (1-a**2/b**2)/(alp**2 * D *t1)**2 * result3
	return result4

def tausr3(*args):  #Minimal approximation
	v,a,b,t1,t2=tausr_base(*args)
	x= a * sqrt(2.)/v/t2
	result1= (b**2-a**2)**2/ (sqrt(2) *v*a*b**2 ) *  SS.iv(0,x)/SS.iv(1,x)
	result2= b**4 * log(b/a) + (b**2-a**2) * (a**2 - 3 *b**2 + 4 * v**2 *t2**2)/4.
	return   result1 + result2/(v**2 * t2 * b**2)

args=[constants[i] for i in ("PC_rp","PC_f","PC_q","PF_sig" ,"PF_n" ,"GRD_nx","GRD_dx","PC_v")]

## Load in data
TS = np.load("../Data/Data_Fig2a.npy")
#print(TS)
xs=np.linspace(0.1,.95,len(TS))
plt.plot(xs,TS)
plt.plot(xs,[tausr1(1-x,*args[1:]) for x in xs])
plt.show()


TS = np.load("../Data/Data_Fig2b.npy")
print TS.shape
from itertools import product
rg=np.linspace(0.,.5,TS.shape[0])
X,Y = np.meshgrid(rg,rg)
xy=np.array([x for x in product(rg,rg)])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.xlabel("F_sigma")
plt.ylabel("C_f")
#plt.zlabel("tau_s")
ax.set_zlim(bottom=0, top=500)
ax.scatter(xy.T[0],xy.T[1],TS[:,:,0].ravel(),c=TS[:,:,0].ravel())
plt.show()
Z=X.copy()
for i in range(Z.shape[0]):
	for j in range(Z.shape[1]):
		Z[i,j]=TS[i,j][0]
fig = plt.figure()
plt.xlabel("F_sigma")
plt.ylabel("C_f")
ax = fig.add_subplot(111)
ax.pcolor(X, Y, Z,  vmin=abs(Z).min(), vmax=500)#abs(Z).max())
plt.show()
