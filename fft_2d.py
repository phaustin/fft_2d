from __future__ import division
from netCDF4 import Dataset
import numpy as np
import math
from scipy import fftpack
from matplotlib import pyplot as plt
#plt.switch_backend('Agg') #batch
plt.switch_backend('MacOSX') #interactive

 
test=Dataset('a17.nc')
tau=test.variables['tau'][:2000,:2000]


nn = tau.shape[0]     # size of each column of the 2D array
mm = tau.shape[1]     # size of each row of the array
m = int(math.floor(mm/2))  #midpoint

scale=0.025  #pixel size in km
x_dens = np.arange(0,(m))+1
x_dens = (x_dens-m)/(mm*scale)
delta_k = 1./scale #1/km
nyquist = delta_k*0.5

fft_tau = fftpack.fft2(tau)
tr_tau = fftpack.fftshift(fft_tau)
e_dens = tr_tau*np.conjugate(tr_tau)/(mm*mm)
e_dens = e_dens.real


plt.close('all')
fig,ax=plt.subplots(2,2)
ax[0,0].set_title('title1')
im1=ax[0,0].imshow(tau)
im2=ax[1,0].imshow(np.log(e_dens))
im3=ax[0,1].hist(tau.ravel())
im4=ax[1,1].hist(np.log(e_dens.ravel()))
plt.draw()
cbar_ax = fig.add_axes([0.45, 0.55, 0.03, 0.3])
fig.colorbar(im1,cax=cbar_ax)
fig.tight_layout()
fig.canvas.draw()
plt.show()

bnstep=2.
nbns = int(round((math.sqrt(2)*mm/bnstep),0)+1)

e_spec = np.zeros(nbns,np.float)
cnt = np.zeros(nbns,np.float)

for i in range(mm):
    if (i%100) == 0:
        print "\t\trow: "+str(i)+" completed"
    for j in range(mm):
        r = math.sqrt(((i+1)-mm/2)**2+((j+1)-mm/2)**2)
        bn = int(math.floor(r/bnstep))
        e_spec[bn]=e_spec[bn]+ np.abs(e_dens[i,j])**2.
        cnt[bn]=cnt[bn]+1

for i in range(nbns):
    if cnt[i]>0:
        e_spec[i]=e_spec[i]/cnt[i]/(4*(math.pi**2))

e_spec=np.sqrt(e_spec)

delta_k=nyquist/nbns
x_ax=np.linspace(delta_k,nyquist,nbns)
        
fig=plt.figure(2)
fig.clf()
ax1=fig.add_subplot(111)
ax1.loglog(x_ax,e_spec)
l0=1.
slope=(-8/3.)
analytic=l0*x_ax**slope
the_line=l0
ax1.loglog(x_ax,analytic,'r-')
fig.tight_layout()
fig.canvas.draw()
plt.show()

        
 
