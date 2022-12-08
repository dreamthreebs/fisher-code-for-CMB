import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

#Assume installed from github using "git clone --recursive https://github.com/cmbant/CAMB.git"
#This file is then in the docs folders
camb_path = os.path.realpath(os.path.join(os.getcwd(),'..'))
sys.path.insert(0,camb_path)
import camb
from camb import model, initialpower
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))

pars=camb.read_ini(os.path.join('/sharefs/alicpt/users/wangyiming25/work/cb/CAMB','inifiles','planck_2018.ini'))
print(pars)
results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
for name in powers: print(name)
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print(totCL.shape)

ombh2_0=pars.ombh2
h=0.00001*ombh2_0
pars.set_cosmology(H0=pars.H0,ombh2=ombh2_0-h,omch2=pars.omch2)
fnresults=camb.get_results(pars)
fnpowers=fnresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fntotCL=fnpowers['total']
pars.set_cosmology(H0=pars.H0,ombh2=ombh2_0+h,omch2=pars.omch2)
fpresults=camb.get_results(pars)
fppowers=fpresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fptotCL=fppowers['total']
ombh2_CLprime=(fptotCL-fntotCL)/(2*h)
pars.set_cosmology(H0=pars.H0,ombh2=ombh2_0,omch2=pars.omch2)
#plt.plot(ls,ombh2_CLprime[:,0])
#plt.show()
print(ombh2_CLprime[500,0])

omch2_0=pars.omch2
h=0.0001*omch2_0
pars.set_cosmology(H0=pars.H0,ombh2=pars.ombh2,omch2=omch2_0-h)
fnresults=camb.get_results(pars)
fnpowers=fnresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fntotCL=fnpowers['total']
pars.set_cosmology(H0=pars.H0,ombh2=pars.ombh2,omch2=omch2_0+h)
fpresults=camb.get_results(pars)
fppowers=fpresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fptotCL=fppowers['total']
omch2_CLprime=(fptotCL-fntotCL)/(2*h)
pars.set_cosmology(H0=pars.H0,ombh2=pars.ombh2,omch2=omch2_0)
#plt.plot(ls,omch2_CLprime[:,0])
#plt.show()
print(omch2_CLprime[500,0])


As_0=pars.InitPower.As
h=0.0001*As_0
pars.InitPower.set_params(As=As_0-h, ns=pars.InitPower.ns, r=pars.InitPower.r)
fnresults=camb.get_results(pars)
fnpowers=fnresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fntotCL=fnpowers['total']
pars.InitPower.set_params(As=As_0+h, ns=pars.InitPower.ns, r=pars.InitPower.r)
fpresults=camb.get_results(pars)
fppowers=fpresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fptotCL=fppowers['total']
As_CLprime=(fptotCL-fntotCL)/(2*h)
pars.InitPower.set_params(As=As_0, ns=pars.InitPower.ns, r=pars.InitPower.r)
#plt.plot(ls,As_CLprime[:,0])
#plt.show()
print(As_CLprime[500,0])


ns_0=pars.InitPower.ns
h=0.005*ns_0
pars.InitPower.set_params(As=pars.InitPower.As, ns=ns_0-h, r=pars.InitPower.r)
fnresults=camb.get_results(pars)
fnpowers=fnresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fntotCL=fnpowers['total']
pars.InitPower.set_params(As=pars.InitPower.As, ns=ns_0+h, r=pars.InitPower.r)
fpresults=camb.get_results(pars)
fppowers=fpresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fptotCL=fppowers['total']
ns_CLprime=(fptotCL-fntotCL)/(2*h)
pars.InitPower.set_params(As=pars.InitPower.As, ns=ns_0, r=pars.InitPower.r)
#plt.plot(ls,ns_CLprime[:,0])
#plt.show()
print(ns_CLprime[500,0])


optical_depth_0=pars.Reion.optical_depth
h=0.0001*optical_depth_0
pars.Reion.set_tau(tau=optical_depth_0-h)
fnresults=camb.get_results(pars)
fnpowers=fnresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fntotCL=fnpowers['total']
pars.Reion.set_tau(tau=optical_depth_0+h)
fpresults=camb.get_results(pars)
fppowers=fpresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fptotCL=fppowers['total']
optical_depth_CLprime=(fptotCL-fntotCL)/(2*h)
pars.Reion.set_tau(tau=optical_depth_0)
#plt.plot(ls,optical_depth_CLprime[:,0])
#plt.show()
print(optical_depth_CLprime[500,0])


H0_0=pars.H0
h=0.0001*H0_0
pars.set_cosmology(H0=H0_0-h,ombh2=pars.ombh2,omch2=pars.omch2)
fnresults=camb.get_results(pars)
fnpowers=fnresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fntotCL=fnpowers['total']
pars.set_cosmology(H0=H0_0+h,ombh2=pars.ombh2,omch2=pars.omch2)
fpresults=camb.get_results(pars)
fppowers=fpresults.get_cmb_power_spectra(pars, CMB_unit='muK')
fptotCL=fppowers['total']
H0_CLprime=(fptotCL-fntotCL)/(2*h)
pars.set_cosmology(H0=H0_0,ombh2=pars.ombh2,omch2=pars.omch2)
#plt.plot(ls,H0_CLprime[:,0])
#plt.show()
print(H0_CLprime[500,0])


ls=np.arange(2601)
fig, ax = plt.subplots(3,2, figsize = (12,12))
ax[0,0].plot(ls[2:],ombh2_CLprime[2:,0])
ax[0,0].set_title('ombh2_CLprime')
ax[0,1].plot(ls[2:], omch2_CLprime[2:,0]);
ax[0,1].set_title('omch2_CLprime')
ax[1,0].plot(ls[2:], As_CLprime[2:,0]);
ax[1,0].set_title('As_CLprime')
ax[1,1].plot(ls[2:], ns_CLprime[2:,0]);
ax[1,1].set_title('ns_CLprime');
ax[2,0].plot(ls[2:], optical_depth_CLprime[2:,0]);
ax[2,0].set_title('optical_depth_CLprime');
ax[2,1].plot(ls[2:], H0_CLprime[2:,0]);
ax[2,1].set_title('H0_CLprime');
plt.show()

FM=np.zeros((6,6))
Derivative={0:ombh2_CLprime[:,0],1:omch2_CLprime[:,0],2:As_CLprime[:,0],3:ns_CLprime[:,0],4:optical_depth_CLprime[:,0],5:H0_CLprime[:,0]}
for i in np.arange(6):
	for j in np.arange(6):
		for l in np.arange(2,2601):
			FM[i,j]+=((2*l+1)*Derivative.get(i)[l]*Derivative.get(j)[l])/(2*((totCL[l,0]+0*totCL[l,0]*(np.exp(-(l**2)*(np.pi/1080)**2/2))**(-2))**2))
print(FM)
F=np.mat(FM)
FI=F.I
FIarray=np.array(FI)
sigma=np.zeros((6,1))
for i in np.arange(6):
	sigma[i]=np.sqrt(FIarray[i,i])
	print(sigma[i])
