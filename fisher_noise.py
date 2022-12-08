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

#get CMB TT EE BB TE power spectrum
pars=camb.read_ini(os.path.join('/sharefs/alicpt/users/wangyiming25/work/cb/CAMB','inifiles','planck_2018.ini'))
print(pars)
results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')
for name in powers: print(name)
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print(totCL.shape)

#Derivative of partial Cl partial ombh2
ls=np.arange(2601)
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
print('ombh2 is ok')

#Derivative of partial Cl partial omch2
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
print('omch2 is ok')

#Derivative of partial Cl partial As
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
print('As is ok')

#Derivative of partial Cl partial ns
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
print('ns is ok')

#Derivative of partial Cl partial tau
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
print('tau is ok')

#Derivative of partial Cl partial H0
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
print('H0 is ok')

#plot the derivative of CL/parameters
ls=np.arange(2601)
fig, ax = plt.subplots(3,2, figsize = (12,12))
ax[0,0].plot(ls[2:],ombh2_CLprime[2:,0])
ax[0,0].set_title('ombh2_CLTTprime')
ax[0,1].plot(ls[2:], omch2_CLprime[2:,0]);
ax[0,1].set_title('omch2_CLTTprime')
ax[1,0].plot(ls[2:], As_CLprime[2:,0]);
ax[1,0].set_title('As_CLTTprime')
ax[1,1].plot(ls[2:], ns_CLprime[2:,0]);
ax[1,1].set_title('ns_CLTTprime');
ax[2,0].plot(ls[2:], optical_depth_CLprime[2:,0]);
ax[2,0].set_title('optical_depth_CLTTprime');
ax[2,1].plot(ls[2:], H0_CLprime[2:,0]);
ax[2,1].set_title('H0_CLTTprime');
plt.show()


fig, ax = plt.subplots(3,2, figsize = (12,12))
ax[0,0].plot(ls[2:],ombh2_CLprime[2:,1])
ax[0,0].set_title('ombh2_CLTEprime')
ax[0,1].plot(ls[2:], omch2_CLprime[2:,1]);
ax[0,1].set_title('omch2_CLTEprime')
ax[1,0].plot(ls[2:], As_CLprime[2:,1]);
ax[1,0].set_title('As_CLTEprime')
ax[1,1].plot(ls[2:], ns_CLprime[2:,1]);
ax[1,1].set_title('ns_CLTEprime');
ax[2,0].plot(ls[2:], optical_depth_CLprime[2:,1]);
ax[2,0].set_title('optical_depth_CLTEprime');
ax[2,1].plot(ls[2:], H0_CLprime[2:,1]);
ax[2,1].set_title('H0_CLTEprime');
plt.show()

fig, ax = plt.subplots(3,2, figsize = (12,12))
ax[0,0].plot(ls[2:],ombh2_CLprime[2:,2])
ax[0,0].set_title('ombh2_CLEEprime')
ax[0,1].plot(ls[2:], omch2_CLprime[2:,2]);
ax[0,1].set_title('omch2_CLEEprime')
ax[1,0].plot(ls[2:], As_CLprime[2:,2]);
ax[1,0].set_title('As_CLEEprime')
ax[1,1].plot(ls[2:], ns_CLprime[2:,2]);
ax[1,1].set_title('ns_CLEEprime');
ax[2,0].plot(ls[2:], optical_depth_CLprime[2:,2]);
ax[2,0].set_title('optical_depth_CLEEprime');
ax[2,1].plot(ls[2:], H0_CLprime[2:,2]);
ax[2,1].set_title('H0_CLEEprime');
plt.show()

fig, ax = plt.subplots(3,2, figsize = (12,12))
ax[0,0].plot(ls[2:],ombh2_CLprime[2:,3])
ax[0,0].set_title('ombh2_CLBBprime')
ax[0,1].plot(ls[2:], omch2_CLprime[2:,3]);
ax[0,1].set_title('omch2_CLBBprime')
ax[1,0].plot(ls[2:], As_CLprime[2:,3]);
ax[1,0].set_title('As_CLBBprime')
ax[1,1].plot(ls[2:], ns_CLprime[2:,3]);
ax[1,1].set_title('ns_CLBBprime');
ax[2,0].plot(ls[2:], optical_depth_CLprime[2:,3]);
ax[2,0].set_title('optical_depth_CLprime');
ax[2,1].plot(ls[2:], H0_CLprime[2:,3]);
ax[2,1].set_title('H0_CLBBprime');
plt.show()
#Noise powerspectrum
NET=350# muk sqrt(s)
N_det=68158
T_obs=6311520# s
f_sky=0.1
Y=0.8
Aperture=0.72# m
frequency_95=95e9# Hz
frequency_150=150e9# Hz
theta_FWHM_95=1.02*3e8/(frequency_95*Aperture)# rad
theta_FWHM_150=1.02*3e8/(frequency_150*Aperture)# rad
sigma2=4*np.pi*f_sky*NET**2/(T_obs*N_det*Y)
Noise=np.ones(2601)
for l in np.arange(2601):
	Noise[l]=sigma2*np.exp(l*(l+1)*theta_FWHM_95**2/(8*np.log(2)))
plt.plot(ls,Noise)
plt.show()

#get fisher matrix
FM=np.zeros((6,6))
Cell=np.zeros((3,3))
Cellprime_i=np.zeros((3,3))
Cellprime_j=np.zeros((3,3))
DerivativeTT={0:ombh2_CLprime[:,0],1:omch2_CLprime[:,0],2:As_CLprime[:,0],3:ns_CLprime[:,0],4:optical_depth_CLprime[:,0],5:H0_CLprime[:,0]}
DerivativeEE={0:ombh2_CLprime[:,1],1:omch2_CLprime[:,1],2:As_CLprime[:,1],3:ns_CLprime[:,1],4:optical_depth_CLprime[:,1],5:H0_CLprime[:,1]}
DerivativeBB={0:ombh2_CLprime[:,2],1:omch2_CLprime[:,2],2:As_CLprime[:,2],3:ns_CLprime[:,2],4:optical_depth_CLprime[:,2],5:H0_CLprime[:,2]}
DerivativeTE={0:ombh2_CLprime[:,3],1:omch2_CLprime[:,3],2:As_CLprime[:,3],3:ns_CLprime[:,3],4:optical_depth_CLprime[:,3],5:H0_CLprime[:,3]}
for i in np.arange(6):
	for j in np.arange(6):
		for l in np.arange(2,2601):
			Cell[0,0]=totCL[l,0]+Noise[l]
			Cell[1,0]=totCL[l,3]
			Cell[0,1]=totCL[l,3]
			Cell[1,1]=totCL[l,1]+Noise[l]
			Cell[2,2]=totCL[l,2]+Noise[l]
			Cellm=np.mat(Cell)
			Cellmatrix_inv=Cellm.I
			Cellprime_i[0,0]=DerivativeTT.get(i)[l]
			Cellprime_i[1,0]=DerivativeTE.get(i)[l]
			Cellprime_i[0,1]=DerivativeTE.get(i)[l]
			Cellprime_i[1,1]=DerivativeEE.get(i)[l]
			Cellprime_i[2,2]=DerivativeBB.get(i)[l]
			Cellprime_i_matrix=np.mat(Cellprime_i)
			Cellprime_j[0,0]=DerivativeTT.get(j)[l]
			Cellprime_j[1,0]=DerivativeTE.get(j)[l]
			Cellprime_j[0,1]=DerivativeTE.get(j)[l]
			Cellprime_j[1,1]=DerivativeEE.get(j)[l]
			Cellprime_j[2,2]=DerivativeBB.get(j)[l]
			Cellprime_j_matrix=np.mat(Cellprime_j)
			Mul=np.dot(Cellmatrix_inv,Cellprime_i_matrix)
			Mult=np.dot(Mul,Cellmatrix_inv)
			Multi=np.dot(Mult,Cellprime_j_matrix)
			FM[i,j]+=(2*l+1)*Multi.trace()/2
			print(l,end=' ')
			if l==2016:
				print(i,j)
		print(FM)
print(FM)
F=np.mat(FM)
FI=F.I
FIarray=np.array(FI)
sigma=np.zeros((6,1))
#get covariance
for i in np.arange(6):
	sigma[i]=np.sqrt(FIarray[i,i])
	print(sigma[i])
