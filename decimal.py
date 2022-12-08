import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
camb_path = os.path.realpath(os.path.join(os.getcwd(),'..'))
sys.path.insert(0,camb_path)
import camb
from camb import model, initialpower
from decimal import Decimal,getcontext 
getcontext().prec = 15
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))
pars = camb.CAMBparams()

#set parameters
pars.set_cosmology(ombh2=0.022032,omch2=0.12038,thetastar=0.0104119,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
pars.InitPower.set_params(ns=0.9619,As=2.21536e-9,nrun=0)
pars.set_for_lmax(lmax=2500,max_eta_k=10000)
#get results TT EE BB TE from number 0 to 3 
results = camb.get_results(pars)
powers=results.get_cmb_power_spectra(pars,lmax=2500,CMB_unit='muK',raw_cl=True)
for name in powers: print(name)
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print(totCL.shape)
print(totCL[100,0])
ls=np.arange(totCL.shape[0])

#get thetastar(one of H0, thetastar, cosmotheta_mc must be given in set_cosmology function)
derived_parameters=results.get_derived_params()
thetastar_0=derived_parameters['thetastar']/100
#get initial value of 6 standard parameters(except thetastar)
ombh2_0=pars.ombh2
omch2_0=pars.omch2
optical_depth_0=pars.Reion.optical_depth
ns_0=pars.InitPower.ns
As_0=pars.InitPower.As
h=0.001
footstep=0.001


#using CubicSpline method to get derivative
#Derivative of partial Cl partial ombh2
x1=ombh2_0-2*h*ombh2_0
x2=ombh2_0-h*ombh2_0
x3=ombh2_0+h*ombh2_0
x4=ombh2_0+2*h*ombh2_0
ombh2_0=Decimal(str(pars.ombh2))
pars.set_cosmology(thetastar=thetastar_0,omch2=omch2_0,ombh2=x1,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y1_results=camb.get_results(pars)
y1_powers=y1_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y1_totCL=y1_powers['total']
pars.set_cosmology(thetastar=thetastar_0,omch2=omch2_0,ombh2=x2,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y2_results=camb.get_results(pars)
y2_powers=y2_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y2_totCL=y2_powers['total']
pars.set_cosmology(thetastar=thetastar_0,omch2=omch2_0,ombh2=x3,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y3_results=camb.get_results(pars)
y3_powers=y3_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y3_totCL=y3_powers['total']
pars.set_cosmology(thetastar=thetastar_0,omch2=omch2_0,ombh2=x4,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y4_results=camb.get_results(pars)
y4_powers=y4_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y4_totCL=y4_powers['total']
x=np.array([Decimal(str(x1)),Decimal(str(x2)),Decimal(str(x3)),Decimal(str(x4))])
fn=np.zeros((2501))
fp=np.zeros((2501))
ombh2_CLprime=np.zeros((2501,4))
for i in np.arange(4):#i is powerspectrum index
    CubicSpline_function_ombh2=[]#at different l
    for l in np.arange(30,2501):
        y=[Decimal(str(y1_totCL[l,i])),Decimal(str(y2_totCL[l,i])),Decimal(str(y3_totCL[l,i])),Decimal(str(y4_totCL[l,i]))]
        cs=CubicSpline(x,y)
        footstep=Decimal('0.001')
        fn[l]=cs((Decimal('1')-footstep)*ombh2_0)
        fp[l]=cs((Decimal('1')+footstep)*ombh2_0)
        ombh2_CLprime[l,i]=float((Decimal(str(fp[l]))-Decimal(str(fn[l])))/(Decimal('2')*footstep))
pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
print('ombh2 is ok')
footstep=0.001
#Derivative of partial Cl partial omch2
x1=(1-2*h)*omch2_0
x2=(1-h)*omch2_0
x3=(1+h)*omch2_0
x4=(1+2*h)*omch2_0
pars.set_cosmology(thetastar=thetastar_0,omch2=x1,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y1_results=camb.get_results(pars)
y1_powers=y1_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y1_totCL=y1_powers['total']
pars.set_cosmology(thetastar=thetastar_0,omch2=x2,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y2_results=camb.get_results(pars)
y2_powers=y2_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y2_totCL=y2_powers['total']
pars.set_cosmology(thetastar=thetastar_0,omch2=x3,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y3_results=camb.get_results(pars)
y3_powers=y3_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y3_totCL=y3_powers['total']
pars.set_cosmology(thetastar=thetastar_0,omch2=x4,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y4_results=camb.get_results(pars)
y4_powers=y4_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y4_totCL=y4_powers['total']
x=np.array([x1,x2,x3,x4])
fn=np.zeros((2501))
fp=np.zeros((2501))
omch2_CLprime=np.zeros((2501,4))
for i in np.arange(4):#i is powerspectrum index
    CubicSpline_function_omch2=[]#at different l
    for l in np.arange(30,2501):
        y=[y1_totCL[l,i],y2_totCL[l,i],y3_totCL[l,i],y4_totCL[l,i]]
        cs=CubicSpline(x,y,bc_type='natural')
        fn[l]=cs((1-footstep)*omch2_0)
        fp[l]=cs((1+footstep)*omch2_0)
    omch2_CLprime[:,i]=(fp-fn)/(2*footstep)
pars.set_cosmology(thetastar=thetastar_0,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
print('omch2 is ok')

#Derivative of partial Cl partial thetastar
x1=(1-2*h)*thetastar_0
x2=(1-h)*thetastar_0
x3=(1+h)*thetastar_0
x4=(1+2*h)*thetastar_0
pars.set_cosmology(thetastar=x1,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y1_results=camb.get_results(pars)
y1_powers=y1_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y1_totCL=y1_powers['total']
pars.set_cosmology(thetastar=x2,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y2_results=camb.get_results(pars)
y2_powers=y2_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y2_totCL=y2_powers['total']
pars.set_cosmology(thetastar=x3,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y3_results=camb.get_results(pars)
y3_powers=y3_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y3_totCL=y3_powers['total']
pars.set_cosmology(thetastar=x4,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
y4_results=camb.get_results(pars)
y4_powers=y4_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y4_totCL=y4_powers['total']
x=[x1,x2,x3,x4]
fn=np.zeros((2501))
fp=np.zeros((2501))
thetastar_CLprime=np.zeros((2501,4))
for i in np.arange(4):#i is powerspectrum index
    CubicSpline_function_thetastar=[]#at different l
    for l in np.arange(30,2501):
        y=[y1_totCL[l,i],y2_totCL[l,i],y3_totCL[l,i],y4_totCL[l,i]]
        CubicSpline_function_thetastar.append(CubicSpline(x,y,bc_type='natural'))
        fn[l]=CubicSpline_function_thetastar[l-30]((1-footstep)*thetastar_0)
        fp[l]=CubicSpline_function_thetastar[l-30]((1+footstep)*thetastar_0)
    thetastar_CLprime[:,i]=(fp-fn)/(2*footstep)
pars.set_cosmology(thetastar=thetastar_0,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
print('thetastar is ok')

#Derivative of partial Cl partial tau
x1=(1-2*h)*optical_depth_0
x2=(1-h)*optical_depth_0
x3=(1+h)*optical_depth_0
x4=(1+2*h)*optical_depth_0
pars.Reion.set_tau(tau=x1)
y1_results=camb.get_results(pars)
y1_powers=y1_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y1_totCL=y1_powers['total']
pars.Reion.set_tau(tau=x2)
y2_results=camb.get_results(pars)
y2_powers=y2_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y2_totCL=y2_powers['total']
pars.Reion.set_tau(tau=x3)
y3_results=camb.get_results(pars)
y3_powers=y3_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y3_totCL=y3_powers['total']
pars.Reion.set_tau(tau=x4)
y4_results=camb.get_results(pars)
y4_powers=y4_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y4_totCL=y4_powers['total']
x=[x1,x2,x3,x4]
fn=np.zeros((2501))
fp=np.zeros((2501))
optical_depth_CLprime=np.zeros((2501,4))
for i in np.arange(4):#i is powerspectrum index
    CubicSpline_function_optical_depth=[]#at different l
    for l in np.arange(30,2501):
        y=[y1_totCL[l,i],y2_totCL[l,i],y3_totCL[l,i],y4_totCL[l,i]]
        CubicSpline_function_optical_depth.append(CubicSpline(x,y,bc_type='natural'))
        fn[l]=CubicSpline_function_optical_depth[l-30]((1-footstep)*optical_depth_0)
        fp[l]=CubicSpline_function_optical_depth[l-30]((1+footstep)*optical_depth_0)
    optical_depth_CLprime[:,i]=(fp-fn)/(2*footstep)
pars.Reion.set_tau(tau=optical_depth_0)
print('tau is ok')

#Derivative of partial Cl partial ns
x1=(1-2*h)*ns_0
x2=(1-h)*ns_0
x3=(1+h)*ns_0
x4=(1+2*h)*ns_0
pars.InitPower.set_params(As=As_0,ns=x1)
y1_results=camb.get_results(pars)
y1_powers=y1_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y1_totCL=y1_powers['total']
pars.InitPower.set_params(As=As_0,ns=x2)
y2_results=camb.get_results(pars)
y2_powers=y2_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y2_totCL=y2_powers['total']
pars.InitPower.set_params(As=As_0,ns=x3)
y3_results=camb.get_results(pars)
y3_powers=y3_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y3_totCL=y3_powers['total']
pars.InitPower.set_params(As=As_0,ns=x4)
y4_results=camb.get_results(pars)
y4_powers=y4_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y4_totCL=y4_powers['total']
x=[x1,x2,x3,x4]
fn=np.zeros((2501))
fp=np.zeros((2501))
ns_CLprime=np.zeros((2501,4))
for i in np.arange(4):#i is powerspectrum index
    CubicSpline_function_ns=[]#at different l
    for l in np.arange(30,2501):
        y=[y1_totCL[l,i],y2_totCL[l,i],y3_totCL[l,i],y4_totCL[l,i]]
        CubicSpline_function_ns.append(CubicSpline(x,y,bc_type='natural'))
        fn[l]=CubicSpline_function_ns[l-30]((1-footstep)*ns_0)
        fp[l]=CubicSpline_function_ns[l-30]((1+footstep)*ns_0)
    ns_CLprime[:,i]=(fp-fn)/(2*footstep)
pars.InitPower.set_params(As=As_0,ns=ns_0)
print('ns is ok')

#Derivative of partial Cl partial As
x1=(1-2*h)*np.log(1e10*As_0)
x2=(1-h)*np.log(1e10*As_0)
x3=(1+h)*np.log(1e10*As_0)
x4=(1+2*h)*np.log(1e10*As_0)
pars.InitPower.set_params(As=np.exp(x1)*1e-10,ns=ns_0)
y1_results=camb.get_results(pars)
y1_powers=y1_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y1_totCL=y1_powers['total']
pars.InitPower.set_params(As=np.exp(x2)*1e-10,ns=ns_0)
y2_results=camb.get_results(pars)
y2_powers=y2_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y2_totCL=y2_powers['total']
pars.InitPower.set_params(As=np.exp(x3)*1e-10,ns=ns_0)
y3_results=camb.get_results(pars)
y3_powers=y3_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y3_totCL=y3_powers['total']
pars.InitPower.set_params(As=np.exp(x4)*1e-10,ns=ns_0)
y4_results=camb.get_results(pars)
y4_powers=y4_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
y4_totCL=y4_powers['total']
x=[x1,x2,x3,x4]
fn=np.zeros((2501))
fp=np.zeros((2501))
As_CLprime=np.zeros((2501,4))
for i in np.arange(4):#i is powerspectrum index
    CubicSpline_function_As=[]#at different l
    for l in np.arange(30,2501):
        y=[y1_totCL[l,i],y2_totCL[l,i],y3_totCL[l,i],y4_totCL[l,i]]
        CubicSpline_function_As.append(CubicSpline(x,y,bc_type='natural'))
        fn[l]=CubicSpline_function_As[l-30]((1-footstep)*np.log(1e10*As_0))
        fp[l]=CubicSpline_function_As[l-30]((1+footstep)*np.log(1e10*As_0))
    As_CLprime[:,i]=(fp-fn)/(2*footstep)
pars.InitPower.set_params(As=As_0,ns=ns_0)
print('As is ok')

#plot the derivative of CL/parameters
#fig, ax = plt.subplots(3,2, figsize = (12,12))
#ax[0,0].loglog(ls,np.absolute(ls*ls*ombh2_CLprime[:,0]))
#ax[0,0].set_title('ombh2_CLTTprime')
#ax[0,1].loglog(ls,np.absolute( ls*ls*omch2_CLprime[:,0]));
#ax[0,1].set_title('omch2_CLTTprime')
#ax[1,0].loglog(ls,np.absolute( ls*ls*As_CLprime[:,0]));
#ax[1,0].set_title('As_CLTTprime')
#ax[1,1].loglog(ls,np.absolute(ls*ls*ns_CLprime[:,0]));
#ax[1,1].set_title('ns_CLTTprime');
#ax[2,0].loglog(ls,np.absolute(ls*ls*optical_depth_CLprime[:,0]));
#ax[2,0].set_title('optical_depth_CLTTprime');
#ax[2,1].loglog(ls,np.absolute(ls*ls*thetastar_CLprime[:,0]));
#ax[2,1].set_title('thetastar_CLTTprime');
#plt.show()
#fig, ax = plt.subplots(3,2, figsize = (12,12))
#ax[0,0].loglog(ls,np.absolute(ls*ls*ombh2_CLprime[:,1]))
#ax[0,0].set_title('ombh2_CLEEprime')
#ax[0,1].loglog(ls,np.absolute( ls*ls*omch2_CLprime[:,1]));
#ax[0,1].set_title('omch2_CLEEprime')
#ax[1,0].loglog(ls,np.absolute( ls*ls*As_CLprime[:,1]));
#ax[1,0].set_title('As_CLEEprime')
#ax[1,1].loglog(ls,np.absolute(ls*ls*ns_CLprime[:,1]));
#ax[1,1].set_title('ns_CLEEprime');
#ax[2,0].loglog(ls,np.absolute(ls*ls*optical_depth_CLprime[:,1]));
#ax[2,0].set_title('optical_depth_CLEEprime');
#ax[2,1].loglog(ls,np.absolute(ls*ls*thetastar_CLprime[:,1]));
#ax[2,1].set_title('thetastar_CLEEprime');
#plt.show()
#fig, ax = plt.subplots(3,2, figsize = (12,12))
#ax[0,0].loglog(ls,np.absolute(ls*ls*ombh2_CLprime[:,3]))
#ax[0,0].set_title('ombh2_CLTEprime')
#ax[0,1].loglog(ls,np.absolute( ls*ls*omch2_CLprime[:,3]));
#ax[0,1].set_title('omch2_CLTEprime')
#ax[1,0].loglog(ls,np.absolute( ls*ls*As_CLprime[:,3]));
#ax[1,0].set_title('As_CLTEprime')
#ax[1,1].loglog(ls,np.absolute(ls*ls*ns_CLprime[:,3]));
#ax[1,1].set_title('ns_CLTEprime');
#ax[2,0].loglog(ls,np.absolute(ls*ls*optical_depth_CLprime[:,3]));
#ax[2,0].set_title('optical_depth_CLTEprime');
#ax[2,1].loglog(ls,np.absolute(ls*ls*thetastar_CLprime[:,3]));
#ax[2,1].set_title('thetastar_CLTEprime');
#plt.show()

#Noise powerspectrum
#NET=350# muk sqrt(s)
#N_det=68158
#T_obs=6311520# s
#f_sky=0.1
#Y=0.8
#Aperture=0.72# m
#frequency_95=95e9# Hz
#frequency_150=150e9# Hz
#theta_FWHM_95=1.02*3e8/(frequency_95*Aperture)# rad
#theta_FWHM_150=1.02*3e8/(frequency_150*Aperture)# rad
#sigma2=4*np.pi*f_sky*NET**2/(T_obs*N_det*Y)
#Noise=np.ones(2601)
#for l in np.arange(2601):
#    Noise[l]=sigma2*np.exp(l*(l+1)*theta_FWHM_95**2/(8*np.log(2)))
#plt.plot(ls,Noise)
#plt.show()

#get derivative
DerivativeTT={0:ombh2_CLprime[:,0],1:omch2_CLprime[:,0],2:thetastar_CLprime[:,0],3:optical_depth_CLprime[:,0],4:ns_CLprime[:,0],5:As_CLprime[:,0]}
DerivativeEE={0:ombh2_CLprime[:,1],1:omch2_CLprime[:,1],2:thetastar_CLprime[:,1],3:optical_depth_CLprime[:,1],4:ns_CLprime[:,1],5:As_CLprime[:,1]}
DerivativeBB={0:ombh2_CLprime[:,2],1:omch2_CLprime[:,2],2:thetastar_CLprime[:,2],3:optical_depth_CLprime[:,2],4:ns_CLprime[:,2],5:As_CLprime[:,2]}
DerivativeTE={0:ombh2_CLprime[:,3],1:omch2_CLprime[:,3],2:thetastar_CLprime[:,3],3:optical_depth_CLprime[:,3],4:ns_CLprime[:,3],5:As_CLprime[:,3]}

#TT fisher matrix
FM=np.zeros((6,6))
Cell=np.zeros((1,1))
Cellprime_i=np.zeros((1,1))
Cellprime_j=np.zeros((1,1))
for i in np.arange(6):
    for j in np.arange(6):
        for l in np.arange(30,2501):
            Cell[0,0]=totCL[l,0]
            Cellm=np.mat(Cell)
            Cellmatrix_inv=Cellm.I
            Cellprime_i[0,0]=DerivativeTT.get(i)[l]
            Cellprime_i_matrix=np.mat(Cellprime_i)
            Cellprime_j[0,0]=DerivativeTT.get(j)[l]
            Cellprime_j_matrix=np.mat(Cellprime_j)
            Mul=np.dot(Cellmatrix_inv,Cellprime_i_matrix)
            Mult=np.dot(Mul,Cellmatrix_inv)
            Multi=np.dot(Mult,Cellprime_j_matrix)
            FM[i,j]+=(2*l+1)*Multi.trace()/2
            if l==2016:
                print(i,j)
        if i==j==3:
            FM[i,j]+=1/(0.013**2)
print(FM)
F=np.mat(FM)
FI=F.I
FIarray=np.array(FI)
sigma=np.zeros((6,1))
#get covariance
for i in np.arange(6):
    sigma[i]=np.sqrt(FIarray[i,i])
    print(sigma[i])
#get conditional errors
conditionalerror=np.zeros((6,1))
for i in np.arange(6):
    conditionalerror[i]=1/(np.sqrt(FM[i,i]))
    print(conditionalerror[i])

#EE fisher matrix
FM_EE=np.zeros((6,6))
for i in np.arange(6):
    for j in np.arange(6):
        for l in np.arange(30,2500):
            Cell[0,0]=totCL[l,1]
            Cellm=np.mat(Cell)
            Cellmatrix_inv=Cellm.I
            Cellprime_i[0,0]=DerivativeEE.get(i)[l]
            Cellprime_i_matrix=np.mat(Cellprime_i)
            Cellprime_j[0,0]=DerivativeEE.get(j)[l]
            Cellprime_j_matrix=np.mat(Cellprime_j)
            Mul=np.dot(Cellmatrix_inv,Cellprime_i_matrix)
            Mult=np.dot(Mul,Cellmatrix_inv)
            Multi=np.dot(Mult,Cellprime_j_matrix)
            FM_EE[i,j]+=(2*l+1)*Multi.trace()/2
        if i==j==3:
            FM_EE+=1/(0.013**2)
print(FM_EE)
F=np.mat(FM_EE)
FI=F.I
FIarray=np.array(FI)
sigma=np.zeros((6,1))
#get covariance
for i in np.arange(6):
    sigma[i]=np.sqrt(FIarray[i,i])
    print(sigma[i])

#TE fisher matrix
FM_TE=np.zeros((6,6))
for i in np.arange(6):
    for j in np.arange(6):
        for l in np.arange(30,2500):
            Cell[0,0]=totCL[l,3]
            Cellm=np.mat(Cell)
            Cellmatrix_inv=Cellm.I
            Cellprime_i[0,0]=DerivativeTE.get(i)[l]
            Cellprime_i_matrix=np.mat(Cellprime_i)
            Cellprime_j[0,0]=DerivativeTE.get(j)[l]
            Cellprime_j_matrix=np.mat(Cellprime_j)
            Mul=np.dot(Cellmatrix_inv,Cellprime_i_matrix)
            Mult=np.dot(Mul,Cellmatrix_inv)
            Multi=np.dot(Mult,Cellprime_j_matrix)
            FM_TE[i,j]+=(2*l+1)*Multi.trace()/2
        if i==j==3:
            FM_TE[i,j]+=1/(0.013**2)
print(FM_TE)
F=np.mat(FM_TE)
FI=F.I
FIarray=np.array(FI)
sigma=np.zeros((6,1))
#get covariance
for i in np.arange(6):
    sigma[i]=np.sqrt(FIarray[i,i])
    print(sigma[i])







