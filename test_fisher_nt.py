#Apply log parameters variation to get stable derivative
import sys, platform, os
#import matplotlib
#from matplotlib import pyplot as plt
import numpy as np
# import threading
# from multiprocessing import Process, cpu_count, Queue, Pool
import time
#from scipy.interpolate import CubicSpline
#import scipy 
camb_path = os.path.realpath(os.path.join(os.getcwd(),'..'))
sys.path.insert(0,camb_path)
import camb
from camb import model, initialpower
#from function import *
#from decimal import Decimal,getcontext 
#getcontext().prec = 15
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))
pars = camb.CAMBparams()

def set_parameters(pars):
    pars.set_cosmology(ombh2=0.022032,omch2=0.12038,thetastar=0.0104119,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
    pars.InitPower.set_params(ns=0.9619,As=2.21536e-9,nrun=0)
    pars.set_for_lmax(lmax=2500,max_eta_k=10000,lens_potential_accuracy=0)
    pars.set_nonlinear_lensing(False)

def get_initial_value(pars,results):
    #get thetastar(one of H0, thetastar, cosmotheta_mc must be given in set_cosmology function)
    derived_parameters=results.get_derived_params()
    thetastar_0=derived_parameters['thetastar']/100
    #get initial value of 6 standard parameters(except thetastar)
    ombh2_0=pars.ombh2
    omch2_0=pars.omch2
    optical_depth_0=pars.Reion.optical_depth
    ns_0=pars.InitPower.ns
    As_0=pars.InitPower.As
    return ombh2_0,omch2_0,thetastar_0,optical_depth_0,ns_0,As_0

def calculate_derivative_and_difference(x1,x2,x3,x4,y1_totCL,y2_totCL,y3_totCL,y4_totCL,dp,difference,initial_value,index):
    from scipy.interpolate import CubicSpline
    for i in np.arange(4):#i is powerspectrum index
        for l in np.arange(30,2501):
            x=np.array([x1,x2,initial_value,x3,x4])
            y=np.array([y1_totCL[l,i],y2_totCL[l,i],totCL[l,i],y3_totCL[l,i],y4_totCL[l,i]])
            cs=CubicSpline(x,y,bc_type='natural')
            dp[l,i,index]=cs.derivative()(initial_value)
            #print('dp is',dp[index])
            if index>=1:
                difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
    return dp,difference
                #print('difference=',difference[index])
               # if difference[l,i,index-1]<1e-5:#here we treat '2' as a mark,telling us the difference following (<1e-4) should not be calculated
               #     difference[l,i,index]=2
               # elif difference[l,i,index-1]==2:
               #     difference[l,i,index]=2

def get_derivative_spectrum(dp,difference,derivative_spectrum):
    for i in np.arange(4):
        for l in np.arange(30,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            derivative_spectrum[l,i]=dp[l,i,minimum_step[0][0]]
            #print('best h is',xcordinate[minimum_step[0][0]],'at l=',l,'on spectrum',i,'difference=',minimum)
    return derivative_spectrum

def get_TT_fisher_matrix(Derivative):
    FM_TT=np.zeros((6,6))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(6):
        for j in np.arange(6):
            for l in np.arange(30,2501):
                Cell[0,0]=totCL[l,0]+Noise[l,0]
                Cellm=np.mat(Cell)
                Cellmatrix_inv=Cellm.I
                Cellprime_i[0,0]=Derivative[l,0,i]
                Cellprime_i_matrix=np.mat(Cellprime_i)
                Cellprime_j[0,0]=Derivative[l,0,j]
                Cellprime_j_matrix=np.mat(Cellprime_j)
                Mul=np.dot(Cellmatrix_inv,Cellprime_i_matrix)
                Mult=np.dot(Mul,Cellmatrix_inv)
                Multi=np.dot(Mult,Cellprime_j_matrix)
                FM_TT[i,j]+=(2*l+1)*Multi.trace()/2
   #             if l==2016:
   #                 print(i,j)
            if i==j==3:
                FM_TT[i,j]+=1/(0.013**2)
   # print(FM_TT)
    FM_TT*=fsky
    F=np.mat(FM_TT)
    FI=F.I
    FIarray=np.array(FI)
    sigma=np.zeros((6,1))
    #get covariance
    for i in np.arange(6):
        sigma[i]=np.sqrt(FIarray[i,i])
        print(sigma[i])
    #get conditional errors
   # conditionalerror=np.zeros((6,1))
   # for i in np.arange(6):
   #     conditionalerror[i]=1/(np.sqrt(FM_TT[i,i]))
   #     print(conditionalerror[i])

def get_EE_fisher_matrix(Derivative):
    FM_EE=np.zeros((6,6))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(6):
        for j in np.arange(6):
            for l in np.arange(30,2500):
                Cell[0,0]=totCL[l,1]+Noise[l,1]
                Cellm=np.mat(Cell)
                Cellmatrix_inv=Cellm.I
                Cellprime_i[0,0]=Derivative[l,1,i]
                Cellprime_i_matrix=np.mat(Cellprime_i)
                Cellprime_j[0,0]=Derivative[l,1,j]
                Cellprime_j_matrix=np.mat(Cellprime_j)
                Mul=np.dot(Cellmatrix_inv,Cellprime_i_matrix)
                Mult=np.dot(Mul,Cellmatrix_inv)
                Multi=np.dot(Mult,Cellprime_j_matrix)
                FM_EE[i,j]+=(2*l+1)*Multi.trace()/2
            if i==j==3:
                FM_EE+=1/(0.013**2)
   # print(FM_EE)
    FM_EE*=fsky
    F=np.mat(FM_EE)
    FI=F.I
    FIarray=np.array(FI)
    sigma=np.zeros((6,1))
    #get covariance
    for i in np.arange(6):
        sigma[i]=np.sqrt(FIarray[i,i])
        print(sigma[i])

def get_TE_fisher_matrix(Derivative):
    FM_TE=np.zeros((6,6))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(6):
        for j in np.arange(6):
            for l in np.arange(30,2501):
                Cell[0,0]=np.sqrt(totCL[l,3]**2+(totCL[l,0]+Noise[l,0])*(totCL[l,1]+Noise[l,0]))
                Cellm=np.mat(Cell)
                Cellmatrix_inv=Cellm.I
                Cellprime_i[0,0]=Derivative[l,3,i]
                Cellprime_i_matrix=np.mat(Cellprime_i)
                Cellprime_j[0,0]=Derivative[l,3,j]
                Cellprime_j_matrix=np.mat(Cellprime_j)
                Mul=np.dot(Cellmatrix_inv,Cellprime_i_matrix)
                Mult=np.dot(Mul,Cellmatrix_inv)
                Multi=np.dot(Mult,Cellprime_j_matrix)
                FM_TE[i,j]+=(2*l+1)*Multi.trace()
            if i==j==3:
                FM_TE[i,j]+=1/(0.013**2)
   # print(FM_TE)
    FM_TE*=fsky 
    F=np.mat(FM_TE)
    FI=F.I
    FIarray=np.array(FI)
    sigma=np.zeros((6,1))
    #get covariance
    for i in np.arange(6):
        sigma[i]=np.sqrt(FIarray[i,i])
        print(sigma[i])

def get_combined_fisher_matrix(Derivative):
    FM=np.zeros((6,6))
    Cell=np.zeros((2,2))
    Cellprime_i=np.zeros((2,2))
    Cellprime_j=np.zeros((2,2))
    for i in np.arange(6):
        for j in np.arange(6):
            for l in np.arange(30,2501):
                Cell[0,0]=totCL[l,0]+Noise[l,0]
                Cell[1,0]=totCL[l,3]
                Cell[0,1]=totCL[l,3]
                Cell[1,1]=totCL[l,1]+Noise[l,1]
    #            Cell[2,2]=totCL[l,2]
                Cellm=np.mat(Cell)
                Cellmatrix_inv=Cellm.I
                Cellprime_i[0,0]=Derivative[l,0,i]
                Cellprime_i[1,0]=Derivative[l,3,i]
                Cellprime_i[0,1]=Derivative[l,3,i]
                Cellprime_i[1,1]=Derivative[l,1,i]
    #            Cellprime_i[2,2]=Derivative[l,2,i]
                Cellprime_i_matrix=np.mat(Cellprime_i)
                Cellprime_j[0,0]=Derivative[l,0,j]
                Cellprime_j[1,0]=Derivative[l,3,j]
                Cellprime_j[0,1]=Derivative[l,3,j]
                Cellprime_j[1,1]=Derivative[l,1,j]
    #            Cellprime_j[2,2]=Derivative[l,2,j]
                Cellprime_j_matrix=np.mat(Cellprime_j)
                Mul=np.dot(Cellmatrix_inv,Cellprime_i_matrix)
                Mult=np.dot(Mul,Cellmatrix_inv)
                Multi=np.dot(Mult,Cellprime_j_matrix)
                FM[i,j]+=(2*l+1)*Multi.trace()/2
   # print(FM)
    FM*=fsky
    F=np.mat(FM)
    FI=F.I
    FIarray=np.array(FI)
    sigma=np.zeros((6,1))
    #get covariance
    for i in np.arange(6):
        sigma[i]=np.sqrt(FIarray[i,i])
        print(sigma[i])

#set parameters
set_parameters(pars)

#get results TT EE BB TE from number 0 to 3 
results = camb.get_results(pars)
powers=results.get_cmb_power_spectra(pars,lmax=2500,CMB_unit='muK',raw_cl=True)
for name in powers: print(name)
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print(totCL.shape)
print(totCL[100,0])
ls=np.arange(totCL.shape[0])

#initial value of parameters
ombh2_0,omch2_0,thetastar_0,optical_depth_0,ns_0,As_0=get_initial_value(pars,results)

#setup for interpolation
length=4
start_footstep=0.05
end_footstep=1e-8
xcordinate=np.geomspace(start_footstep,end_footstep,length)
#ombh2_CLprime=np.zeros((2501,4))
#omch2_CLprime=np.zeros((2501,4))
#thetastar_CLprime=np.zeros((2501,4))
#optical_depth_CLprime=np.zeros((2501,4))
#ns_CLprime=np.zeros((2501,4))
#As_CLprime=np.zeros((2501,4))
time_0=time.perf_counter()
'''
    if we want to apply multiprocessing on six or more derivative calculation,there are several steps to be refered to:
    1. Make your calculation a function
    2. Transfer several child processes to run these function  
    3. If you have some data which could be used in the following program, use queue in multiprocessing package. q.put() and q.get()
    4. I do not know why if I don't terminate the child processes, the main process will not continue. So I terminate all child processes
    5. Because different child processes run their processes at the same time, the data you want to return will not be chronologically, so here I return data as a dictionary and use dict=(**x,**y) to combine these data together
    6. Pool in multiprocessing package may be more useful and reasonable. However, errors merged if Pool.get() method are used. That might come from some program mistakes in my codes
'''
#using CubicSpline method to get derivative
#Derivative of partial Cl partial ombh2
def ombh2_prime(q,pars):
    index=0
    dp=np.zeros((2501,4,length))
    difference=np.ones((2501,4,length))
    ombh2_CLprime=np.zeros((2501,4))
    for h in xcordinate:
        x1=ombh2_0-2*h*ombh2_0
        x2=ombh2_0-h*ombh2_0
        x3=ombh2_0+h*ombh2_0
        x4=ombh2_0+2*h*ombh2_0
        pars.set_cosmology(thetastar=thetastar_0,ombh2=x1,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
        y1_results=camb.get_results(pars)
        y1_powers=y1_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y1_totCL=y1_powers['total']
        pars.set_cosmology(thetastar=thetastar_0,ombh2=x2,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
        y2_results=camb.get_results(pars)
        y2_powers=y2_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y2_totCL=y2_powers['total']
        pars.set_cosmology(thetastar=thetastar_0,ombh2=x3,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
        y3_results=camb.get_results(pars)
        y3_powers=y3_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y3_totCL=y3_powers['total']
        pars.set_cosmology(thetastar=thetastar_0,ombh2=x4,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
        y4_results=camb.get_results(pars)
        y4_powers=y4_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y4_totCL=y4_powers['total']
        dp,difference=calculate_derivative_and_difference(x1,x2,x3,x4,y1_totCL,y2_totCL,y3_totCL,y4_totCL,dp,difference,ombh2_0,index)
        print('loop=',index)
        index+=1
    ombh2_CLprime=get_derivative_spectrum(dp,difference,ombh2_CLprime)
    pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
    print('ombh2 is ok')
    #return True
    d={}
    d['ombh2_CLprime']=ombh2_CLprime
    q.put(d,block=False)
#ombh2_CLprime=ombh2_prime(pars)
##Derivative of partial Cl partial omch2
def omch2_prime(q,pars):
    index=0
    dp=np.zeros((2501,4,length))
    difference=np.ones((2501,4,length))
    omch2_CLprime=np.zeros((2501,4))
    for h in xcordinate:
        x1=omch2_0-2*h*omch2_0
        x2=omch2_0-h*omch2_0
        x3=omch2_0+h*omch2_0
        x4=omch2_0+2*h*omch2_0
        pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=x1,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
        y1_results=camb.get_results(pars)
        y1_powers=y1_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y1_totCL=y1_powers['total']
        pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=x2,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
        y2_results=camb.get_results(pars)
        y2_powers=y2_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y2_totCL=y2_powers['total']
        pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=x3,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
        y3_results=camb.get_results(pars)
        y3_powers=y3_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y3_totCL=y3_powers['total']
        pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=x4,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
        y4_results=camb.get_results(pars)
        y4_powers=y4_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y4_totCL=y4_powers['total']
        dp,difference=calculate_derivative_and_difference(x1,x2,x3,x4,y1_totCL,y2_totCL,y3_totCL,y4_totCL,dp,difference,omch2_0,index)
        print('loop=',index)
        index+=1
    omch2_CLprime=get_derivative_spectrum(dp,difference,omch2_CLprime)
    pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
    print('omch2 is ok')
    #return omch2_CLprime
    d={}
    d['omch2_CLprime']=omch2_CLprime
    q.put(d,block=False)
#omch2_CLprime=omch2_prime(pars)

#Derivative of partial Cl partial thetastar
def thetastar_prime(q,pars):
    index=0
    dp=np.zeros((2501,4,length))
    difference=np.ones((2501,4,length))
    thetastar_CLprime=np.zeros((2501,4))
    for h in xcordinate:
        x1=(1-2*h)*thetastar_0*100
        x2=(1-h)*thetastar_0*100
        x3=(1+h)*thetastar_0*100
        x4=(1+2*h)*thetastar_0*100
        pars.set_cosmology(thetastar=x1/100,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0,theta_H0_range=(1,200))
        y1_results=camb.get_results(pars)
        y1_powers=y1_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y1_totCL=y1_powers['total']
        pars.set_cosmology(thetastar=x2/100,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0,theta_H0_range=(1,200))
        y2_results=camb.get_results(pars)
        y2_powers=y2_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y2_totCL=y2_powers['total']
        pars.set_cosmology(thetastar=x3/100,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0,theta_H0_range=(1,200))
        y3_results=camb.get_results(pars)
        y3_powers=y3_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y3_totCL=y3_powers['total']
        pars.set_cosmology(thetastar=x4/100,omch2=omch2_0,ombh2=ombh2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0,theta_H0_range=(1,200))
        y4_results=camb.get_results(pars)
        y4_powers=y4_results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
        y4_totCL=y4_powers['total']
        dp,difference=calculate_derivative_and_difference(x1,x2,x3,x4,y1_totCL,y2_totCL,y3_totCL,y4_totCL,dp,difference,100*thetastar_0,index)
        print('loop=',index)
        index+=1
    thetastar_CLprime=get_derivative_spectrum(dp,difference,thetastar_CLprime)
    pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
    print('thetastar is ok')
    #return thetastar_CLprime
    d={}
    d['thetastar_CLprime']=thetastar_CLprime
    q.put(d,block=False)
#thetastar_CLprime=thetastar_prime(pars)


#Derivative of partial Cl partial tau
def optical_depth_prime(q,pars):
    index=0
    dp=np.zeros((2501,4,length))
    difference=np.ones((2501,4,length))
    optical_depth_CLprime=np.zeros((2501,4))
    for h in xcordinate:
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
        dp,difference=calculate_derivative_and_difference(x1,x2,x3,x4,y1_totCL,y2_totCL,y3_totCL,y4_totCL,dp,difference,optical_depth_0,index)
        print('loop=',index)
        index+=1
    optical_depth_CLprime=get_derivative_spectrum(dp,difference,optical_depth_CLprime)
    pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
    print('tau is ok')
    #return optical_depth_CLprime
    d={}
    d['optical_depth_CLprime']=optical_depth_CLprime
    q.put(d,block=False)
#optical_depth_CLprime=optical_depth_prime(pars)


#Derivative of partial Cl partial ns
def ns_prime(q,pars):
    index=0
    dp=np.zeros((2501,4,length))
    difference=np.ones((2501,4,length))
    ns_CLprime=np.zeros((2501,4))
    for h in xcordinate:
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
        dp,difference=calculate_derivative_and_difference(x1,x2,x3,x4,y1_totCL,y2_totCL,y3_totCL,y4_totCL,dp,difference,ns_0,index)
        print('loop=',index)
        index+=1
    ns_CLprime=get_derivative_spectrum(dp,difference,ns_CLprime)
    pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
    print('ns is ok')
    #return ns_CLprime
    d={}
    d['ns_CLprime']=ns_CLprime
    q.put(d,block=False)
#ns_CLprime=ns_prime(pars)

#Derivative of partial Cl partial As
def As_prime(q,pars):
    index=0
    dp=np.zeros((2501,4,length))
    difference=np.ones((2501,4,length))
    As_CLprime=np.zeros((2501,4))
    for h in xcordinate:
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
        dp,difference=calculate_derivative_and_difference(x1,x2,x3,x4,y1_totCL,y2_totCL,y3_totCL,y4_totCL,dp,difference,np.log(1e10*As_0),index)
        print('loop=',index)
        index+=1
    As_CLprime=get_derivative_spectrum(dp,difference,As_CLprime)
    pars.set_cosmology(thetastar=thetastar_0,ombh2=ombh2_0,omch2=omch2_0,tau=0.0925,omk=0.0,YHe=0.24,mnu=0)
    print('As is ok')
    #return list(As_CLprime)
    # d={}
    # d['As_CLprime']=As_CLprime
    # q.put(d,block=False)

##get derivative
##DerivativeTT={0:ombh2_CLprime[:,0],1:omch2_CLprime[:,0],2:thetastar_CLprime[:,0],3:optical_depth_CLprime[:,0],4:ns_CLprime[:,0],5:As_CLprime[:,0]}
##DerivativeEE={0:ombh2_CLprime[:,1],1:omch2_CLprime[:,1],2:thetastar_CLprime[:,1],3:optical_depth_CLprime[:,1],4:ns_CLprime[:,1],5:As_CLprime[:,1]}
##DerivativeBB={0:ombh2_CLprime[:,2],1:omch2_CLprime[:,2],2:thetastar_CLprime[:,2],3:optical_depth_CLprime[:,2],4:ns_CLprime[:,2],5:As_CLprime[:,2]}
##DerivativeTE={0:ombh2_CLprime[:,3],1:omch2_CLprime[:,3],2:thetastar_CLprime[:,3],3:optical_depth_CLprime[:,3],4:ns_CLprime[:,3],5:As_CLprime[:,3]}
#Derivative=np.zeros((2501,4,6))
#Derivative[:,:,0]=ombh2_CLprime
#Derivative[:,:,1]=omch2_CLprime
#Derivative[:,:,2]=thetastar_CLprime
#Derivative[:,:,3]=optical_depth_CLprime
#Derivative[:,:,4]=ns_CLprime
#Derivative[:,:,5]=As_CLprime

##TT fisher matrix
#get_TT_fisher_matrix(Derivative)
##EE fisher matrix
#get_EE_fisher_matrix(Derivative)
##TE fisher matrix
#get_TE_fisher_matrix(Derivative)
##combined fisher matrix
#get_combined_fisher_matrix(Derivative)
##combined fisher matrix
##get_combined_fisher_matrix(Derivative)




