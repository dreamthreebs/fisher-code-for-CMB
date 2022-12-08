import fileinput
import numpy as np
import re
import os
import sys
import time
from matplotlib import pyplot as plt

def get_thetastarmc_from_files():
    filename="output.txt"
    with fileinput.FileInput(filename) as f:
        for line in f :
            if line.startswith('100 theta (CosmoMC)'):
                thetastarmc=re.findall(r"\d+\.?\d*",line[line.find("="):])[0]
                return float(thetastarmc)

def dls2cls(dls):
    cls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls

def set_hubble(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('hubble'):
            print('hubble = '+str(new_value))
        else:
            print(line.strip())

def set_ombh2(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('ombh2'):
            print('ombh2 = '+str(new_value))
        else:
            print(line.strip())

def set_omch2(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('omch2'):
            print('omch2 = '+str(new_value))
        else:
            print(line.strip())

def set_As(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('scalar_amp(1)'):
            print('scalar_amp(1) = '+str(new_value))
        else:
            print(line.strip())

def set_ns(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('scalar_spectral_index(1)'):
            print('scalar_spectral_index(1) = '+str(new_value))
        else:
            print(line.strip())

def set_optical_depth(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('re_optical_depth'):
            print('re_optical_depth = '+str(new_value))
        else:
            print(line.strip())

def set_DM_Pann(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('DM_Pann'):
            print('DM_Pann = '+str(new_value))
        else:
            print(line.strip())

def set_DM_Gamma(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('DM_Gamma'):
            print('DM_Gamma = '+str(new_value))
        else:
            print(line.strip())

def set_DM_mass(new_value):
    for line in fileinput.input('test.ini',inplace=True):
        if line.startswith('DM_mass'):
            print('DM_mass = '+str(new_value))
        else:
            print(line.strip())





#initial value of six fundamental parameters
ombh2_0 = 0.02242
omch2_0 = 0.11933
As_0 = 2.105209331337507e-09 #scalar_amp(1)
ns_0 = 0.9665 #scalar_spectral_index(1)
optical_depth_0 = 0.0561 # re_optical_depth
hubble_0 = 67.66
thetastarmc_0=1.040997
DM_mass_0=1e-4
DM_Pann_0=0
DM_Gamma_0=0
params_num=7

print('thetastarmc_0= ',thetastarmc_0)
#set initial value of fundamental Model
set_ns(ns_0)
set_ombh2(ombh2_0)
set_omch2(omch2_0)
set_optical_depth(optical_depth_0)
set_As(As_0)
set_hubble(hubble_0)
set_DM_mass(DM_mass_0)
set_DM_Pann(DM_Pann_0)
set_DM_Gamma(DM_Gamma_0)

def DM_Pann_prime(params_value):
    set_DM_mass(params_value)
    set_DM_Gamma(0)
    set_DM_Pann(0)
    DM_Pann_0=0
    length=30
    start_footstep=1e-26
    end_footstep=1e-31
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    index=0
    dp=np.zeros((2501,3,length))
    difference=np.ones((2501,3,length))
    DM_Pann_CLprime=np.zeros((2501,3))
    ls=np.arange(2501)
    x1=DM_Pann_0
    set_DM_Pann(x1)
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
    dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    clsn=dls2cls(dlsn)
    
    for h in xcordinate:
        x2=DM_Pann_0+h
        set_DM_Pann(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        dlsp=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        clsp=dls2cls(dlsp)
        dp[:,:,index]=(clsp-clsn)/h
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,2501):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        print('loop=',index)
        index+=1
    '''
        show dp(derivative against params) at difference footstep to find stable point
    '''
    # take footstep at smallest difference as DM_Pann_CLprime
    sum=0
    min_step_mat=np.zeros((2501,3))
    delta_dls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l,i]=minimum_step[0][0]
            DM_Pann_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(2499*3))
    print(best_h)
    set_DM_Pann(xcordinate[best_h])
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])
    dlsd=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    delta_dls=dlsd-dlsn
    DM_Pann_CLprime=dp[:,:,best_h]
    set_DM_Pann(0)
    return DM_Pann_CLprime

def DM_Gamma_prime(params_value):
    set_DM_Pann(0)
    set_DM_mass(params_value)
    set_DM_Gamma(0)
    DM_Gamma_0=0
    length=30
    start_footstep=1e-24
    end_footstep=1e-31
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    index=0
    dp=np.zeros((2501,3,length))
    difference=np.ones((2501,3,length))
    DM_Gamma_CLprime=np.zeros((2501,3))
    ls=np.arange(2501)
    x1=DM_Gamma_0
    set_DM_Gamma(x1)
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
    dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    clsn=dls2cls(dlsn)
    
    for h in xcordinate:
        x2=DM_Gamma_0+h
        set_DM_Gamma(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        dlsp=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        clsp=dls2cls(dlsp)
        dp[:,:,index]=(clsp-clsn)/h
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,2501):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        print('loop=',index)
        index+=1
    '''
        show dp(derivative against params) at different footstep to find stable point
    '''
    # take footstep at smallest difference as DM_Gamma_CLprime
    sum=0
    min_step_mat=np.zeros((2501,3))
    delta_dls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            DM_Gamma_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(2499*3))
    print(best_h)
    set_DM_Gamma(xcordinate[best_h])
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])
    dlsd=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    delta_dls=dlsd-dlsn
    DM_Gamma_CLprime=dp[:,:,best_h]
    set_DM_Gamma(0)
    return DM_Gamma_CLprime


def thetastarmc_prime():
    ls=np.arange(2501)
    length=30
    start_footstep=2e-1
    end_footstep=1e-4
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    thetastar_tmp=np.zeros((length*2,1))
    hubble_tmp=np.zeros((length*2,1))
    dp=np.zeros((2501,3,length))
    difference=np.ones((2501,3,length))
    thetastarmc_CLprime=np.zeros((2501,3))
    index=0
    for h in xcordinate:
        x1=hubble_0-h*hubble_0
        x2=hubble_0+h*hubble_0
        hubble_tmp[index]=x1
        hubble_tmp[2*length-index-1]=x2
        set_hubble(x1)
        os.system('./camb test.ini | grep "100 theta" >output.txt')
        thetastar_tmp[index]=get_thetastarmc_from_files()
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_hubble(x2)
        os.system('./camb test.ini | grep "100 theta" >output.txt')
        thetastar_tmp[2*length-index-1]=get_thetastarmc_from_files()
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(thetastar_tmp[2*length-index-1]-thetastar_tmp[index])
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,2501):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((2501,3))
    delta_dls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            thetastarmc_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(2499*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    thetastarmc_CLprime=dp[:,:,10]
    set_hubble(hubble_0)
    return thetastarmc_CLprime

# thetastarmc_CLprime=thetastarmc_prime()

#derivative against ombh2
def ombh2_prime():
    ls=np.arange(2501)
    length=30
    start_footstep=1e-1
    end_footstep=1e-7
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((2501,3,length))
    difference=np.ones((2501,3,length))
    ombh2_CLprime=np.zeros((2501,3))
    index=0
    for h in xcordinate:
        x1=ombh2_0-h*ombh2_0
        x2=ombh2_0+h*ombh2_0
        set_ombh2(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_ombh2(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*ombh2_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,2501):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((2501,3))
    delta_dls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            ombh2_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(2499*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    ombh2_CLprime=dp[:,:,best_h]
    set_ombh2(ombh2_0)
    return ombh2_CLprime

# ombh2_CLprime=ombh2_prime()

def omch2_prime():
    ls=np.arange(2501)
    length=30
    start_footstep=1e-1
    end_footstep=1e-7
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((2501,3,length))
    difference=np.ones((2501,3,length))
    omch2_CLprime=np.zeros((2501,3))
    index=0
    for h in xcordinate:
        x1=omch2_0-h*omch2_0
        x2=omch2_0+h*omch2_0
        set_omch2(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_omch2(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*omch2_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,2501):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((2501,3))
    delta_dls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            omch2_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(2499*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    omch2_CLprime=dp[:,:,best_h]
    set_omch2(omch2_0)
    return omch2_CLprime

# omch2_CLprime=omch2_prime()

def optical_depth_prime():
    ls=np.arange(2501)
    length=30
    start_footstep=1e-1
    end_footstep=1e-7
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((2501,3,length))
    difference=np.ones((2501,3,length))
    optical_depth_CLprime=np.zeros((2501,3))
    index=0
    for h in xcordinate:
        x1=optical_depth_0-h*optical_depth_0
        x2=optical_depth_0+h*optical_depth_0
        set_optical_depth(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_optical_depth(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*optical_depth_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,2501):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((2501,3))
    delta_dls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            optical_depth_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(2499*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    optical_depth_CLprime=dp[:,:,best_h]
    set_optical_depth(optical_depth_0)
    return optical_depth_CLprime

# optical_depth_CLprime=optical_depth_prime()

def ns_prime():
    ls=np.arange(2501)
    length=30
    start_footstep=1e-1
    end_footstep=1e-7
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((2501,3,length))
    difference=np.ones((2501,3,length))
    ns_CLprime=np.zeros((2501,3))
    index=0
    for h in xcordinate:
        x1=ns_0-h*ns_0
        x2=ns_0+h*ns_0
        set_ns(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_ns(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*ns_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,2501):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((2501,3))
    delta_dls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            ns_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(2499*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    ns_CLprime=dp[:,:,best_h]
    set_ns(ns_0)
    return ns_CLprime

# ns_CLprime=ns_prime()

def As_prime():
    ls=np.arange(2501)
    length=30
    start_footstep=1e-1
    end_footstep=1e-7
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((2501,3,length))
    difference=np.ones((2501,3,length))
    As_CLprime=np.zeros((2501,3))
    index=0
    for h in xcordinate:
        x1=1e-10*np.exp((1-h)*np.log(1e10*As_0))
        x2=1e-10*np.exp((1+h)*np.log(1e10*As_0))
        set_As(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_As(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(np.log(1e10*x2)-np.log(1e10*x1))
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,2501):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((2501,3))
    delta_dls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            As_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(2499*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    As_CLprime=dp[:,:,best_h]
    set_As(As_0)
    return As_CLprime

# As_CLprime=As_prime()

def DM_mass_prime():
    ls=np.arange(2501)
    length=30
    start_footstep=9e-1
    end_footstep=1e-2
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((2501,3,length))
    difference=np.ones((2501,3,length))
    DM_mass_CLprime=np.zeros((2501,3))
    index=0
    for h in xcordinate:
        x1=(1-h)*DM_mass_0
        x2=(1+h)*DM_mass_0
        set_DM_mass(x1)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_n=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_n=dls2cls(DL_n)
        set_DM_mass(x2)
        os.system('./camb test.ini')
        test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
        insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
        DL_p=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
        CL_p=dls2cls(DL_p)
        dp[:,:,index]=(CL_p-CL_n)/(2*h*DM_mass_0)
        for i in np.arange(3):#i is powerspectrum index
            for l in np.arange(2,2501):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        index+=1
        print(index)
    sum=0
    min_step_mat=np.zeros((2501,3))
    delta_dls=np.zeros((2501,3))
    for i in np.arange(3):
        for l in np.arange(2,2501):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l][i]=minimum_step[0][0]
            DM_mass_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            sum+=minimum_step[0][0]
            print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(2499*3))
    print(best_h) # for this parameter, 10 as best_h is better, around 1.45e-2 
    DM_mass_CLprime=dp[:,:,best_h]
    set_DM_mass(DM_mass_0)
    return dp,ls,xcordinate,min_step_mat,DM_mass_CLprime

# dp,ls,xcordinate,min_step_mat,DM_mass_CLprime=DM_mass_prime()

def get_TT_fisher_matrix(Derivative):
    FM_TT=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(30,620):
                Cell[0,0]=clsn[l,0]+Noisecls[l,0]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=Derivative[l,0,i]
                Cellprime_j[0,0]=Derivative[l,0,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM_TT[i,j]+=(2*l+1)*np.trace(Multi)/2
   #             if l==2016:
   #                 print(i,j)
            # if i==j==3:
                # FM_TT[i,j]+=1/(0.013**2)
   # print(FM_TT)
    FM_TT*=fsky
    FI=np.linalg.inv(FM_TT)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        # print(sigma[i])
    return sigma

def get_EE_fisher_matrix(Derivative):
    FM_EE=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(30,620):
                Cell[0,0]=clsn[l,1]+Noisecls[l,1]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=Derivative[l,1,i]
                Cellprime_j[0,0]=Derivative[l,1,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM_EE[i,j]+=(2*l+1)*np.trace(Multi)/2
            # if i==j==3:
                # FM_EE+=1/(0.013**2)
   # print(FM_EE)
    FM_EE*=fsky
    FI=np.linalg.inv(FM_EE)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        # print(sigma[i])
    return sigma

def get_TE_fisher_matrix(Derivative):
    FM_TE=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(30,620):
                Cell[0,0]=np.sqrt(clsn[l,2]**2+(clsn[l,0]+Noisecls[l,0])*(clsn[l,1]+Noisecls[l,1]))
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=Derivative[l,2,i]
                Cellprime_j[0,0]=Derivative[l,2,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM_TE[i,j]+=(2*l+1)*np.trace(Multi)
            # if i==j==3:
                # FM_TE[i,j]+=1/(0.013**2)
   # print(FM_TE)
    FM_TE*=fsky
    FI=np.linalg.inv(FM_TE)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        # print(sigma[i])
    return sigma

def get_combined_fisher_matrix(Derivative):
    FM=np.zeros((params_num,params_num))
    Cell=np.zeros((2,2))
    Cellprime_i=np.zeros((2,2))
    Cellprime_j=np.zeros((2,2))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(30,620):
                Cell[0,0]=clsn[l,0]+Noisecls[l,0]
                Cell[1,0]=clsn[l,2]
                Cell[0,1]=clsn[l,2]
                Cell[1,1]=clsn[l,1]+Noisecls[l,1]
    #            Cell[2,2]=totCL[l,2]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=Derivative[l,0,i]
                Cellprime_i[1,0]=Derivative[l,2,i]
                Cellprime_i[0,1]=Derivative[l,2,i]
                Cellprime_i[1,1]=Derivative[l,1,i]
    #            Cellprime_i[2,2]=Derivative[l,2,i]
                Cellprime_j[0,0]=Derivative[l,0,j]
                Cellprime_j[1,0]=Derivative[l,2,j]
                Cellprime_j[0,1]=Derivative[l,2,j]
                Cellprime_j[1,1]=Derivative[l,1,j]
    #            Cellprime_j[2,2]=Derivative[l,2,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM[i,j]+=(2*l+1)*np.trace(Multi)/2
   # print(FM)
    FM*=fsky
    FI=np.linalg.inv(FM)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        # print(sigma[i])
    return sigma

def cangjs_Noise_level():
    ls=np.arange(2501)
    flag='all'
    Noise=np.zeros((2501,3))
    CMB_temperature=2.725#K
    sensitivity_tem_95=1.45#muK  arcmin
    sensitivity_tem_150=1.45#muK arcmin
    sensitivity_pol_95=2.06#muK  arcmin
    sensitivity_pol_150=2.06#muK arcmin
    thetaFWHM_95=15.37*0.000291#rad
    thetaFWHM_150=9.73*0.000291#rad
    omega_minus_one_95_tem=(sensitivity_tem_95*0.000291)**2
    omega_minus_one_95_pol=(sensitivity_pol_95*0.000291)**2
    omega_minus_one_150_tem=(sensitivity_tem_150*0.000291)**2
    omega_minus_one_150_pol=(sensitivity_pol_150*0.000291)**2
    for i in np.arange(4):
        for l in np.arange(2501):
            if flag=='95Hz':
                if i==0:
                    Noise[l,i]=omega_minus_one_95_tem*np.exp((l*(l+1)*(thetaFWHM_95)**2)/(8*np.log(2)))
                if i==1:
                    Noise[l,i]=omega_minus_one_95_pol*np.exp((l*(l+1)*(thetaFWHM_95)**2)/(8*np.log(2)))
            if flag=='150Hz':
                if i==0:
                    Noise[l,i]=omega_minus_one_150_tem*np.exp((l*(l+1)*(thetaFWHM_150)**2)/(8*np.log(2)))
                if i==1:
                    Noise[l,i]=omega_minus_one_150_pol*np.exp((l*(l+1)*(thetaFWHM_150)**2)/(8*np.log(2)))
            if flag=='all':
                if i==0:
                    Noise[l,i]=1/(1/(omega_minus_one_95_tem*np.exp((l*(l+1)*(thetaFWHM_95)**2)/(8*np.log(2))))+1/(omega_minus_one_150_tem*np.exp((l*(l+1)*(thetaFWHM_150)**2)/(8*np.log(2)))))
                if i==1:
                    Noise[l,i]=1/(1/(omega_minus_one_95_pol*np.exp((l*(l+1)*(thetaFWHM_95)**2)/(8*np.log(2))))+1/(omega_minus_one_150_pol*np.exp((l*(l+1)*(thetaFWHM_150)**2)/(8*np.log(2)))))
    return Noise

def Ali_Noise_level():
    ali_noise_dls=np.loadtxt('Ali_noise.dat',usecols=(1,2))
    print(ali_noise_dls)
    insert_first_thirty_rows=np.zeros((30,2))
    half_Noise=np.insert(ali_noise_dls,0,insert_first_thirty_rows,axis=0)
    insert_last_many_rows=np.zeros((2501-len(half_Noise),2))
    Noised=np.insert(half_Noise,len(half_Noise),insert_last_many_rows,axis=0)
    Noise=noise_dls2cls(Noised)
    return Noise

def noise_dls2cls(dls):
    cls=np.zeros((2501,2))
    for i in np.arange(2):
        for l in np.arange(30,620):
            cls[l,i]=dls[l,i]*(2*np.pi)/(l*(l+1))
    return cls


def initial_totCL():
    os.system('./camb test.ini')
    test_scalCls=np.loadtxt('test_scalCls.dat',usecols=(1,2,3))#only read TT,EE,TE column data
    insert_first_two_rows=np.array([[0,0,0],[0,0,0]])#make dls range from 0 to lmax not 2(derived by camb)
    dlsn=np.insert(test_scalCls,0,insert_first_two_rows,axis=0)
    clsn=dls2cls(dlsn)
    return clsn


"""
    main function
"""

'''
    calculate the derivatives
'''

clsn=initial_totCL()

# DM_Pann_CLprime=DM_Pann_prime(DM_mass_0)
# DM_Gamma_CLprime=DM_Gamma_prime(DM_mass_0)
# ombh2_CLprime=ombh2_prime()
# omch2_CLprime=omch2_prime()
# thetastarmc_CLprime=thetastarmc_prime()
# optical_depth_CLprime=optical_depth_prime()
# ns_CLprime=ns_prime()
# As_CLprime=As_prime()
# np.savez('derivative_data',DM_Pann_CLprime=DM_Pann_CLprime,DM_Gamma_CLprime=DM_Gamma_CLprime,ombh2_CLprime=ombh2_CLprime,omch2_CLprime=omch2_CLprime,thetastarmc_CLprime=thetastarmc_CLprime,optical_depth_CLprime=optical_depth_CLprime,ns_CLprime=ns_CLprime,As_CLprime=As_CLprime)

'''
    load data from npz file
'''

npzdata=np.load('derivative_data.npz')
DM_Pann_CLprime=npzdata['DM_Pann_CLprime']
DM_Gamma_CLprime=npzdata['DM_Gamma_CLprime']
ombh2_CLprime=npzdata['ombh2_CLprime']
omch2_CLprime=npzdata['omch2_CLprime']
thetastarmc_CLprime=npzdata['thetastarmc_CLprime']
optical_depth_CLprime=npzdata['optical_depth_CLprime']
ns_CLprime=npzdata['ns_CLprime']
As_CLprime=npzdata['As_CLprime']



'''
    fisher matrix -> one sigma error
'''

Derivative_Pann=np.zeros((2501,3,params_num))
Derivative_Pann[:,:,0]=ombh2_CLprime
Derivative_Pann[:,:,1]=omch2_CLprime
Derivative_Pann[:,:,2]=thetastarmc_CLprime
Derivative_Pann[:,:,3]=optical_depth_CLprime
Derivative_Pann[:,:,4]=ns_CLprime
Derivative_Pann[:,:,5]=As_CLprime
Derivative_Pann[:,:,6]=DM_Pann_CLprime # choose one of clprime

Derivative_Gamma=np.zeros((2501,3,params_num))
Derivative_Gamma[:,:,0]=ombh2_CLprime
Derivative_Gamma[:,:,1]=omch2_CLprime
Derivative_Gamma[:,:,2]=thetastarmc_CLprime
Derivative_Gamma[:,:,3]=optical_depth_CLprime
Derivative_Gamma[:,:,4]=ns_CLprime
Derivative_Gamma[:,:,5]=As_CLprime
Derivative_Gamma[:,:,6]=DM_Gamma_CLprime


'''
    setup for Noise level and fsky
'''

# Noisecls=np.zeros((2501,3)) # zero noise level
# fsky=1.0

# Noisecls=Ali_Noise_level()
# fsky=0.37

##################################################

# # one sigma error
# sigma_TT=get_TT_fisher_matrix(Derivative_Pann)
# sigma_EE=get_EE_fisher_matrix(Derivative_Pann)
# sigma_TE=get_TE_fisher_matrix(Derivative_Pann)
# sigma_combined=get_combined_fisher_matrix(Derivative_Pann)
# print('Pann_sigma_TT= ',sigma_TT)
# print('Pann_sigma_EE= ',sigma_EE)
# print('Pann_sigma_TE= ',sigma_TE)
# print('Pann_sigma_combined= ',sigma_combined)
# sigma_TT=get_TT_fisher_matrix(Derivative_Gamma)
# sigma_EE=get_EE_fisher_matrix(Derivative_Gamma)
# sigma_TE=get_TE_fisher_matrix(Derivative_Gamma)
# sigma_combined=get_combined_fisher_matrix(Derivative_Gamma)
# print('Gamma_sigma_TT= ',sigma_TT)
# print('Gamma_sigma_EE= ',sigma_EE)
# print('Gamma_sigma_TE= ',sigma_TE)
# print('Gamma_sigma_combined= ',sigma_combined)


'''
    fisher matrix at different noise level
'''

Noisecls_0=Ali_Noise_level()
# Noisecls_0=cangjs_Noise_level()
fsky=0.37
noise_low_fold=0
noise_high_fold=100000
noise_length=10
nsfold=np.linspace(noise_low_fold,noise_high_fold,noise_length)
sigma_of_Gamma=np.zeros((noise_length,4))
sigma_of_Pann=np.zeros((noise_length,4))
print("calculation on different noise level")
for i in range(len(nsfold)):
    print(i)
    Noisecls=Noisecls_0*nsfold[i]
    sigma_of_Gamma[i,0]=get_TT_fisher_matrix(Derivative_Gamma)[params_num-1]
    sigma_of_Gamma[i,1]=get_EE_fisher_matrix(Derivative_Gamma)[params_num-1]
    sigma_of_Gamma[i,2]=get_TE_fisher_matrix(Derivative_Gamma)[params_num-1]
    sigma_of_Gamma[i,3]=get_combined_fisher_matrix(Derivative_Gamma)[params_num-1]
    sigma_of_Pann[i,0]=get_TT_fisher_matrix(Derivative_Pann)[params_num-1]
    sigma_of_Pann[i,1]=get_EE_fisher_matrix(Derivative_Pann)[params_num-1]
    sigma_of_Pann[i,2]=get_TE_fisher_matrix(Derivative_Pann)[params_num-1]
    sigma_of_Pann[i,3]=get_combined_fisher_matrix(Derivative_Pann)[params_num-1]
plt.figure(1)
plt.semilogy(nsfold,sigma_of_Pann[:,0],color='b')
plt.semilogy(nsfold,sigma_of_Pann[:,1],color='g')
plt.semilogy(nsfold,sigma_of_Pann[:,2],color='y')
plt.semilogy(nsfold,sigma_of_Pann[:,3],color='r')
plt.legend(['TT','EE','TE','combined'])
plt.xlabel('noise fold compare to standard noise level')
plt.ylabel('one sigma error for dm_pann at different noise level')
plt.savefig('error_at_dif_noise_Pann.png',dpi=150)
plt.close()
plt.figure(2)
plt.semilogy(nsfold,sigma_of_Gamma[:,0],color='b')
plt.semilogy(nsfold,sigma_of_Gamma[:,1],color='g')
plt.semilogy(nsfold,sigma_of_Gamma[:,2],color='y')
plt.semilogy(nsfold,sigma_of_Gamma[:,3],color='r')
plt.legend(['TT','EE','TE','combined'])
plt.xlabel('noise fold compare to standard noise level')
plt.ylabel('one sigma error for dm_gamma at different noise level')
plt.savefig('error_at_dif_noise_Gamma.png',dpi=150)
plt.close()

'''
    fisher matrix at different fsky
'''

# Noisecls=Ali_Noise_level()
# fsky_start=0.1
# fsky_end=1
# fsky_length=20
# sigma_of_Gamma=np.zeros((fsky_length,4))
# sigma_of_Pann=np.zeros((fsky_length,4))
# fs=np.linspace(fsky_start, fsky_end, fsky_length)
# print("calculation on different fsky:")
# for i in range(len(fs)):
#     fsky=fs[i]
#     print(i)
#     sigma_of_Gamma[i,0]=get_TT_fisher_matrix(Derivative_Gamma)[params_num-1]
#     sigma_of_Gamma[i,1]=get_EE_fisher_matrix(Derivative_Gamma)[params_num-1]
#     sigma_of_Gamma[i,2]=get_TE_fisher_matrix(Derivative_Gamma)[params_num-1]
#     sigma_of_Gamma[i,3]=get_combined_fisher_matrix(Derivative_Gamma)[params_num-1]
#     sigma_of_Pann[i,0]=get_TT_fisher_matrix(Derivative_Pann)[params_num-1]
#     sigma_of_Pann[i,1]=get_EE_fisher_matrix(Derivative_Pann)[params_num-1]
#     sigma_of_Pann[i,2]=get_TE_fisher_matrix(Derivative_Pann)[params_num-1]
#     sigma_of_Pann[i,3]=get_combined_fisher_matrix(Derivative_Pann)[params_num-1]
# plt.figure(1)
# plt.semilogy(fs,sigma_of_Pann[:,0],color='b')
# plt.semilogy(fs,sigma_of_Pann[:,1],color='g')
# plt.semilogy(fs,sigma_of_Pann[:,2],color='y')
# plt.semilogy(fs,sigma_of_Pann[:,3],color='r')
# plt.legend(['TT','EE','TE','combined'])
# plt.xlabel('fsky')
# plt.ylabel('one sigma error for dm_pann')
# plt.savefig('error_at_dif_fsky_pann.png',dpi=150)
# plt.close()
# plt.figure(2)
# plt.semilogy(fs,sigma_of_Gamma[:,0],color='b')
# plt.semilogy(fs,sigma_of_Gamma[:,1],color='g')
# plt.semilogy(fs,sigma_of_Gamma[:,2],color='y')
# plt.semilogy(fs,sigma_of_Gamma[:,3],color='r')
# plt.legend(['TT','EE','TE','combined'])
# plt.xlabel('fsky')
# plt.ylabel('one sigma error for dm_gamma')
# plt.savefig('error_at_dif_fsky_gamma.png',dpi=150)
# plt.close()



# To check if hubble constant have linear relationship with thetastarmc
# plt.plot(hubble_tmp,thetastar_tmp/thetastarmc_0)
# plt.show()

# plt.plot(xcordinate,dp[200,1,:])
# plt.show()
