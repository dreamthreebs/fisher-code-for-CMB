import sys, platform, os
import numpy as np
import time
# import matplotlib
# from matplotlib import pyplot as plt
#Assume installed from github using "git clone --recursive https://github.com/cmbant/CAMB.git"
#This file is then in the docs folders
camb_path = os.path.realpath(os.path.join(os.getcwd(),'..'))
sys.path.insert(0,camb_path)
import camb
from camb import model, initialpower
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))

lmin=2
lmax=2500

#get CMB TT EE BB TE power spectrum
pars=camb.read_ini(os.path.join('/Users/wangyiming/Documents/work/CAMB','inifiles','planck_2018.ini'))
pars.set_accuracy(AccuracyBoost=2.0)
pars.max_l_tensor=1000
# pars.DoLensing=False
As_0=pars.InitPower.As
ns_0=pars.InitPower.ns
r_0=0.02
nt_0=0

pars.WantTensors= True
pars.set_for_lmax(lmax=4000,lens_potential_accuracy=0)
pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=nt_0)
print(pars)
results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
for name in powers: print(name)
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print(totCL.shape)
ls = np.arange(totCL.shape[0])
mat_len = lmax+1

def ten_sca_r_prime(pars,nt_0,r_0):
    r_0=r_0
    length=5
    start_footstep=1e-10
    end_footstep=1e-2
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((mat_len,4,length))
    difference=np.ones((mat_len,4,length))
    ten_sca_r_CLprime=np.zeros((mat_len,4))
    ls=np.arange(mat_len)
    index=0
    for h in xcordinate:
        x1=(1-h)*r_0
        x2=(1+h)*r_0
        pars.InitPower.set_params(As=As_0,ns=ns_0,r=x1,nt=nt_0)
        results_back=camb.get_results(pars)
        powers_back=results_back.get_cmb_power_spectra(pars,lmax=2500,CMB_unit='muK',raw_cl=True)
        cls_back=powers_back['total']
        pars.InitPower.set_params(As=As_0,ns=ns_0,r=x2,nt=nt_0)
        results_forward=camb.get_results(pars)
        powers_forward=results_forward.get_cmb_power_spectra(pars,lmax=2500,CMB_unit='muK',raw_cl=True)
        cls_forward=powers_forward['total']
        dp[:,:,index]=(cls_forward-cls_back)/(2*h*r_0)
        for i in np.arange(4):
            for l in np.arange(2,mat_len):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        print('loop=',index)
        index+=1
    sum=0
    min_step_mat=np.zeros((mat_len,4))
    delta_dls=np.zeros((mat_len,4))
    for i in np.arange(4):
        for l in np.arange(2,mat_len):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l,i]=minimum_step[0][0]
            ten_sca_r_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            if i==2 and 2<=l<=1000: 
                sum+=minimum_step[0][0]
                print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(999))
    print(best_h)
    ten_sca_r_CLprime=dp[:,:,best_h]
    pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=nt_0)
    return dp,ls,xcordinate,min_step_mat,ten_sca_r_CLprime

def nt_prime_forward(pars):
    nt_0=0
    length=10
    start_footstep=1e-10
    end_footstep=1e-2
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((mat_len,4,length))
    difference=np.ones((mat_len,4,length))
    nt_CLprime=np.zeros((mat_len,4))
    ls=np.arange(mat_len)
    index=0
    for h in xcordinate:
        x2=nt_0+h
        pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=x2)
        results_forward=camb.get_results(pars)
        powers_forward=results_forward.get_cmb_power_spectra(pars,lmax=2500,CMB_unit='muK',raw_cl=True)
        cls_forward=powers_forward['total']
        dp[:,:,index]=(cls_forward-totCL)/h
        for i in np.arange(4):
            for l in np.arange(2,mat_len):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        print('loop=',index)
        index+=1
    sum=0
    min_step_mat=np.zeros((mat_len,4))
    delta_dls=np.zeros((mat_len,4))
    for i in np.arange(4):
        for l in np.arange(2,mat_len):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l,i]=minimum_step[0][0]
            nt_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            if i==2 and 2<=l<=1000: 
                sum+=minimum_step[0][0]
                # print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(999))
    print(best_h)
    nt_CLprime=dp[:,:,best_h]
    pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=nt_0)
    return nt_CLprime

def nt_prime(pars,nt_0,r_0):
    nt_0=nt_0
    length=5
    start_footstep=1e-10
    end_footstep=1e-2
    xcordinate=np.geomspace(start_footstep,end_footstep,length)
    dp=np.zeros((mat_len,4,length))
    difference=np.ones((mat_len,4,length))
    nt_CLprime=np.zeros((mat_len,4))
    ls=np.arange(mat_len)
    index=0
    for h in xcordinate:
        x1=nt_0-h
        x2=nt_0+h
        pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=x1)
        results_back=camb.get_results(pars)
        powers_back=results_back.get_cmb_power_spectra(pars,lmax=2500,CMB_unit='muK',raw_cl=True)
        cls_back=powers_back['total']
        pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=x2)
        results_forward=camb.get_results(pars)
        powers_forward=results_forward.get_cmb_power_spectra(pars,lmax=2500,CMB_unit='muK',raw_cl=True)
        cls_forward=powers_forward['total']
        dp[:,:,index]=(cls_forward-cls_back)/(2*h)
        for i in np.arange(4):
            for l in np.arange(2,mat_len):
                if index>=1:
                    if dp[l,i,index-1]!=0.0: 
                        difference[l,i,index]=np.absolute((dp[l,i,index-1]-dp[l,i,index]))/(np.absolute(dp[l,i,index-1]))
                    else:
                        difference[l,i,index]=1.0
        print('loop=',index)
        index+=1
    sum=0
    min_step_mat=np.zeros((mat_len,4))
    delta_dls=np.zeros((mat_len,4))
    for i in np.arange(4):
        for l in np.arange(2,mat_len):
            minimum=np.amin(difference[l,i])
            minimum_step=np.where(difference[l,i]==minimum)
            min_step_mat[l,i]=minimum_step[0][0]
            nt_CLprime[l,i]=dp[l,i,minimum_step[0][0]]
            if i==2 and 2<=l<=1000: 
                sum+=minimum_step[0][0]
                # print('best h is',xcordinate[minimum_step],'at l=',l,'on spectrum',i,'difference=',minimum)
    best_h=round(sum/(999))
    print(best_h)
    nt_CLprime=dp[:,:,best_h]
    pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=nt_0)
    return dp,ls,xcordinate,min_step_mat,nt_CLprime


def get_BB_fisher_matrix(Derivative,params_num):
    FM_BB=np.zeros((params_num,params_num))
    Cell=np.zeros((1,1))
    Cellprime_i=np.zeros((1,1))
    Cellprime_j=np.zeros((1,1))
    for i in np.arange(params_num):
        for j in np.arange(params_num):
            for l in np.arange(50,1000):
                Cell[0,0]=totCL[l,2]+Noisecls[l,2]
                Cell_inv=np.linalg.inv(Cell)
                Cellprime_i[0,0]=Derivative[l,2,i]
                Cellprime_j[0,0]=Derivative[l,2,j]
                Mul=np.matmul(Cell_inv,Cellprime_i)
                Mult=np.matmul(Mul,Cell_inv)
                Multi=np.matmul(Mult,Cellprime_j)
                FM_BB[i,j]+=(2*l+1)*np.trace(Multi)/2
   #             if l==2016:
   #                 print(i,j)
            # if i==j==3:
                # FM_BB[i,j]+=1/(0.013**2)
   # print(FM_BB)
    FM_BB*=fsky
    FI=np.linalg.inv(FM_BB)
    sigma=np.zeros((params_num,1))
    #get covariance
    for i in np.arange(params_num):
        sigma[i]=np.sqrt(FI[i,i])
        print(sigma[i])
    return sigma

"""
    main function
"""
### useless
# nt_CLprime=nt_prime_forward(pars)

'''
    calculate one sigma error at nt=0, r=0.02
'''

# npzdata=np.load('tmp.npz')
# Derivative_nt=npzdata['Derivative_nt']

# time_0=time.time()

# dp,ls,xcordinate,min_step_mat,nt_CLprime=nt_prime(pars,-0.025,0.2)
# dp,ls,xcordinate,min_step_mat,ten_sca_r_CLprime=ten_sca_r_prime(pars,-0.025,0.2)

# fsky=1
# params_num=2
# Noisecls=np.zeros((2501,4))
# Derivative_nt=np.zeros((2501,4,params_num))
# Derivative_nt[:,:,0]=nt_CLprime
# Derivative_nt[:,:,1]=ten_sca_r_CLprime
# pars.InitPower.set_params(As=As_0,ns=ns_0,r=0.2,nt=-0.0025)
# results = camb.get_results(pars)
# powers =results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
# totCL=powers['total']
# sigma=get_BB_fisher_matrix(Derivative_nt,params_num) 
#
# print('calculation time= ', time.time()-time_0)

'''
    calculate one sigma error for nt at different nt fiducial values when r=0.02
'''

# r_0=0.02
# nt_len=5
# sigma_nt=np.zeros((nt_len))
# sigma_r=np.zeros((nt_len))
# nt_set=np.linspace(-3,3,nt_len)
# for i in range(nt_len):
#     nt_0=nt_set[i]
#     dp,ls,xcordinate,min_step_mat,nt_CLprime=nt_prime(pars, nt_0, r_0)
#     # plt.figure(1)
#     # plt.loglog(xcordinate,np.absolute(dp[200,2,:]))
#     # plt.xlabel('nt')
#     # plt.ylabel('dp')
#     # plt.figure(2)
#     # plt.semilogx(ls[2:1000],ls[2:1000]*ls[2:1000]*nt_CLprime[2:1000,2])
#     # plt.xlabel('ls')
#     # plt.ylabel('nt_CLprime')
#     # plt.figure(3)
#     # plt.plot(ls[2:1000],min_step_mat[2:1000,2])
#     # plt.xlabel('ls')
#     # plt.ylabel('min_step_mat')
#     # plt.show()
#     dp,ls,xcordinate,min_step_mat,ten_sca_r_CLprime=ten_sca_r_prime(pars, nt_0, r_0)
#     # plt.figure(1)
#     # plt.loglog(xcordinate,np.absolute(dp[200,2,:]))
#     # plt.xlabel('r')
#     # plt.ylabel('dp')
#     # plt.figure(2)
#     # plt.semilogx(ls[2:1000],ls[2:1000]*ls[2:1000]*nt_CLprime[2:1000,2])
#     # plt.xlabel('ls')
#     # plt.ylabel('ten_sca_r_CLprime')
#     # plt.figure(3)
#     # plt.plot(ls[2:1000],min_step_mat[2:1000,2])
#     # plt.xlabel('ls')
#     # plt.ylabel('min_step_mat')
#     # plt.show()
#     fsky=1
#     params_num=2
#     Noisecls=np.zeros((2501,4))
#     Derivative=np.zeros((2501,4,params_num))
#     Derivative[:,:,0]=nt_CLprime
#     Derivative[:,:,1]=ten_sca_r_CLprime
#     pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=nt_0)
#     results = camb.get_results(pars)
#     powers =results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
#     totCL=powers['total']
#     sigma_nt[i]=get_BB_fisher_matrix(Derivative,params_num)[0] 
#     sigma_r[i]=get_BB_fisher_matrix(Derivative, params_num)[1]
# print(sigma_nt)
# print(sigma_r)
# np.savez('sigma',sigma_nt=sigma_nt,sigma_r=sigma_r)


'''
    calculate sigma at different tensor-scalar ratio when nt=0
'''

# nt_0=0
# r_len=10
# sigma_nt=np.zeros((r_len))
# sigma_r=np.zeros((r_len))
# r_set=np.linspace(1e-3,1e-1,r_len)
# for i in range(r_len):
#     r_0=r_set[i]
#     dp,ls,xcordinate,min_step_mat,nt_CLprime=nt_prime(pars, nt_0, r_0)
#     # plt.figure(1)
#     # plt.loglog(xcordinate,np.absolute(dp[200,2,:]))
#     # plt.xlabel('nt')
#     # plt.ylabel('dp')
#     # plt.figure(2)
#     # plt.semilogx(ls[2:1000],ls[2:1000]*ls[2:1000]*nt_CLprime[2:1000,2])
#     # plt.xlabel('ls')
#     # plt.ylabel('nt_CLprime')
#     # plt.figure(3)
#     # plt.plot(ls[2:1000],min_step_mat[2:1000,2])
#     # plt.xlabel('ls')
#     # plt.ylabel('min_step_mat')
#     # plt.show()
#     dp,ls,xcordinate,min_step_mat,ten_sca_r_CLprime=ten_sca_r_prime(pars, nt_0, r_0)
#     # plt.figure(1)
#     # plt.loglog(xcordinate,np.absolute(dp[200,2,:]))
#     # plt.xlabel('r')
#     # plt.ylabel('dp')
#     # plt.figure(2)
#     # plt.semilogx(ls[2:1000],ls[2:1000]*ls[2:1000]*nt_CLprime[2:1000,2])
#     # plt.xlabel('ls')
#     # plt.ylabel('ten_sca_r_CLprime')
#     # plt.figure(3)
#     # plt.plot(ls[2:1000],min_step_mat[2:1000,2])
#     # plt.xlabel('ls')
#     # plt.ylabel('min_step_mat')
#     # plt.show()
#     fsky=1
#     params_num=2
#     Noisecls=np.zeros((2501,4))
#     Derivative=np.zeros((2501,4,params_num))
#     Derivative[:,:,0]=nt_CLprime
#     Derivative[:,:,1]=ten_sca_r_CLprime
#     pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=nt_0)
#     results = camb.get_results(pars)
#     powers =results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
#     totCL=powers['total']
#     sigma_nt[i]=get_BB_fisher_matrix(Derivative,params_num)[0] 
#     sigma_r[i]=get_BB_fisher_matrix(Derivative, params_num)[1]
# print(sigma_nt)
# print(sigma_r)

# np.savez('sigma_1',sigma_nt=sigma_nt,sigma_r=sigma_r)

# np.savez('derivative_data',nt_CLprime=nt_CLprime,ten_sca_r_CLprime=ten_sca_r_CLprime)

'''
    fisher matrix -> one sigma error
'''

# npzdata=np.load('derivative_data.npz')
# nt_CLprime=npzdata['nt_CLprime']
# ten_sca_r_CLprime=npzdata['ten_sca_r_CLprime']



# fsky=1
# params_num=2
# Noisecls=np.zeros((2501,4))
# Derivative=np.zeros((2501,4,params_num))
# Derivative[:,:,0]=nt_CLprime
# Derivative[:,:,1]=ten_sca_r_CLprime
# sigma=get_BB_fisher_matrix(Derivative,params_num) 


'''
    calculate delta r/r and delta nt/nt at different r, nt and plot it 
'''


# r_len=20
# nt_len=20
# sigma_nt=np.zeros((r_len,nt_len))
# sigma_r=np.zeros((r_len,nt_len))
# r_set=np.linspace(1e-3,5e-2,r_len)
# nt_set=np.linspace(-0.04,0.04,nt_len)
# for j in range(nt_len):
#     for i in range(r_len):
#         r_0=r_set[i]
#         nt_0=nt_set[j]
#         dp,ls,xcordinate,min_step_mat,nt_CLprime=nt_prime(pars, nt_0, r_0)
#         # plt.figure(1)
#         # plt.loglog(xcordinate,np.absolute(dp[200,2,:]))
#         # plt.xlabel('nt')
#         # plt.ylabel('dp')
#         # plt.figure(2)
#         # plt.semilogx(ls[2:1000],ls[2:1000]*ls[2:1000]*nt_CLprime[2:1000,2])
#         # plt.xlabel('ls')
#         # plt.ylabel('nt_CLprime')
#         # plt.figure(3)
#         # plt.plot(ls[2:1000],min_step_mat[2:1000,2])
#         # plt.xlabel('ls')
#         # plt.ylabel('min_step_mat')
#         # plt.show()
#         dp,ls,xcordinate,min_step_mat,ten_sca_r_CLprime=ten_sca_r_prime(pars, nt_0, r_0)
#         # plt.figure(1)
#         # plt.loglog(xcordinate,np.absolute(dp[200,2,:]))
#         # plt.xlabel('r')
#         # plt.ylabel('dp')
#         # plt.figure(2)
#         # plt.semilogx(ls[2:1000],ls[2:1000]*ls[2:1000]*nt_CLprime[2:1000,2])
#         # plt.xlabel('ls')
#         # plt.ylabel('ten_sca_r_CLprime')
#         # plt.figure(3)
#         # plt.plot(ls[2:1000],min_step_mat[2:1000,2])
#         # plt.xlabel('ls')
#         # plt.ylabel('min_step_mat')
#         # plt.show()
#         fsky=1
#         params_num=2
#         Noisecls=np.zeros((2501,4))
#         Derivative=np.zeros((2501,4,params_num))
#         Derivative[:,:,0]=nt_CLprime
#         Derivative[:,:,1]=ten_sca_r_CLprime
#         pars.InitPower.set_params(As=As_0,ns=ns_0,r=r_0,nt=nt_0)
#         results = camb.get_results(pars)
#         powers =results.get_cmb_power_spectra(pars,lmax=2500, CMB_unit='muK',raw_cl=True)
#         totCL=powers['total']
#         sigma_nt[i,j]=get_BB_fisher_matrix(Derivative,params_num)[0] 
#         sigma_r[i,j]=get_BB_fisher_matrix(Derivative, params_num)[1]
#         print('this is loop:i=',i,'j=',j)
# print(sigma_nt)
# print(sigma_r)

# np.savez('sigma_15',sigma_nt=sigma_nt,sigma_r=sigma_r)

'''
    plot data
'''

#def log_tick_formatter(val, pos=None):
#    return r"$10^{{{:.0f}}}$".format(val)

def plot_3d_figure(r_len,nt_len,r_min,r_max,nt_min,nt_max,npzfile):
    from matplotlib import pyplot as plt
    r_len=r_len
    nt_len=nt_len
    r_set=np.linspace(r_min,r_max,r_len)
    nt_set=np.linspace(nt_min,nt_max,nt_len)
    npzdata=np.load(npzfile)
    sigma_nt=npzdata['sigma_nt']
    sigma_r=npzdata['sigma_r']
    X,Y=np.meshgrid(r_set,nt_set,indexing='ij')
    
    #ax = plt.axes(projection='3d')
    #ax.plot_surface(X, Y, np.log10(sigma_nt), rstride=1, cstride=1,
    #                cmap='viridis', edgecolor='none')
    #plt.show()
    
    #plot delta nt / nt
    a=sigma_nt.transpose()
    b=a/nt_set[:,None]
    c=b.transpose()
    ax = plt.axes(projection='3d')
    ax.set_xlabel('r')
    ax.set_ylabel('nt')
    ax.set_zlabel('delta nt/ nt in log10 scale')
    ax.plot_surface(X, Y, np.log10(np.absolute(c)), rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
    # ax.set_zlim(0.2, 1.5)
    plt.show()
    #plot delta r/ r
    ax = plt.axes(projection='3d')
    ax.set_xlabel('r')
    ax.set_ylabel('nt')
    ax.set_zlabel('delta r/ r in log10 scale')
    ax.plot_surface(X, Y, np.log10(sigma_r/r_set[:,None]), rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
    plt.show()

# plot_3d_figure(20,20,1e-3,5e-2,-0.04,0.04,'sigma_12.npz')
# plot_3d_figure(20,20,1e-3,5e-2,-0.04,0.04,'sigma_15.npz')
# plot_3d_figure(40,40,1e-3,5e-2,-0.04,0.04,'sigma_14.npz')

def plot_heat_map(r_len,nt_len,r_min,r_max,nt_min,nt_max,npzfile):
    from matplotlib import pyplot as plt
    r_len=r_len
    nt_len=nt_len
    r_set=np.linspace(r_min,r_max,r_len)
    nt_set=np.linspace(nt_min,nt_max,nt_len)
    npzdata=np.load(npzfile)
    sigma_nt=npzdata['sigma_nt']
    sigma_r=npzdata['sigma_r']
    X,Y=np.meshgrid(r_set,nt_set,indexing='ij')
    
    #ax = plt.axes(projection='3d')
    #ax.plot_surface(X, Y, np.log10(sigma_nt), rstride=1, cstride=1,
    #                cmap='viridis', edgecolor='none')
    #plt.show()
    
    #plot delta nt / nt
    a=sigma_nt.transpose()
    b=a/nt_set[:,None]
    d=b.transpose()
    fig, ax = plt.subplots()
    
    c = ax.pcolormesh(X, Y, np.log10(np.absolute(d)), cmap='viridis',shading='auto')#  'RdBu'
    ax.set_title('delta nt/ nt in log10 scale')
    # set the limits of the plot to the limits of the data
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig.colorbar(c, ax=ax)
    ax.set_xlabel('r')
    ax.set_ylabel('nt')
    # ax.set_zlabel('delta nt / nt in log10 scale')
    plt.show()

    fig, ax = plt.subplots()
    c = ax.pcolormesh(X, Y, np.log10(np.absolute(sigma_r/r_set[:,None])), cmap='viridis',shading='auto')#  'RdBu'
    ax.set_title('delta r/ r in log10 scale')
    # set the limits of the plot to the limits of the data
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig.colorbar(c, ax=ax)
    ax.set_xlabel('r')
    ax.set_ylabel('nt')
    # ax.set_zlabel('delta r/ r in log10 scale')
    plt.show()

    # ax = plt.axes(projection='3d')
    # ax.set_xlabel('r')
    # ax.set_ylabel('nt')
    # ax.set_zlabel('delta r/ r in log10 scale')
    # ax.plot_surface(X, Y, np.log10(sigma_r/r_set[:,None]), rstride=1, cstride=1,
    #                 cmap='viridis', edgecolor='none')
    # plt.show()

# plot_heat_map(20,20,1e-3,2e-2,-0.04,0.04,'sigma_13.npz')
# plot_heat_map(20,20,1e-3,5e-2,-0.04,0.04,'sigma_12.npz')
# plot_heat_map(20,20,1e-3,5e-2,-0.04,0.04,'sigma_15.npz')
# plot_heat_map(40,40,1e-3,5e-2,-0.04,0.04,'sigma_14.npz')














# r_len=5
# nt_len=80
# r_set=np.linspace(1e-3,2e-1,r_len)
# nt_set=np.linspace(-0.025,0.025,nt_len)
# npzdata=np.load('sigma_10.npz')
# sigma_nt=npzdata['sigma_nt']
# sigma_r=npzdata['sigma_r']
#X,Y=np.meshgrid(r_set,nt_set,indexing='ij')

#ax = plt.axes(projection='3d')
#ax.plot_surface(X, Y, np.log10(sigma_nt), rstride=1, cstride=1,
#                cmap='viridis', edgecolor='none')
#plt.show()

##plot delta nt / nt
#X,Y=np.meshgrid(r_set,nt_set,indexing='ij')
#a=sigma_nt.transpose()
#b=a/nt_set[:,None]
#c=b.transpose()
#ax = plt.axes(projection='3d')
#ax.set_xlabel('r')
#ax.set_ylabel('nt')
#ax.set_zlabel('delta nt/ nt in log10 scale')
#ax.plot_surface(X, Y, np.log10(np.absolute(c)), rstride=1, cstride=1,
#                cmap='viridis', edgecolor='none')
#plt.show()
##plot delta r/ r
#ax = plt.axes(projection='3d')
#ax.set_xlabel('r')
#ax.set_ylabel('nt')
#ax.set_zlabel('delta r/ r in log10 scale')
#ax.plot_surface(X, Y, np.log10(sigma_r/r_set[:,None]), rstride=1, cstride=1,
#                cmap='viridis', edgecolor='none')
#plt.show()

'''
    plot totCL and unlensedCL --
'''
# from matplotlib import pyplot as plt
# fig, ax = plt.subplots(2,2, figsize = (12,12))
# ax[0,0].plot(ls,totCL[:,0], color='k')
# ax[0,0].plot(ls,unlensedCL[:,0], color='r')
# ax[0,0].set_title('TT')
# ax[0,1].plot(ls[2:], 1-unlensedCL[2:,0]/totCL[2:,0]);
# ax[0,1].set_title(r'$\Delta TT$')
# ax[1,0].plot(ls,totCL[:,1], color='k')
# ax[1,0].plot(ls,unlensedCL[:,1], color='r')
# ax[1,0].set_title(r'$EE$')
# ax[1,1].plot(ls,totCL[:,3], color='k')
# ax[1,1].plot(ls,unlensedCL[:,3], color='r')
# ax[1,1].set_title(r'$TE$');
# for ax in ax.reshape(-1): ax.set_xlim([2,2500]);


# from matplotlib import pyplot as plt
# fig, ax = plt.subplots(2,2, figsize = (12,12))
# ax[0,0].plot(ls,ls*ls*totCL[:,0], color='k')
# ax[0,0].plot(ls,ls*ls*unlensedCL[:,0], color='r')
# ax[0,0].set_title('TT')
# ax[0,1].plot(ls[2:], 1-unlensedCL[2:,0]/totCL[2:,0]);
# ax[0,1].set_title(r'$\Delta TT$')
# ax[1,0].plot(ls,ls*ls*totCL[:,1], color='k')
# ax[1,0].plot(ls,ls*ls*unlensedCL[:,1], color='r')
# ax[1,0].set_title(r'$EE$')
# ax[1,1].plot(ls,ls*ls*totCL[:,3], color='k')
# ax[1,1].plot(ls,ls*ls*unlensedCL[:,3], color='r')
# ax[1,1].set_title(r'$TE$');
# for ax in ax.reshape(-1): ax.set_xlim([2,2500]);

