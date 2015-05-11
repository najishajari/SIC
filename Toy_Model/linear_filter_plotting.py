# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from linear_filter import linear_filter_SIC,linear_filter_CoM,linear_filter_decreasing_impulse
import numpy as np
from time import time
from matplotlib.ticker import MaxNLocator # added 
dropbox_path='/is/ei/naji/Dropbox/Winter Semester 2014/Master Thesis/'
cluster_path='/agbs/cluster/naji/'
#manuscript_path='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Manuscripts/'
manuscript_path='/Users/Naji/Dropbox/Winter Semster 2014/Master Thesis/Manuscripts/'
tableau20 = np.asarray([(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] )/255.
             
def hist_plot(ksim_size=1000,FO=5,BO=5,series_length=10000,AR_amp_upperbound=.1,noise_amp=.0, win_length=500):
    """The function plots the histogram of the SDR for both directions and their difference for IIR filters with given orders
    input:
        FO:                 refer to linear_filter_SIC()
        BO:                 refer to linear_filter_SIC()
        series_length:      refer to linear_filter_SIC()
        AR_amp_upperbound:  refer to linear_filter_SIC()
        noise_amp:          refer to linear_filter_SIC()
        ksim_size:          number of trials, or more precisely different filters that are tested
    """
    
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 22}
    plt.rc('font', **font)
    #plt.rc('text', usetex=True)
    Dir=np.empty((ksim_size,2,2))
    Dir=linear_filter_SIC(FO=FO,BO=BO,ksim_size=ksim_size,x_size=series_length,AR_amp_upperbound=.10,noise_amp=.0,win_length=win_length)
    fig,axarr=plt.subplots(2,1)
    plt.subplots_adjust(hspace=0.001)
    plt.subplots_adjust(wspace=0.001)
    axarr[0].set_xticklabels([])
    xticklabels = axarr[0].get_xticklabels() 
    plt.setp(xticklabels, visible=False)
    nbins = len(axarr[0].get_xticklabels()) # added 
    axarr[0].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='lower')) # added 
    fig.set_size_inches(10,5.5)
    axarr[0].hist(Dir[:,0,1],alpha=0.5,bins=50,weights=np.ones_like(Dir[:,0,1])/len(Dir[:,0,1]),label=r'$\rho_{\textbf{X}\to \textbf{Y}}$')
    
    axarr[0].hist(Dir[:,1,0],alpha=0.5,bins=50,weights=np.ones_like(Dir[:,1,0])/len(Dir[:,1,0]),label=r'$\rho{\textbf{Y}\to \textbf{X}}$')
    axarr[0].set_xlim(-.6,2.6)
    axarr[0].legend(loc=2,prop={'size':15})

    axarr[1].hist(Dir[:,0,1]-Dir[:,1,0],bins=50,color='r',alpha=.5,weights=np.ones_like(Dir[:,0,1]-Dir[:,1,0])/len(Dir[:,0,1]-Dir[:,1,0]),label=r'$\rho_{\textbf{X}\to \textbf{Y}}-\rho_{\textbf{Y}\to \textbf{X}}$')
    axarr[1].set_xlim(-.6,2.6)
    axarr[1].legend(loc=1,prop={'size':15})
    #plt.tight_layout()
    plt.rc('text', usetex=True)
    fig.subplots_adjust(hspace=0)
    print (np.mean(Dir[:,1,0]<1))
    plt.savefig(manuscript_path+'Linear_Filters/Figures/Delta_Hist-BO=FO=5_nonoise.svg',dpi=500,transparent=True)#, bbox_inches='tight')
    print ("Performance: ",np.mean(Dir[:,0,1]-Dir[:,1,0]>0))
    plt.show()

def CoM_plot(ksim_size=1000):
    Dirs=np.empty((4,ksim_size,2,2))
    for i in range(4):
        print 5*2**i
        Dirs[i]=linear_filter_CoM(FO=1,BO=5*2**i,ksim_size=1000,x_size=10000,AR_amp_upperbound=0.2,noise_amp=0.)
    fig,axarr=plt.subplots(2,2)
    fig.set_size_inches(14.5,10.5)
    axarr[0,0].hist(Dirs[0,:,0,1],bins=70)
    axarr[0,0].set_xlim(-.15,.15)
    axarr[0,0].set_xlabel(r'$FO=5$')
    
    axarr[1,0].hist(Dirs[1,:,0,1],bins=70)
    axarr[1,0].set_xlabel(r'$FO=10$')
    axarr[1,0].set_xlim(-.15,.15)
    
    axarr[0,1].hist(Dirs[2,:,0,1],bins=70)
    axarr[0,1].set_xlim(-.15,.15)
    axarr[0,1].set_xlabel(r'$FO=20$')
    
    axarr[1,1].hist(Dirs[3,:,0,1],bins=70)
    axarr[1,1].set_xlim(-.15,.15)
    axarr[1,1].set_xlabel(r'$FO=40$')
    
    axarr[0,0].set_ylabel(r'$\tilde{\Delta}^\infty_{X\to Y}$')
    axarr[0,0].set_ylim(0,70)
    axarr[0,1].set_ylabel(r'$\tilde{\Delta}^\infty_{X\to Y}$')
    axarr[0,1].set_ylim(0,70)
    axarr[1,0].set_ylabel(r'$\tilde{\Delta}^\infty_{X\to Y}$')
    axarr[1,0].set_ylim(0,70)
    axarr[1,1].set_ylabel(r'$\tilde{\Delta}^\infty_{X\to Y}$')
    axarr[1,1].set_ylim(0,70)

    plt.savefig(manuscript_path+'Linear Filters/Figures/test-CoM-BO=FO=10,nonoise.eps',dpi=500,transparent=True)
    plt.show()
def plot_performance_dim_change(trials=15,xsize=2500,ksim_size=1000,base_order=2,order_increment=2):
    plt.rcParams.update({'font.size': 30})
    plt.rc('text', usetex=True)
    fig=plt.figure()
    fig.set_size_inches(14.5,10.5)
    Dirs1=np.empty((trials,ksim_size,2,2))
    Dirs2=np.empty((trials,ksim_size,2,2))
    Dirs3=np.empty((trials,ksim_size,2,2))
    for i in xrange(trials):
        Dirs1[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=base_order+order_increment*(i+1)-1,BO=base_order*(i+1)-1,AR_amp_upperbound=0.1,noise_amp=0)
        Dirs2[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=0,BO=base_order+order_increment*(i+1)-1,AR_amp_upperbound=0.1,noise_amp=0)
        Dirs3[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=base_order+order_increment*(i+1)-1,BO=1,AR_amp_upperbound=0.1,noise_amp=0)
    suc_prop1=np.mean(Dirs1[:,:,0,1]>Dirs1[:,:,1,0],axis=1)
    suc_prop2=np.mean(Dirs2[:,:,0,1]>Dirs2[:,:,1,0],axis=1)
    suc_prop3=np.mean(Dirs3[:,:,0,1]>Dirs3[:,:,1,0],axis=1)
    
    sample_size=ksim_size
    #calculating the error bars using Wilson score test
    z=1.96
    below_err1=suc_prop1-(1./(1.+(1./(sample_size))*z**2))*(suc_prop1+(1./(2*sample_size))*z**2-z*np.sqrt(1./sample_size*suc_prop1*(1.-suc_prop1)+(1./(4*sample_size**2))*z**2))
    top_err1=(1./(1.+(1./(sample_size))*z**2))*(suc_prop1+(1./(2*sample_size))*z**2+z*np.sqrt(1./sample_size*suc_prop1*(1.-suc_prop1)+(1./(4*sample_size**2))*z**2))-suc_prop1

    below_err2=suc_prop2-(1./(1.+(1./(sample_size))*z**2))*(suc_prop2+(1./(2*sample_size))*z**2-z*np.sqrt(1./sample_size*suc_prop2*(1.-suc_prop2)+(1./(4*sample_size**2))*z**2))
    top_err2=(1./(1.+(1./(sample_size))*z**2))*(suc_prop2+(1./(2*sample_size))*z**2+z*np.sqrt(1./sample_size*suc_prop2*(1.-suc_prop2)+(1./(4*sample_size**2))*z**2))-suc_prop2

    below_err3=suc_prop3-(1./(1.+(1./(sample_size))*z**2))*(suc_prop3+(1./(2*sample_size))*z**2-z*np.sqrt(1./sample_size*suc_prop3*(1.-suc_prop3)+(1./(4*sample_size**2))*z**2))
    top_err3=(1./(1.+(1./(sample_size))*z**2))*(suc_prop3+(1./(2*sample_size))*z**2+z*np.sqrt(1./sample_size*suc_prop3*(1.-suc_prop3)+(1./(4*sample_size**2))*z**2))-suc_prop3

    #plotting the error bars comparing two cases with equal FO and BO vs. FO=0
    # and plotting the error bars comparing two cases with equal FO and BO vs. BO=0
    
    plt.errorbar(np.arange(trials)*order_increment+base_order,suc_prop1,yerr=[below_err1,top_err1],color='b',ecolor='b',fmt='--o',label=r"$FO(\mathcal{S})=BO(\mathcal{S})=FO(\mathcal{S}')=BO(\mathcal{S}')$")
    plt.errorbar(np.arange(trials)*order_increment+base_order,suc_prop2,yerr=[below_err2,top_err2],color='r',ecolor='r',fmt='--o',label=r"$FO(\mathcal{S})=FO(\mathcal{S}')=0$")
    plt.errorbar(np.arange(trials)*order_increment+base_order,suc_prop3,yerr=[below_err3,top_err3],color='g',ecolor='g',fmt='--o',label=r"$BO(\mathcal{S})=BO(\mathcal{S}')=0$")
    output_name='order=[2,20]_x_size=10000_ksim_size=1000_noise=null_fw_dim=0_vs_fw_dim=bw_dim'

    #labeling the axes and legends
    plt.xlabel(r'${\rm Order}$')
    plt.ylabel(r'$\%-{\rm Success}$')    
    plt.legend(loc=4,prop={'size':30})
    #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4, mode="expand", borderaxespad=0., prop={'size':20})    
    plt.xlim(0,22)
    plt.ylim(.5,1.1)
    plt.savefig(manuscript_path+'Linear_Filters/Figures/'+output_name+'.eps',dpi=100,transparent=True)
    
    plt.show()
    return Dirs1,Dirs2,Dirs3
def plot_performance_dim_change_from_file(trials=1000,order_increment=5,param=None):
    """Deprecated!!!"""
    plt.rcParams.update({'font.size': 20})
    plt.rc('text', usetex=True)
    fig=plt.figure()
    #fig.set_size_inches(14.5,10.5)
    #reading data from file
    if param==0:
        f_1=open(dropbox_path+'Linear Filters/Dimension Change/order=[2,20],x_size=10000,ksim_size=1000,noise=null,fw_dim=0')
        f_2=open(dropbox_path+'Linear Filters/Dimension Change/order=[2,20],x_size=10000,ksim_size=1000,noise=null,fw_dim=bw_dim')
    elif param==1:
        f_1=open(dropbox_path+'Linear Filters/Dimension Change/order=[2,20],x_size=10000,ksim_size=1000,noise=null,bw_dim=0')
        f_2=open(dropbox_path+'Linear Filters/Dimension Change/order=[2,20],x_size=10000,ksim_size=1000,noise=null,fw_dim=bw_dim')
    txt_data_1=f_1.read()
    suc_prop_1=np.asarray([float(x) for x in txt_data_1.strip().split(';')])

    #plotting of performance for FO \neq BO
    NAI=1.96*np.sqrt(1./trials*suc_prop_1*(1-suc_prop_1))
    e1=plt.errorbar(np.arange(1,17)*order_increment,suc_prop_1,yerr=NAI,color='b',ecolor='b',fmt='--o')


    txt_data_2=f_2.read()
    suc_prop_2=np.asarray([float(x) for x in txt_data_2.strip().split(';')])

    #plotting of performance for FO=BO
    NAI=1.96*np.sqrt(1./trials*suc_prop_2*(1-suc_prop_2))
    e2=plt.errorbar(np.arange(1,17)*order_increment,suc_prop_2,yerr=NAI,color='r',ecolor='r',fmt='--o')
    
    plt.xlim(order_increment-1,order_increment*(17+1)+1)
    plt.ylabel(r'$\%$-Success')    
    for b in e1[1]:
        b.set_clip_on(False)
    for b in e2[1]:
        b.set_clip_on(False)
    
    if param==0:
        plt.legend([e1,e2],[r'$FO(\mathcal{S})=0$',r'$FO(\mathcal{S})=BO(\mathcal{S})$'],loc=4)
        plt.xlabel(r'$BO(\mathcal{S})$')
        plt.savefig(manuscript_path+'Linear_Filters/Figures/order=[2,20],x_size=10000,ksim_size=1000,noise=null,fw_dim=0 vs fw_dim=bw_dim.eps',dpi=500,transparent=True)
        plt.show()
    elif param==1:
        plt.legend([e1,e2],[r'$BO(\mathcal{S})=0$',r'$FO(\mathcal{S})=BO(\mathcal{S})$'],loc=4)
        plt.xlabel(r'$FO(\mathcal{S})$')
        plt.savefig(manuscript_path+'Linear_Filters/Figures/order=[2,20],x_size=10000,ksim_size=1000,noise=null,bw_dim=0 vs fw_dim=bw_dim.eps',dpi=500,transparent=True)
        plt.show()
def plot_performance_noise_change(trials=500,xsize=5000,ksim_size=100,order=10,noise_step=.05):
    plt.rcParams.update({'font.size': 30})
    print "salam"
    plt.rc('text', usetex=True)
    fig=plt.figure()
    Dirs_1=np.empty((trials,ksim_size,2,2))
    Dirs2=np.empty((trials,ksim_size,2,2))
    Dirs3=np.empty((trials,ksim_size,2,2))
    Dirs4=np.empty((trials,ksim_size,2,2))
    Dirs5=np.empty((trials,ksim_size,2,2))
    Dirs6=np.empty((trials,ksim_size,2,2))
    Dirs_7=np.empty((trials,ksim_size,2,2))
    Dirs_8=np.empty((trials,ksim_size,2,2))
    for i in xrange(trials):
        Dirs_1[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=order,BO=order,AR_amp_upperbound=0.1,noise_amp=noise_step*i,in_noise=False,out_noise=True,denoise=True)
        #Dirs2[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=1,BO=order,AR_amp_upperbound=0.1,noise_amp=noise_step*i,in_noise=False,out_noise=True,denoise=True)
        Dirs3[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=order,BO=order,AR_amp_upperbound=0.1,noise_amp=noise_step*i,in_noise=True,out_noise=True,denoise=True)
        #Dirs4[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=1,BO=order,AR_amp_upperbound=0.1,noise_amp=noise_step*i,in_noise=True,out_noise=True,denoise=True)
        Dirs5[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=order,BO=order,AR_amp_upperbound=0.1,noise_amp=noise_step*i,in_noise=False,out_noise=True,denoise=False)
        #Dirs6[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=1,BO=order,AR_amp_upperbound=0.1,noise_amp=noise_step*i,in_noise=False,out_noise=True,denoise=False)
        Dirs_7[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=order,BO=order,AR_amp_upperbound=0.1,noise_amp=noise_step*i,in_noise=True,out_noise=True,denoise=False)
        #Dirs_8[i]=linear_filter_SIC(ksim_size=ksim_size,x_size=xsize,FO=1,BO=order,AR_amp_upperbound=0.1,noise_amp=noise_step*i,in_noise=True,out_noise=True,denoise=False)
    output_name='order=5_x_size=10000_ksim_size=1000_noise=[0,0_6    ,0_02]_fw_dim=bw_dim'
    suc_prop_1=np.mean(Dirs_1[:,:,0,1]>Dirs_1[:,:,1,0],axis=1)
    
    z=1.96
    below_err1=suc_prop_1-(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop_1+(1./(2*ksim_size))*z**2-z*np.sqrt(1./ksim_size*suc_prop_1*(1.-suc_prop_1)+(1./(4*ksim_size**2))*z**2))
    top_err1=(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop_1+(1./(2*ksim_size))*z**2+z*np.sqrt(1./ksim_size*suc_prop_1*(1.-suc_prop_1)+(1./(4*ksim_size**2))*z**2))-suc_prop_1

    plt.errorbar(np.arange(trials)*noise_step,suc_prop_1,yerr=[below_err1,top_err1],color=tableau20[0],ecolor=tableau20[0],fmt='-o',label=r'${\rm Output\ noise\ with\ denoising}$')

    suc_prop3=np.mean(Dirs3[:,:,0,1]>Dirs3[:,:,1,0],axis=1)
    below_err3=suc_prop3-(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop3+(1./(2*ksim_size))*z**2-z*np.sqrt(1./ksim_size*suc_prop3*(1.-suc_prop3)+(1./(4*ksim_size**2))*z**2))
    top_err3=(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop3+(1./(2*ksim_size))*z**2+z*np.sqrt(1./ksim_size*suc_prop3*(1.-suc_prop3)+(1./(4*ksim_size**2))*z**2))-suc_prop3

    plt.errorbar(np.arange(trials)*noise_step,suc_prop3,yerr=[below_err3,top_err3],color=tableau20[4],ecolor=tableau20[4],fmt='-o',label=r'${\rm Input\ noise\ with\ denoising}$')


    suc_prop5=np.mean(Dirs5[:,:,0,1]>Dirs5[:,:,1,0],axis=1)
    
    below_err5=suc_prop5-(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop5+(1./(2*ksim_size))*z**2-z*np.sqrt(1./ksim_size*suc_prop5*(1.-suc_prop5)+(1./(4*ksim_size**2))*z**2))
    top_err5=(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop5+(1./(2*ksim_size))*z**2+z*np.sqrt(1./ksim_size*suc_prop5*(1.-suc_prop5)+(1./(4*ksim_size**2))*z**2))-suc_prop5

    plt.errorbar(np.arange(trials)*noise_step,suc_prop5,yerr=[below_err5,top_err5],color=tableau20[8],ecolor=tableau20[8],fmt='-o',label=r'${\rm Output\ noise\ w/o\ denoising}$')

#    suc_prop2=np.mean(Dirs2[:,:,0,1]>Dirs2[:,:,1,0],axis=1)
#    below_err2=suc_prop2-(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop2+(1./(2*ksim_size))*z**2-z*np.sqrt(1./ksim_size*suc_prop2*(1.-suc_prop2)+(1./(4*ksim_size**2))*z**2))
#    top_err2=(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop2+(1./(2*ksim_size))*z**2+z*np.sqrt(1./ksim_size*suc_prop2*(1.-suc_prop2)+(1./(4*ksim_size**2))*z**2))-suc_prop2
#
#    plt.errorbar(np.arange(trials)*noise_step,suc_prop2,yerr=[below_err2,top_err2],color=tableau20[2],ecolor=tableau20[2],fmt='--o',label=r'$FO(\mathcal{S})=0$')

    
#    suc_prop4=np.mean(Dirs4[:,:,0,1]>Dirs4[:,:,1,0],axis=1)
#    below_err4=suc_prop4-(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop4+(1./(2*ksim_size))*z**2-z*np.sqrt(1./ksim_size*suc_prop4*(1.-suc_prop4)+(1./(4*ksim_size**2))*z**2))
#    top_err4=(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop4+(1./(2*ksim_size))*z**2+z*np.sqrt(1./ksim_size*suc_prop4*(1.-suc_prop4)+(1./(4*ksim_size**2))*z**2))-suc_prop4
#
#    plt.errorbar(np.arange(trials)*noise_step,suc_prop4,yerr=[below_err4,top_err4],color=tableau20[6],ecolor=tableau20[6],fmt='--o',label=r'$FO(\mathcal{S})=0,{\rm\ with\ input\ noise}$')
#    
#
#    suc_prop6=np.mean(Dirs6[:,:,0,1]>Dirs6[:,:,1,0],axis=1)
#    below_err6=suc_prop6-(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop6+(1./(2*ksim_size))*z**2-z*np.sqrt(1./ksim_size*suc_prop6*(1.-suc_prop6)+(1./(4*ksim_size**2))*z**2))
#    top_err6=(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop6+(1./(2*ksim_size))*z**2+z*np.sqrt(1./ksim_size*suc_prop6*(1.-suc_prop6)+(1./(4*ksim_size**2))*z**2))-suc_prop6
#
#    plt.errorbar(np.arange(trials)*noise_step,suc_prop6,yerr=[below_err6,top_err6],color=tableau20[10],ecolor=tableau20[10],fmt='--o',label=r'$FO(\mathcal{S})=0, {\rm w/o\ denoising}$')

    suc_prop_7=np.mean(Dirs_7[:,:,0,1]>Dirs_7[:,:,1,0],axis=1)
    below_err_7=suc_prop_7-(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop_7+(1./(2*ksim_size))*z**2-z*np.sqrt(1./ksim_size*suc_prop_7*(1.-suc_prop_7)+(1./(4*ksim_size**2))*z**2))
    top_err_7=(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop_7+(1./(2*ksim_size))*z**2+z*np.sqrt(1./ksim_size*suc_prop_7*(1.-suc_prop_7)+(1./(4*ksim_size**2))*z**2))-suc_prop_7

    plt.errorbar(np.arange(trials)*noise_step,suc_prop_7,yerr=[below_err_7,top_err_7],color=tableau20[2],ecolor=tableau20[2],fmt='-o',label=r'${\rm Input\ \&\ output\ noise\ w/o\ denoising}$')

#    suc_prop_8=np.mean(Dirs_8[:,:,0,1]>Dirs_8[:,:,1,0],axis=1)
#    below_err_8=suc_prop_8-(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop_8+(1./(2*ksim_size))*z**2-z*np.sqrt(1./ksim_size*suc_prop_8*(1.-suc_prop_8)+(1./(4*ksim_size**2))*z**2))
#    top_err_8=(1./(1.+(1./(2*ksim_size))*z**2))*(suc_prop_8+(1./(2*ksim_size))*z**2+z*np.sqrt(1./ksim_size*suc_prop_8*(1.-suc_prop_8)+(1./(4*ksim_size**2))*z**2))-suc_prop_8
#
#    plt.errorbar(np.arange(trials)*noise_step,suc_prop_8,yerr=[below_err_8,top_err_8],color=tableau20[10],ecolor=tableau20[10],fmt='--o',label=r'$FO(\mathcal{S})=0, {\rm w/o\ denoising}$')
    
            
    ###reporting an average of variance
    #plt.plot(trials*[np.mean(Dirs_1[:,:,0,0])],np.arange(trials)*noise_step,label=r'${\rm input\ variance}$')
    #plt.plot(trials*[np.mean(Dirs_1[:,:,1,1])],np.arange(trials)*noise_step,label=r'${\rm output\ variance}$')

    #f=open('/is/ei/naji/Dropbox/Winter Semester 2014/Master Thesis/Manuscripts/Linear_Filters/Figures'+output_name,'w')
    #print (np.mean(Dirs[:,:,0,1]>Dirs[:,:,1,0],axis=1),file=f)
    #f.close()
    plt.xlim(-0.02)
    plt.xlabel(r'$\sigma$')
    plt.ylabel(r'$\%$-Success')    
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0., prop={'size':20})    
    plt.gcf().set_size_inches(15,11)
    #plt.tight_latyout()
    plt.rcParams['xtick.major.pad']='20'
    plt.savefig(dropbox_path+'Manuscripts/Linear_Filters/Figures/'+output_name+'.eps',dpi=500,transparent=True,bbox_inches='tight')
    plt.show()
#plot_performance_dim_change(param=1)
#CoM_plot()
t=time()
#plot_performance_noise_change()
hist_plot(ksim_size=1000,BO=1,FO=30,series_length=100,win_length=50)
#plot_performance_noise_change()
#hist_plot()
#plot_performance_dim_change()
#plot_performance_dim_change(trials=10,ksim_size=1000)
#plot_performance_noise_change(trials=25,order=5,noise_step=0.1)
#plot_performance_dim_change(trials=10,order_increment=2)
#plot_performance_dim_change()
#plot_performance_dim_change()
#plot_performance_dim_change(trials=10)
#hist_plot()
print time()-t