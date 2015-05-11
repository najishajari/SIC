import numpy as np
from matplotlib import pyplot as plt
res_SIC=np.empty((32,32,300))
res_GC=np.empty((32,32,300))

cluster_path='/agbs/cluster/naji/Linear Filters/real_data'
from matplotlib.colors import BoundaryNorm

for i in range(32):
    for j in range(32):
        f=open(cluster_path+'/vvp01/2006-4-9_17-29-30/output/'+str(i)+','+str(j)+'.txt')
        f_GC=open(cluster_path+'/vvp01/2006-4-9_17-29-30/GC_output/'+str(i)+','+str(j)+'.txt')
        for k in range(300):
            line= f.readline()
            line_GC=f_GC.readline()
            line=line.split(';')
            line_GC=line_GC.split(';')
            
            #res[i,j]= np.abs(np.float64(line)[4])<np.abs(np.float64(line)[5])
            res_SIC[i,j,k]= np.float64(line)[0]>np.float64(line)[1]
            #res_2[i,j,k]= np.abs(np.float64(line)[2])<np.abs(np.float64(line)[3])
            res_GC[i,j,k]= np.float64(line_GC)[0]<np.float64(line_GC)[1]
            
manuscript_path='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Manuscripts/'

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}
fig=plt.figure()
plt.rc('font', **font)
plt.rc('text',usetex=True)
#fig.set_size_inches(12.5,16.5)
res_mean=np.squeeze(np.apply_over_axes(np.mean,res_SIC,[0,1]))
res_GC_mean=np.squeeze(np.apply_over_axes(np.mean,res_GC,[0,1]))
#plt.plot(np.squeeze(np.apply_over_axes(np.mean,res_2,[0,1])))
#plt.plot(np.squeeze(np.apply_over_axes(np.mean,res_3,[0,1])))
z=1.96
sample_size=1024
trials=res_mean.shape[0]
print (1./(1.+(1./(2*sample_size))*z**2))
below_err=res_mean-(1./(1.+(1./(2*sample_size))*z**2))*(res_mean+(1./(2*sample_size))*z**2-z*np.sqrt(1./sample_size*res_mean*(1.-res_mean)+(1./(4*sample_size**2))*z**2))
#print below_err
top_err=(1./(1.+(1./(2*sample_size))*z**2))*(res_mean+(1./(2*sample_size))*z**2+z*np.sqrt(1./sample_size*res_mean*(1.-res_mean)+(1./(4*sample_size**2))*z**2))-res_mean
e=plt.errorbar(np.arange(2,(trials+1)*2,2),res_mean,yerr=[below_err,top_err],color='b',ecolor='#990000',fmt='-o',label=r'${\rm SIC}$')
plt.plot([0.5]*trials*2,'--')
below_err_GC=res_GC_mean-(1./(1.+(1./(2*sample_size))*z**2))*(res_GC_mean+(1./(2*sample_size))*z**2-z*np.sqrt(1./sample_size*res_GC_mean*(1.-res_GC_mean)+(1./(4*sample_size**2))*z**2))
top_err_GC=(1./(1.+(1./(2*sample_size))*z**2))*(res_GC_mean+(1./(2*sample_size))*z**2+z*np.sqrt(1./sample_size*res_GC_mean*(1.-res_GC_mean)+(1./(4*sample_size**2))*z**2))-res_GC_mean

e_GC=plt.errorbar(np.arange(2,(trials+1)*2,2),res_GC_mean,yerr=[below_err_GC,top_err_GC],color='r',ecolor='#990000',fmt='-o',label=r'${\rm WG}$')
plt.xlim(-1,300)


plt.xlabel(r'${\rm Time(s)}$')
plt.ylabel(r'$\%-{\rm CA}_3\to {\rm CA}_1$')    
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0., prop={'size':20})    

#fig,ax=plt.subplots(1,10 ,sharey=True)
#fig.set_size_inches(12.5,5.5)

#fig_GC,ax_GC=plt.subplots(1,10,sharey=True)
#fig_GC.set_size_inches(12.5,5.5)

#p=plt.imshow(res,extent=[0,32,0,32],origin='lower',alpha=1,aspect='1',interpolation='none')

# define the colormap
cmap = plt.get_cmap('PuOr')

# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
plt.savefig(manuscript_path+'Linear_Filters/Figures/SIC-WG-vvp01_2006-4-9_17-29-30.eps',transparent=True,dpi=500)


#define the bins and normalize and forcing 0 to be part of the colorbar!
#bounds = np.linspace(-np.max(np.abs(res)),np.max(np.abs(res)),24)
#idx=np.searchsorted(bounds,0)
#bounds=np.insert(bounds,idx,0)
#norm = BoundaryNorm(bounds, cmap.N)        
#for k in np.arange(10):
#    p=ax[k].imshow(res[:,:,k],extent=[0,32,0,32], aspect=1,interpolation='none',origin='lower',cmap=cmap)
#    p_GC=ax_GC[k].imshow(res_GC[:,:,k],extent=[0,32,0,32], aspect=1,interpolation='none',origin='lower',cmap=cmap)
    #ax[k].set_xlim(0,32)
    #ax_GC[k].set_xlim(0,32)
    #p.set_xticklabels(())
    #p.gca().set_xlim([0,32])
    #p_GC.gca().set_xlim([0,32])
    #plt.colorbar(p,ax=ax)
#print ("cor:",np.mean(np.abs(res-res_GC)))
#fig.subplots_adjust(hspace=0.,wspace=0.) 
#fig_GC.subplots_adjust(hspace=0.,wspace=0.) 
#print res
#plt.tight_layout()
plt.show()