from matplotlib import pyplot as plt
import numpy as np
res=np.empty((32,18))
cluster_path='/agbs/cluster/naji/Linear Filters/real_data'
from matplotlib.colors import BoundaryNorm

for i in range(32):
    for j in range(18):
        f=open(cluster_path+'/ec016/output/'+str(i)+','+str(j)+'.txt')
        line= f.readline()
        line=line.split(';')
        #res[i,j]= np.abs(np.float64(line)[2])<np.abs(np.float64(line)[3])
        res[i,j]= np.float64(line)[4]>np.float64(line)[5]

print np.mean(res)
fig_1,ax_1=plt.subplots(1,1)
#p=plt.imshow(res,extent=[0,32,0,32],origin='lower',alpha=1,aspect='1',interpolation='none')

# define the colormap
cmap = plt.get_cmap('PuOr')

# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

#p=ax_1.imshow(res,extent=[0,32,0,32],origin='lower', interpolation='none',cmap=cmap)    


# define the bins and normalize and forcing 0 to be part of the colorbar!
bounds = np.linspace(-np.max(np.abs(res)),np.max(np.abs(res)),12)
idx=np.searchsorted(bounds,0)
bounds=np.insert(bounds,idx,0)
norm = BoundaryNorm(bounds, cmap.N)        
p=ax_1.imshow(res,extent=[0,18,0,32], interpolation='none',cmap=cmap,norm=norm)
plt.colorbar(p,ax=ax_1)

#print res
plt.show()