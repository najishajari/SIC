from os import listdir
from matplotlib import pyplot as plt
import numpy as np
from scipy.io import wavfile

#file_path='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/'
#file_path='/Users/Naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/'
tableau20 = np.asarray([(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] )/255.
             

lengths=[8,16,32,64,128]
file_counter=0
def welch_performance(palce_name,song_name):
    counter=0
    final_res=[]
    for length in lengths:
        cluster_path='/agbs/cluster/naji/Linear Filters/Echo/out/'+song_name+'/'+palce_name+'/'+str(length)+'/'
        file_names=(listdir(cluster_path))
        #final_res=np.empty((len(lengths),len(file_names)))
        res=[]
        for file_name in file_names:
            f=open(cluster_path+file_name)
            f.readline()
            res.append([int(file_name[:-4]),float(f.readline())])
            #counter+=1 
            f.close()
        res=np.asarray(res)
        sorted_args=np.argsort(res[:,0])
        res=res[sorted_args]
        #print res
        final_res.append(res)
        #counter+=1
    return np.asarray(final_res)
    
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

final_res=[]
place_names=['Room','Hall']
song_names=['Lacrimosa', 'Winter']
input_dir='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/Sounds/'
fig,ax=plt.subplots(len(place_names),len(song_names))
color_idx=2
for place_idx in range(len(place_names)):
    for song_name_idx in range(len(song_names)):
        res=welch_performance(place_names[place_idx],song_names[song_name_idx])
        #plt.errorbar(res[:,0],res[:,1],color=tableau20[file_counter+2],fmt='-', label=r'$'+str(length)+'$')
        color_idx=2
        for length_idx in range(len(lengths)):
            print res[length_idx].shape
            s=wavfile.read(input_dir+song_names[song_name_idx]+'/Original/'+song_names[song_name_idx]+'.wav')
            half_length_index=s[1].shape[0]/((lengths[length_idx]+1)*1000)
            #print half_length_index,res[length_idx].shape
            #half_length_index=lengths[length_idx].shape[0]/2
            z=1.96
            length=lengths[length_idx]
            suc_prop=res[length_idx][:half_length_index,1]
            below_err=suc_prop-(1./(1.+(1./(2*length))*z**2))*(suc_prop+(1./(2*length))*z**2-z*np.sqrt(1./length*suc_prop*(1.-suc_prop)+(1./(4*length**2))*z**2))
            top_err=(1./(1.+(1./(2*length))*z**2))*(suc_prop+(1./(2*length))*z**2+z*np.sqrt(1./length*suc_prop*(1.-suc_prop)+(1./(4*length**2))*z**2))-suc_prop
            
            if np.min(suc_prop-below_err)<0.5:
                print np.min(suc_prop-below_err),lengths[length_idx],song_names[song_name_idx]
            ax[place_idx,song_name_idx].errorbar(res[length_idx][:half_length_index,0]*500,suc_prop*100.,color=tableau20[color_idx+3],fmt='-', label=r'$'+str(lengths[length_idx])+'$')
            ax[place_idx,song_name_idx].set_xscale('log')
            ax[place_idx,song_name_idx].set_ylim(0.5,1.1)
            ax[place_idx,song_name_idx].set_title(song_names[song_name_idx]+', '+place_names[place_idx])
            ax[place_idx,song_name_idx].set_ylabel(r'${\rm Performance\ \%}$')
            ax[place_idx,song_name_idx].set_xlabel(r'${\rm Window\ Length}$')
            color_idx+=1
        
plt.rcParams['xtick.major.pad']='8'
plt.subplots_adjust(left=None, bottom=0.5, right=None, top=0.9, wspace=None, hspace=None)
leg=ax[0,0].legend(bbox_to_anchor=(0.55, 1.22, 1., .102), loc=3, ncol=5, mode="expand", borderaxespad=0.)
#leg.get_frame().set_alpha(0.)    
#leg.get_frame().set_edgecolor('black')
plt.tight_layout()
plt.savefig('/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Paper/drafts/Spectral_Independence_Criterion_SIC/Figures/echo_fig.eps',transparent=True,bbox_inches='tight')
plt.show()
#plt.xlim(-50,1500)
#plt.ylim(.5,1.1)
#plt.savefig('/is/ei/naji/Dropbox/Room.eps')
#plt.show()

