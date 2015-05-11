from scipy.io import wavfile
from matplotlib import pyplot as plt
def split_wave(num=11,song_name=None,input_path=None,output_path=None,file_name=None):
    """function for splitting a wave file into desired number of pieces"""
    print input_path+file_name+'.wav'
    print output_path
    s=wavfile.read(input_path+file_name+'.wav')
    raw_len=s[1].shape[0]/(num*s[0])
    print s[1].shape[0],s[0]
    #for i in range(num+1):
        #wavfile.write(output_path+song_name+'_'+str(i+1),s[0],s[1][raw_len*i*s[0]:raw_len*(i+1)*s[0]])

song_name='Lacrimosa'
file_name='Lacrimosa'

#song_name='Winter'
#file_name='Winter'

#input_dir='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/Sounds/'+song_name+'/Experiments/Hall/3/'
#output_dir='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/Sounds/'+song_name+'/Experiments/Hall/3/Segments/'

input_dir='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/Sounds/'+song_name+'/Original/'
output_dir='/is/ei/naji/Dropbox/Winter Semster 2014/Master Thesis/Programming/Echo Test/Sounds/'+song_name+'/Original/Segments/'

split_wave(num=6,song_name=song_name,input_path=input_dir,file_name=file_name,output_path=output_dir)
        
    

