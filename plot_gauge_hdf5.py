import matplotlib.pyplot as plt
import h5py
from pylab import *

class nf(float):
     def __repr__(self):
         str = '%.1f' % (self.__float__(),)
         if str[-1]=='0':
             return '%.0f' % self.__float__()
         else:
             return '%.1f' % self.__float__()

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
       
if __name__=="__main__":
 
    outdir='New_AS_output_pgftri_fullwave_pure'
    hdffile='Alaska_Jagurs_post_simul.h5'
    #sim=sys.argv[1]
    fop=h5py.File(outdir+'/'+'waveform_9g_233.h5','r')    
    frw=h5py.File(outdir+'/'+'waveform_11g_233.h5','r')
    gname= {101:'Alitak',7101:'Alitak',  801:'Kodiak', 7801:'Kodiak',
	 901:'Langara', 201:'Seward', 601:'Yakutat', 7601:'Yakutat',
	701:'Unalaska', 46402:46402,46403:46403, 46407:46407, 
	46408:46408, 46409:46409, 46410:46410, 46411:46411, 46419:46419}

    plot_window={101:[50,200], 7101:[50,220], 201:[10,200],\
    601:[20,150], 7601:[10,150], 801:[10,200], 7801:[10,200],\
    901:[20,180], 46402:[30,120], 46403:[30,80], 46408:[50,150],\
    46409:[0,40], 46410:[0,60], 46411:[150,250], 46419:[100,200]}	
    
    glist=[46402, 46403, 46408, 46410, 46419, 46411]#, 7101, 201, 7601, 7801, 901 ]  
 
    print glist
    fig=plt.figure(figsize=(8, 6))
    n=len(glist)
    nc=2
    nr=int(ceil(n/float(nc)))
    plt.clf()
    fig.subplots_adjust(hspace=.4, wspace=.2)


    for k in range(len(glist)):
	ax=fig.add_subplot(nr,nc,k+1)
	plt.title('%s'%gname[glist[k]])
	tw=plot_window[glist[k]]
	
	who=frw["observed/%s/data"%glist[k]][...]
	time=frw["time/%s/data"%glist[k]][...]
	#ma = (t >= tw[0]) & (t <= tw[1])
	#time = t[ma]
	#who = who[ma]    
	plt.plot(time, who, '--b',  label='observed')	
	
	whc=frw["computed/%s/data"%(glist[k])][...]
	  
	plt.plot(time, whc, '-r', label='non-optimal')
	whc=fop["computed/%s/data"%(glist[k])][...]
 
	plt.plot(time, whc, '-g', linewidth=2, label='optimal')
	#plt.grid(True)

	simpleaxis(ax)
	ti=time[0]
	tf=time[-1]
	xtk=np.linspace(ti,tf,5)
	plt.xticks(xtk.round())
	ytk=np.linspace(round(whc.min(),2),round(whc.max(),2),5)
	ytk=[round(y,2) for y in ytk]
	plt.yticks(ytk)

	    	
	if glist[k]==46411:
	    plt.legend(loc=3, borderaxespad=0., fontsize=10)


    #plt.legend(loc=0, bbox_to_anchor=(1., .65), borderaxespad=0., fontsize=12)
     
    fig.text(0.5, 0.02, 'time (min)', ha='center', va='center')
    fig.text(0.06, 0.5, 'water height (cm)', ha='center', va='center', rotation='vertical')	
    #fig.suptitle(title)	
    plt.savefig('compare_tg',bbox_inches='tight', pad_inches=0.1)


