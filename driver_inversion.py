'''
This is written by Dr. Jakir Hossen
Date: Feb, 2018
'''
import sys, os
import h5py
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal as sg
import math as mt
from scipy.io import netcdf

# These are associated with source inversion: LSQ and GFTRI

from Support_Inversion import Fault_parameters, Data, LSQ_method, \
    Pure_gftri, Output



#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Main Program &&&&&&&&&&&&&&&&&&&&&&&&

class main(Fault_parameters, Data, LSQ_method, Pure_gftri, Output):
    '''
    The derivedClass main inherits from different 
    BaseClass: Fault_parameters, Data, LSQ_method, GFTRI, AS, Output
    '''
    # constructor
    def __init__(self, glist=None,
		    plt_glist=None,
		    hdf_pfile=None, 
		    hdf_gfile=None, 
		    hdf_ofile=None, 
		    dem_file = None,
		    inv_time=None,
		    gloc=None, 
		    remove_grid=None, 
		    waveoption=None, 
		    stn_name=None,
		    tw=None,
		    method = None,
		    interp_rate=None,
		    source_dim=None):

	Fault_parameters.__init__(self, hdf_pfile, remove = remove_grid)

	Data.__init__(self, glist, 
			    hdf_gfile, 
			    hdf_ofile,
			    waveoption=waveoption, 
			    interp_rate=interp_rate,
			    tw=tw)
			    
	self.inv_time=inv_time 
	self.method=method
	self.dem_file=dem_file
	print 'Total inversion time: %s'%inv_time
	
	A, b, time_param=self.compute_system()	# Read GFs and observation
	self.A_in=A
	self.b_in=b
	self.time_parin=time_param
	self.gglist=glist
	self.update(A, b, time_param)
	
	# ************** Least squares method ************************
	if method == 'lsq':
	    LSQ_method.__init__(self)

	elif method == 'pgftri':		  
	    tshift=self.shift_param()		# time shift for pure GFTRI
	    Pure_gftri.__init__(self, gloc=gloc, ind_shift=tshift)
	else:
	    print 'Please, provide method: lsq, pgftri (pure) or cgftri (complete)'
	    sys.exit(0)
	    
			    
	Output.__init__(self, gauge_name=stn_name, 
			    source_dim=source_dim)

    def update(self, A=None, b=None, time_param=None, glist=None):

	if glist is not None:
	    self.glist=glist

	if A is None:
	    self.A, self.b, self.time_param=self.compute_system()
	else:    
	    self.A=A
	    self.b=b
	    self.time_param=time_param

    
    
    # ================ create A and b for Ax=b  ========================		    
    def compute_system(self):
	
	A_computed, b_observed, time_param=self.read_hdf5file()
	A=np.array(A_computed)
	ba=np.array(b_observed)
	
	return A, ba, time_param

    def compute_modelerror(self, amp):
	A=self.A_in
	ba=self.b_in
	gglist=self.gglist
	bb=ba[:,None]
	ym=amp[:,None]
	Am=np.mat(A)
	yy=np.array(Am*ym)	
	misfit_err=np.linalg.norm(bb-yy)/float(len(gglist))
	print 'error:',misfit_err 
	    
	return misfit_err
	
    # ======== calcuate simulated waveform using solution: amplitude ====
    def compute_waveform(self):

	if self.method == 'lsq': 	# least squares method
	    
	    amp=self.LSQ_inversion(option='damping')
	    
	elif self.method == 'pgftri':	# partial GFTRI    
	    amp=self.amp_pure_GFTRI()
	    
	else:
	    print 'Please, provide method: lsq, pgftri (pure) or cgftri (complete)'
	    sys.exit(0)
	A=self.A
	ym=amp[:,None]
	Am=np.mat(A)
	y=np.array(Am*ym)	    

	return amp, y
	
    # ====================== Time shift for Pure GFTRI =================
    def shift_param(self):
	lon,lat, elev = self.read_dem(self.dem_file)
	ind_shift=[]
	cols=self.num_sbr
	for i in range(cols):
	    cntr=self.sbfs[i][0]
	    ix=np.argmin(abs(lon-cntr[0]))
	    iy=np.argmin(abs(lat-cntr[1]))
	    elevi=np.mean(abs(elev[iy-5:iy+5,ix-5:ix+5]))
	    trad=2000./np.sqrt(9.8*elevi)	    
	    indnt=0#int(round(trad/self.interp_rate))
	    ind_shift.append(indnt)
	return ind_shift



    # ================= inversion without AS method =====================	
        
    def sourceInversion(self, outdir):

	amp, y = self.compute_waveform()

	ext=int(self.inv_time/60.0)  # ext- only to save files 
	ng=len(self.glist)    

	self.plot_slip(amp, filename=outdir+'/'+'slip_%sg'%ng)

	self.create_full_source(amp, \
	    filename=outdir+'/'+'source_model_%sg_%s'%(ng, ext))	

	self.plot_gauge(tg=y, obs=self.b, time_param=self.time_param, 
	    filename=outdir+'/'+'waveform_%sg_%s'%(ng, ext))


# &&&&&&&&&&&&&&&&&&&&&&   main program  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	
if __name__=="__main__":
    
	if len(sys.argv) <= 2:
	    print 'Please, provide method: lsq, pgftri, or cgftri'
	    print 'and waveoption: fullwave, fwave (First wave), tw (time window)'
	        
	# ================= options and directory ======================
	

	waveoption='tw' 	#sys.argv[1]
	method='pgftri'	#sys.argv[2]

	data_dir='./Input/'

	# =========================== hdf5 files ======================= 
	hdf_pfile=data_dir+'fault_param_25km_Jagurs.h5'
	hdf_gfile=data_dir+'Alaska_GFs_25km_Jagurs.h5'
	rem_grd=[0, 10,20, 30, 40, 50, 9,19,29,39,49,59]

	hdf_ofile=data_dir+'observation_Alaska.h5'

	glist=[46402, 46403, 46408, 46407, 46410, 46411, 46419]

	outdir='Output'
	
	if not os.path.exists(outdir):
	    os.makedirs(outdir) 

	# ================ gauge information ============================
	
	gauge_name= {101:'Alitak',7101:'Alitak',  801:'Kodiak', 7801:'Kodiak',
	     901:'Langara', 201:'Seward', 601:'Yakutat', 7601:'Yakutat',
	    701:'Unalaska', 46402:46402,46403:46403, 46407:46407, 
	    46408:46408, 46409:46409, 46410:46410, 46411:46411, 46419:46419}

	time_window={101:[80,200], 7101:[60,150], 201:[60,150],\
	601:[100,150], 7601:[60,150], 801:[80,200], 7801:[80,200],\
	901:[80,200], 46402:[60,120], 46403:[20,80], 46408:[80,140],\
	46409:[1,30], 46410:[10,70], 46407:[150,210], 46411:[180,240], 46419:[110,170]}	
	
	plot_window={101:[110,180], 7101:[110,180], 201:[100,160],\
	601:[90,150], 7601:[90,150], 801:[110,200], 7801:[110,200],\
	901:[120,180], 46402:[60,120], 46403:[30,80], 46408:[100,150],\
	46409:[0,40], 46410:[0,60], 46411:[200,240], 46419:[130,180]}	

	
	ng=len(glist)

	stg=np.loadtxt(data_dir+'stations.txt', skiprows=1) 
	gloc={}
	for i in range(len(stg)):
		gloc[int(stg[i,2])]=(stg[i,1],stg[i,0])  # (lon, lat)
	#print gloc
	# =====================  time information ======================
	
	final_time=4*3600	


	dem_file=data_dir+'gebco_Alaska.grd'
	# ================== main source inversion =====================
	
	# ---------------- removed source patches ------------------
		
	#print rem_grd
	#sys.exit(0)
	src_dim=[208, 214, 54, 58]
	mobj = main(glist=glist, 
		    hdf_pfile = hdf_pfile, 
		    hdf_gfile = hdf_gfile, 
		    hdf_ofile = hdf_ofile, 
		    dem_file = dem_file,
		    inv_time = final_time, 
		    waveoption = waveoption, 
		    stn_name = gauge_name,
		    tw = time_window,
		    method = method,
		    remove_grid=rem_grd,
		    gloc=gloc,
		    interp_rate=30,
		    source_dim=src_dim)	

	print 'number of considered source grid:', mobj.num_sbr
	
	# ***************** source Inversion *********************

	mobj.sourceInversion(outdir)
	
	
	
	


        
