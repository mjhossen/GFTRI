import sys, os
import h5py
import numpy as np
from pylab import ceil
from pyproj import Geod
from scipy import interpolate
from scipy.interpolate import griddata
from scipy.io import netcdf
from matplotlib.colors import Normalize
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Polygon
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import rgb2hex
import math as mt


# &&&&&&&&&&&&&&&&&&&&&&&&& Fault Information &&&&&&&&&&&&&&&&&&&&&&&&&&

class Fault_parameters:

    def __init__(self, hdf_file, remove=None):
	if remove is None:
	    self.remove=[9999]
	else:
	    self.remove=remove
	     
	sbfs, sbf_ind=self.read_parameter(hdf_file)  # read fault param 
	self.sbfs=sbfs
	self.sbf_ind=sbf_ind
	self.num_sbr=len(sbfs)
	 
    def read_parameter(self, hdf_file):
	
	f=h5py.File(hdf_file,'r')
	self.hypo=f["fault/parameters/hypo_center"][...]
	fault_length=f["fault/parameters"].attrs['fault_length']
	fault_width=f["fault/parameters"].attrs["fault_width"]
	self.nx=f["fault/parameters"].attrs["sbf_nx"]
	self.ny=f["fault/parameters"].attrs["sbf_ny"]
	self.sbf_dx=f["fault/parameters"].attrs["sbf_dx"]
	self.sbf_dy=f["fault/parameters"].attrs["sbf_dy"]
        self.strike=f["fault/parameters"].attrs["strike"]
        self.rake=f["fault/parameters"].attrs["rake"]
        self.dip=f["fault/parameters"].attrs["dip"]
	#print nx, ny
	cols=self.nx*self.ny
	sbfs={}
	sbf_ind={}
	k=0
	for i in range(cols):
	    if i not in self.remove:
		sbf_ind[k]=i
		sbfs[k]=[f["fault/sbf%s/center"%i][...],f["fault/sbf%s/polygon"%i][...]]
		k=k+1
	f.close()

	return sbfs, sbf_ind

    def distance_from_hypo_formula(self):
	from obspy.geodetics.base import gps2dist_azimuth as dist_cal
	hypo=self.hypo
	dist=[]
	for i in range(len(self.sbfs)):
	    cntr=self.sbfs[i][0]
	    lon1=cntr[0]
	    lat1=cntr[1]
	    lon2=hypo[0]
	    lat2=hypo[1]
	    
	    d1=dist_cal(lat1, lon1, lat2, lon2)
	    d=d1[0]/1000.0
	    if i==306:
		print i, lon1, lat1, lon2, lat2, d
	    dist.append(d)
	return dist

# &&&&&&&&&&&&&&&&&&&&&&&&&& Initial Data &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class Data:
    '''
    To create data structure from simulation results
    '''
    def __init__(self, glist, 
		hdf_gfile, 
		hdf_ofile, 
		weights=None, 
		sbf_dist=None, 
		interp_rate=None,
		waveoption=None,
		tw=None):

	self.glist=glist
	self.gfile=hdf_gfile
	self.ofile=hdf_ofile  
	self.interp_rate=interp_rate
	self.tw=tw
	self.waveoption=waveoption
	
    def read_hdf5file(self,inv_time=None,  \
			glist=None, plot_tg=False):
	import peakutils	# this is for peak detection in the waveform    
	from scipy.signal import find_peaks
	
	hdf_gfile=self.gfile
	hdf_ofile=self.ofile
	interp_rate=self.interp_rate
	
	if glist is None:
	    glist=self.glist

	if inv_time is None:
	    inv_time=self.inv_time

	cols=self.num_sbr
	sbf_ind=self.sbf_ind

	f=h5py.File(hdf_ofile,'r')
	obs_plot=False
	b_observed=[]

	time_param={}
	# Read observation data given in the glist. But for this we choose two options: 
	#1. provide time window for individual observations or 2. consider the same time window for all
	# One important think is that observation does not depend on the choice of source grid point but it is important
	# while reading GFs.
	
	for k in range(len(glist)):
	    observed=f["observed/%s/data"%glist[k]][...]
	    to=observed[:,0]
	    who=observed[:,1]
	    dt=interp_rate
	    t_0, t_f=1.0, inv_time
	    # ============= determine the time window ==================
	    
	    indf= np.argmin(abs(to-inv_time))

	    if self.waveoption=='fullwave':
		indmax=indf
		indmin=1			# from initial time
		#obs_plot=True
		t0=to[indmin]
		tf=to[indmax]
		if tf>t_f:
		    tf=t_f
	    elif self.waveoption=='fwave':
		# for TG we use peak detection to find the first peak
		# and DART we juse find the index for maximum

		ind, _=find_peaks(who, height=0)
		#print indexes
		imax=np.argmin(abs(who-max(who[ind])))
		    
		#print imax
		dst=to[10]-to[9]		
		ntr=int(30*60/dst)  	# add 20 min more data to the maximum peak
		indmax=min(indf,imax+2*ntr)	# constraints within final time
		indmin=max(1,imax-2*ntr)
		#indmin=1
		obs_plot=True
		t0=to[indmin]
		tf=to[indmax]
		if tf>t_f:
		    tf=t_f

	    elif self.waveoption == 'tw':
		t0, tf=self.tw[glist[k]]
		t0=t0*60
		tf=tf*60
		# only use the part which contains the first wave.
		ma=to<=t0
		who[ma]=0
		ma=to>=tf
		who[ma]=0

		print glist[k], t0, tf
		if tf>t_f:
		    tf=t_f
		obs_plot=False
		
	    else:
		print 'please, provide wave option: First wave, Full wave or tw'
	    if plot_tg:
		indmax=indf 		# for upto final time
		indmin=0
		obs_plot=False
	    
	    # ================= Automatic Time window ===================
	    print glist[k], t0, tf
	    nt=int(round((tf-t0))/dt)+1
	    t_new=np.linspace(t0,tf,nt)
	    f_obs=interpolate.interp1d(to,who)
	    obs_new=f_obs(t_new) 
	    
	    if obs_plot:
		plt.clf()
		plt.plot(t_new, obs_new)
		plt.savefig('fig_%s'%glist[k])

	    time_param[glist[k]]=[t_new[0], t_new[-1], nt, dt]
	    
	    for ii in range(len(obs_new)):
		b_observed.append(obs_new[ii])
	f.close()	    
	
	fg=h5py.File(hdf_gfile,'r')
	A_computed=np.zeros((len(b_observed),cols))
	#print A_computed.shape

	# Read the GFs corresponding to the accounted source grid point. 
	# Use source index stored in sbf_ind array
	
	for i in range(cols):
	    #print sbf_ind[i]
	    Nobsp=0
	    for j in range(len(glist)):
		computed=fg["computed/gf%s/%s/data"%(sbf_ind[i],glist[j])][...]
		tc=computed[:,0]
		whc=computed[:,1]
		tc=tc[::2]
		whc=whc[::2]
		t_0, t_f, nt, dt=time_param[glist[j]]
		t_new=np.linspace(t_0, t_f, nt)
		#print t_f, tc[-1]
		f_com=interpolate.interp1d(tc,whc)
		c_new=f_com(t_new) 
		A_computed[Nobsp:Nobsp+nt,i]=c_new[:]
		Nobsp+=nt

	fg.close()
	
	return A_computed, b_observed, time_param

    def read_dem(self, filename):

	ncf = netcdf.netcdf_file(filename,'r')
	if hasattr(ncf,'cellsize'):
	    # Set UTM grid params from netcdf header
	    cellsz = ncf.cellsize[0]
	    nrows  = ncf.nrows[0]
	    ncols  = ncf.ncols[0]
	    xll    = ncf.xllcorner[0]   # Easting of lower left corner
	    yll    = ncf.yllcorner[0]  # Northing of lower left corner
	    xur    = xll + (ncols-1)*cellsz
	    Yur    = yll + (nrows-1)*cellsz
	    x = np.linspace(ncf.xllcorner[0],ncf.xllcorner[0] + \
				(ncols-1)*cellsz,ncols)
	    y = np.linspace(ncf.yllcorner[0],ncf.yllcorner[0] + \
				(nrows-1)*cellsz,nrows)
	    zone   = ncf.zone[0]
	    zone = 54
	    zdat   = np.flipud(ncf.variables['elevation'][:])
	    # Made from GMT?
	elif 'x_range' in ncf.variables:
	    xrng = ncf.variables['x_range']
	    yrng = ncf.variables['y_range']
	    zrng = ncf.variables['z_range']
	    dxdy = ncf.variables['spacing']
	    dim = ncf.variables['dimension']
	    lnll = xrng[0]
	    lnur = xrng[1]
	    ltll = yrng[0]
	    ltur = yrng[1]
	    nx=dim[0]
	    ny=dim[1]
	    zdat = ncf.variables['z'][:].reshape(ny,nx)
	    zdat = np.flipud(zdat)
	    x = np.linspace(lnll,lnur,nx)
	    y = np.linspace(ltll,ltur,ny)
	elif ncf.Conventions == 'COARDS/CF-1.0':
	    if 'x' in ncf.variables:
		x = ncf.variables['x'][:]
		y = ncf.variables['y'][:]
	    elif 'lon' in ncf.variables:
		x = ncf.variables['lon'][:]
		y = ncf.variables['lat'][:]
	    else:
		print 'No x or lon found in COARDS netcdf file %s' %  filename
		return(-1)
	    zdat=ncf.variables['z'][:]
	else:
	    print 'File %s is not  ANUGA or GMT' % filename
	    return(-1)
	ncf.close()
	
	return (x,y,zdat)

# &&&&&&&&&&&&&&&&&&&&&&&&&&& Partial GFTRI &&&&&&&&&&&&&&&&&&&&&&&&&&&&

class Pure_gftri:

    def __init__(self, gloc=None, ind_shift=None):
	
	if gloc is not None:
	    self.gloc=gloc
	else:
	    print "Warning! gloc could be a problem"
	self.ind_shift=ind_shift    
   
    def distance_obs_sr(self, cntr, gpos):
	#from obspy.geodetics.base import gps2dist_azimuth as dist_cal
	from geopy.distance import geodesic as dist_cal	
	lon1=cntr[0]
	lat1=cntr[1]
	lon2=gpos[0]
	lat2=gpos[1]	
	d1=dist_cal((lat1, lon1), (lat2, lon2))
	#d=d1[0]/1000.0

	return d1	

    # ========================== Pure GFTRI ============================	
    def amp_pure_GFTRI(self):

	'''
	In this method we apply time reverse imaging method to calculate amplitude 
	of each source 
	'''
	
	glist=self.glist
	ng=len(glist)
	cols=self.num_sbr
	amp=np.zeros(cols)

	for k in range(cols):
	    gf=self.A[:,k]
	    obs=self.b
	    cntr=self.sbfs[k][0]
	    if self.ind_shift is None:
		ind_nt=0
	    else:
		ind_nt = self.ind_shift[k]
	    Nobs=0
	    ngc=0
	    for i in range(ng):
		t0,tf,nt,dt=self.time_param[glist[i]]
		g=gf[Nobs:Nobs+nt]
		fcong=np.convolve(g,g[::-1])


		#print gnorm
		o=obs[Nobs:Nobs+nt]
		Nobs+=nt
		gpos=self.gloc[glist[i]]
		#print gpos, cntr
		dist_gsr=self.distance_obs_sr(cntr,gpos)
		if dist_gsr<= 1800.:
		    atn_factor=2.  
		else:
		    atn_factor=2.#.8
		    
		#print glist[i], atn_factor, dist_gsr
		
		gnorm=fcong[nt-1-ind_nt]*atn_factor	
		

		if gnorm==0:
		    print 'At station: %s, gnorm is %s'%(glist[i], gnorm)
		    continue		
			
		if dist_gsr <= 50.0:
		    #print 'Neglecting gauge:%s located within %d km or gnorm is zero'%(glist[i], dist_gsr)
		    continue
		else:
		    ngc=ngc+1
		    
		fcon=np.convolve(g,o[::-1])
		amp_val=fcon[nt-1-ind_nt]/gnorm # amp=(1/G^TG)G^Td(-t)
		if ngc==1:
		    cum_fcon=amp_val
		else:
		    cum_fcon+=amp_val

	    cum_fcon=cum_fcon/float(ngc)  # averaging as TRI is self averaging
	    amp[k] = cum_fcon

	return amp



# &&&&&&&&&&&&&&&&&&&& Least Squares inversion &&&&&&&&&&&&&&&&&&&&&&&&&

class LSQ_method:
    '''
    This is to find best damping and smoothing parameters 
    '''
    def __init__(self):
	pass

    def LSQ_inversion(self, option=None):
	 	
	self.b=self.b[:,None]
	A=self.A
	b=self.b

	# ------------ smoothing matrix---------------------------------
	
	P=self.smooth_matrix_sbfs()
	print P

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.spy(P, markersize=5)
	plt.savefig('smooth_matrix_rad')
	num_tw=1
	prow,pcol=P.shape	
	I=np.identity(pcol*num_tw)
	z=np.zeros((pcol*num_tw,1), float)
	    
	P_mw=np.zeros((prow,pcol*num_tw))	
	zz=np.zeros((prow,1), float)
	for ii in range(num_tw):
	    P_mw[:,ii*pcol:(ii+1)*pcol]=P[:,:]	
	    
	print P_mw.shape
	
	    
	# ---------------- estimate weights and apply on A and b --------------

	weights, sval=self.estimate_weights()
	 
	A_weighted, b_weighted=self.weighted_matrix(weights)
	A=np.array(A_weighted)
	b=np.array(b_weighted)
	print 'size: A,b: ', A.shape,b.shape
	# ---------------- Estimate parameter for damping and smoothing---------------
	print 'option:',option
	# ---------------- Estimate parameter for damping and smoothing---------------
	if option=='damping':
	    lmd_d=self.damping_grid_search(A,b, I, z,iter=2)
	    print 'parameter:',lmd_d
	    B=np.vstack((A, lmd_d*I))
	    d=np.vstack((b,z))
	elif option=='smoothing':
	    lmd_d, lmd_s=self.smoothing_grid_search(A,b, P_mw, zz,iter=2)
	    B=np.vstack((A, lmd_d*I, lmd_s*P_mw))
	    d=np.vstack((b,z,zz))
	    print 'parameter: damping-%s and smoothing-%s'%(lmd_d, lmd_s)

	elif option=='both':
	    lmd_d, lmd_s=self.grid_search(A,b, I,z, P_mw, zz,iter=3)
	    B=np.vstack((A, lmd_d*I, lmd_s*P_mw))
	    d=np.vstack((b,z,zz))
	    print 'parameter: damping-%s and smoothing-%s'%(lmd_d, lmd_s)

	else:
	    print 'no regularization constraints is used'
	    B=A
	    d=b
		
	# ----------------- solve the system with estimate parameter------------------
	
	amp=np.linalg.lstsq(B, d[:,0])[0]
	return amp

	
    def smooth_matrix_sbfs(self):
	from obspy.geodetics.base import gps2dist_azimuth as dist_cal
	cols=self.num_sbr
	print cols
	P=np.zeros((cols, cols),float)
	nx=self.nx
	ny=self.ny
	nmax=max(nx,ny)+1
	sbfs=self.sbfs
	for i in range(cols):
	    lon1, lat1=sbfs[i][0]
	    lb=max(0,i-nmax)
	    ub=min(cols, i+nmax)
	    for l in range(lb,ub):
		lon2, lat2=sbfs[l][0]
		d1=dist_cal(lat1, lon1, lat2, lon2)
		d=d1[0]/1000.0
		if d <=self.sbf_dx+2:
		    P[i,l]=1
		    if i==l:
			P[i,i]=-4
	return P
	
    def estimate_weights(self):
	A=self.A
	b=self.b
	glist=self.glist
	time_param=self.time_param
	c=np.linalg.lstsq(A, b[:,0]) 
	slip=c[0]
	sval=c[3]
	ax=np.dot(A,slip)
	fit=ax[:,None]
	res=b-fit
	wgt={}
	Nobsp=0
	for i in range(len(glist)):
	    t0, tf, nt, dt=time_param[glist[i]]
	    wgt[glist[i]]=1.0/np.std(res[Nobsp:Nobsp+nt,0])
	    Nobsp+=nt
	return wgt, sval
	
    def weighted_matrix(self, weights):
	A=self.A
	b=self.b
	glist=self.glist
	time_param=self.time_param
	rows, cols=A.shape
	A_weighted=np.zeros_like(A)
	for i in range(cols):
	    Nobsp=0
	    for j in range(len(glist)):
		t0,tf, nt, dt=time_param[glist[j]]		    
		A_weighted[Nobsp:Nobsp+nt,i]=A[Nobsp:Nobsp+nt,i]*weights[glist[j]]
		Nobsp+=nt
		
	b_weighted=np.zeros_like(b)
	
	Nobsp=0
	for j in range(len(glist)):
	    t0,tf,nt, dt=time_param[glist[j]]		    
	    b_weighted[Nobsp:Nobsp+nt,0]=b[Nobsp:Nobsp+nt,0]*weights[glist[j]]
	    Nobsp+=nt
	print Nobsp
	
	return A_weighted, b_weighted
	    
    def grid_search(self, A,b, I,z, P,zz,iter=None):
	time_param=self.time_param
	glist=self.glist

	if iter==None:
	    iter=2
	#lmd1=np.logspace(-4,2,10)
	#lmd1=np.logspace(-3,2,10)
	lmd1=np.linspace(0.01,10,8)
	lamda_damp=lmd1#.tolist()+lmd2.tolist()
	lamda_smooth=lamda_damp#lmd1.tolist()+lmd2.tolist()
	print lamda_damp
	
	best_model=[]
	for k in range(iter):
	    model=[]
	    if k!=0:
		print 'iter, parameters:', k, lmd_d, lmd_s
		lamda_damp=np.linspace(lmd_d-.5*lmd_d, lmd_d+0.5*lmd_d,5)
		lamda_smooth=np.linspace(lmd_s-.5*lmd_s, lmd_s+0.5*lmd_s,5)	
		if lmd_d < 1:
		    lamda_damp=np.linspace(lmd_d-.9*lmd_d, lmd_d+1,5)
		if lmd_s<1:
		    lamda_smooth=np.linspace(lmd_s-.9*lmd_s, lmd_s+1,5)	  
	    for lmd_d in lamda_damp:
		for lmd_s in lamda_smooth:
		    B=np.vstack((A, lmd_d*I, lmd_s*P))
		    d=np.vstack((b,z,zz))	
		    c=np.linalg.lstsq(B, d[:,0])
		    slip=c[0]
		    ym=slip[:,None]
		    Am=np.mat(A)
		    #print Am.shape, ym.shape
		    fit=np.array(Am*ym)   
		    res=b-fit
		    chisq = np.zeros(len(glist))
		    Nobsp=0
		    for i in range(len(glist)):
			t0, tf, nt, dt=time_param[glist[i]]
			tr_res = res[Nobsp:Nobsp+nt]
			tr_res-=np.mean(tr_res) 
 			chisq[i] = np.sum(tr_res**2)
			Nobsp+=nt			
		    print 'chisq, Nobs, p1, p1',np.sum(chisq), Nobsp, lmd_d, lmd_s
		    chi=np.sum(chisq)
		    model.append([lmd_d,lmd_s,chi, c[2]])
		    if chi>1.5*Nobsp:
			break
	    model_a=np.array(model)
	    ind_best=np.argmin(abs(model_a[:,2]-Nobsp))
	    lmd_d=model_a[ind_best,0]
	    lmd_s=model_a[ind_best,1]
	    chi=model_a[ind_best,2]
	    rank=model_a[ind_best,3]
	    best_model.append([lmd_d,lmd_s,chi])
	    print lmd_d, lmd_s, rank, chi
	#-------------complete all iterations -------------------------    
	best_a=np.array(best_model)
	ind_best=np.argmin(abs(best_a[:,2]-Nobsp))
	lmd_d=best_a[ind_best,0]
	lmd_s=best_a[ind_best,1]
	chi=best_a[ind_best,2]
 
	return lmd_d,lmd_s

    def damping_grid_search(self, A,b, I,z,iter=None):
	time_param=self.time_param
	glist=self.glist

	if iter==None:
	    iter=2
	lmd1=np.logspace(-4,1,10)
	lamda_damp=lmd1
	#lamda_damp=np.linspace(.0001,10,50) # for refine search without condition on sv
	print lamda_damp
	
	for k in range(iter):
	    model=[]
	    if k!=0:
		print 'iter, parameters:', k, lmd_d
		lamda_damp=np.linspace(lmd_d-.8*lmd_d, lmd_d+0.9*lmd_d,5)
	    for lmd_d in lamda_damp:
		    B=np.vstack((A, lmd_d*I))
		    d=np.vstack((b,z))	
		    c=np.linalg.lstsq(B, d[:,0])  # for refine search without condition on sv
		    #c=np.linalg.lstsq(B, d[:,0],rcond=.005)
		    slip=c[0]
		    ym=slip[:,None]
		    Am=np.mat(A)
		    #print Am.shape, ym.shape
		    fit=np.array(Am*ym)   
		    res=b-fit
		    chisq = np.zeros(len(glist))
		    Nobsp=0
		    for i in range(len(glist)):
			t0, tf, nt, dt=time_param[glist[i]]
			tr_res = res[Nobsp:Nobsp+nt]
			tr_res-=np.mean(tr_res) 
 			chisq[i] = np.sum(tr_res**2)
			Nobsp+=nt			
		    chi=np.sum(chisq)
		    model.append([lmd_d, chi, c[2]])
	    model_a=np.array(model)
	    x=model_a[:,0]
	    chi=model_a[:,1]
	    plt.clf()
	    plt.plot(x,abs(chi-Nobsp))
	    plt.savefig('best_model_damping_%s'%k)
	    ind_best=np.argmin(abs(chi-Nobsp))
	    lmd_d=model_a[ind_best,0]
	    cval=model_a[ind_best,1]
	    rank=model_a[ind_best,2]
	    
	    print 'lamda:%s, rank:%s, chisquare:%s'%(lmd_d, rank, cval)
	return lmd_d

    def smoothing_grid_search(self, A,b, P,z,iter=None):
	time_param=self.time_param
	glist=self.glist

	if iter==None:
	    iter=2
	lmd1=np.logspace(-4,1,10)
	#lmd2=np.linspace(0.2,50,8)
	lamda_smooth=lmd1#.tolist()+lmd2.tolist()

	print lamda_smooth
	lmd_d=0.0
	
	for k in range(iter):
	    model=[]
	    if k!=0:
		print 'iter, parameters:', k, lmd_s
		lamda_smooth=np.linspace(lmd_s-.8*lmd_s, lmd_s+0.9*lmd_s,5)
	    for lmd_s in lamda_smooth:
		    B=np.vstack((A, lmd_s*P))
		    d=np.vstack((b,z))	
		    c=np.linalg.lstsq(B, d[:,0])
		    slip=c[0]
		    ym=slip[:,None]
		    Am=np.mat(A)
		    #print Am.shape, ym.shape
		    fit=np.array(Am*ym)   
		    res=b-fit
		    chisq = np.zeros(len(glist))
		    Nobsp=0
		    for i in range(len(glist)):
			t0, tf, nt, dt=time_param[glist[i]]
			tr_res = res[Nobsp:Nobsp+nt]
			tr_res-=np.mean(tr_res) 
 			chisq[i] = np.sum(tr_res**2)
			Nobsp+=nt			
		    chi=np.sum(chisq)
		    model.append([lmd_s, chi, c[2]])
	    model_a=np.array(model)
	    x=model_a[:,0]
	    chi=model_a[:,1]
	    plt.clf()
	    plt.plot(x,abs(chi-Nobsp))
	    plt.savefig('best_model_smoothing_%s'%k)
	    ind_best=np.argmin(abs(chi-Nobsp))
	    lmd_s=model_a[ind_best,0]
	    cval=model_a[ind_best,1]
	    rank=model_a[ind_best,2]
	    
	    print 'lamda:%s, mu:%s, rank:%s, chisquare:%s'%(lmd_d, lmd_s, rank, cval)
	return lmd_d, lmd_s

# &&&&&&&&&&&&&&&&& Result plot &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

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
    
class Output:
    '''
    plot the output or create output file  
    '''
    def __init__(self, gauge_name=None, source_dim=None):
	self.gauge_name=gauge_name
	if source_dim is not None:
	    self.source_dim=source_dim 	 

    def create_hdffile(self, tg=None, obs=None, time_duration=None, hdf_file=None, tw_plot=None, time_param=None):
	glist=self.glist
	if time_param is None:
	    time_window=self.time_param
	else:
	    time_window=time_param
	n=len(glist)
	Nobsp=0

	f=h5py.File(hdf_file,'w')
	    
	for i in range(n):
	    tw=time_window[glist[i]]
	    time=np.linspace(tw[0],tw[1],tw[2])/60.
	    obs_h=obs[Nobsp:Nobsp+tw[2]]*100 
	    tg_h=tg[Nobsp:Nobsp+tw[2]]*100   
	    Nobsp=Nobsp+tw[2]
	    if tw_plot is not None:
		t0,tf=tw_plot[glist[i]]
	    	mask=(time>=t0) & (time<=tf)
	    	t_new=time[mask]
	    	o_new=obs_h[mask] 
		tg_new=tg_h[mask]
	    else:
		t_new=time
	        o_new=obs_h
             	tg_new=tg_h
	       

	    owf=zip(t_new, o_new)
	    oset=f.create_dataset("observed/%s/data"%(glist[i]),data=owf)
	    cwf=zip(t_new, tg_new)
	    cset=f.create_dataset("computed/%s/data"%(glist[i]),data=cwf)
	f.close()	

    def plot_gauge(self, tg=None, obs=None, weight=None, filename=None, time_param=None, hdffile=False):
	import math as m
	glist=self.glist
	if time_param is None:
	    time_window=self.time_param
	else:
	    time_window=time_param
	n=len(glist)
	gauge_name=self.gauge_name
	
	if hdffile:
	    f=h5py.File(filename+'.h5','w')	
	# number of rows and columns for subplot

	if n<6:
	    nc=n/2
	    fig=plt.figure(figsize=(8, 6))
	else:
	    nc=2
	    fig=plt.figure(figsize=(8, 6))
	nr=int(ceil(n/float(nc)))
	plt.clf()

	fig.subplots_adjust(hspace=.35, wspace=.15)
	Nobsp=0
	sval=[]
	for i in range(n):
	    tw=time_window[glist[i]]
	    TD=tw[1]#min(tw[1], 6000)
	    time=np.linspace(tw[0],TD,tw[2])/60.0
	    ax=fig.add_subplot(nr,nc,i+1)

	    gl=gauge_name[glist[i]]
	    #print gl, glist[i]
	    obs_h=obs[Nobsp:Nobsp+tw[2]]*100
	    print 'max wave height', obs_h.max()
	    tg_h=tg[Nobsp:Nobsp+tw[2]]*100
	    
	    # to create hdf5 file
	    if hdffile:
		f.create_dataset("computed/%s/data"%(glist[i]),data=tg_h)
		f.create_dataset("time/%s/data"%(glist[i]),data=time)
		f.create_dataset("observed/%s/data"%(glist[i]),data=obs_h)
		
	    indmax=np.argmin(abs(obs_h-obs_h.max()))
	    diff_err=abs(obs_h[indmax]-tg_h.max())#[indmax])
	    score=(1-diff_err/abs(obs_h[indmax]))*100	 
	    sval.append(score) 
	    if len(str(glist[i])) <5:
		gls=str(glist[i])
		if gls[0]=='7': gls=gls[1:] # To put the origina code. 
			#In the inversion, we use some stations (very close the harbour) 
			#which are shifted towards the ocean to obtain accurate tsunami simulation result 
		plt.text(time[13],max(obs_h.max(), tg_h.max()), '%s (%s)'%(gl,gls))
	    else:
		plt.text(time[13],max(obs_h.max(), tg_h.max()), '%s'%(gl))
	    plt.plot(time, tg_h,'-r',label='Computed',linewidth=2)
	    if obs.any()!=None:
		ax.hold(True)
		ax.plot(time,obs_h, '--k',label='Observed')
		ax.hold(False)
	    Nobsp=Nobsp+tw[2]
	    simpleaxis(ax)
	    xtk=np.linspace(tw[0],TD,5)/60.0
	    plt.xticks(xtk.round())
	    if abs(obs_h.max()-obs_h.min())>=.6:
		dwh=.2
	    elif abs(obs_h.max()-obs_h.min())<.6 and abs(obs_h.max()-obs_h.min())>=.1:
		dwh=.05
	    else:
		dwh=.01		
	    #ytk=np.arange(round(tg_h.min(),2),round(tg_h.max(),2),dwh)
	    tgm=m.ceil(tg_h.max())
	    ytk=np.linspace(-tgm,tgm,3)
	    #ytk=[round(y,1) for y in ytk]
	    plt.yticks(ytk)
	    if gl=='sanf':
		plt.legend(loc=3, fontsize=9)		

	plt.legend(loc=0, bbox_to_anchor=(1.2, .6), borderaxespad=0., fontsize=10)
	# Set common labels
	fig.text(0.5, 0.05, 'Time (min)', ha='center', va='center', fontsize=13)
	fig.text(0.075, 0.5, 'Water height (cm)', ha='center', va='center', rotation='vertical', fontsize=13)	
	#fig.suptitle(title)	
	plt.savefig(filename+'.png',bbox_inches='tight', pad_inches=0.1)
	print 'score value is: %5.1f'%np.mean(sval)
	
	if hdffile:
	    f.close()
	    
    def drawmap(self):
	
	lllon,urlon,lllat,urlat=self.source_dim
	m = Basemap(lllon,lllat,urlon,urlat,projection='merc',resolution='l')
	m.drawmapboundary(fill_color='aqua') 
	# fill continents, set lake color same as ocean color. 
	m.drawcoastlines()
	m.fillcontinents(color='coral',lake_color='aqua')
	if urlon - lllon < 2.5:
	  dlon = 0.5
	elif urlon - lllon < 5.:
	  dlon =1.0
	else:
	  dlon = 2.
	m.drawparallels(np.arange(-90.,90.,1.0),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,dlon),labels=[0,0,0,1])
	return m
	    
    def write_netcdf(self,filename,rlon,rlat,z,title=None):
	#
	from Scientific.IO.NetCDF import NetCDFFile
	#nc = Dataset(filename,'w',format='NETCDF3_CLASSIC')
	print 'netcdf filename: ', filename
	print rlon[0], rlon[-1], rlat[0], rlat[-1], z.min(), z.max()
	print len(rlon), len(rlat), z.shape
	#nc = netcdf.netcdf_file(filename,'w')
	nc = NetCDFFile(filename,'w')
	if title is None:
	    title=''
	nc.title = title
	nc.source = ''
	nc.createDimension('side',2)
	nc.createDimension('xysize',len(rlon)*len(rlat))
	y_range = nc.createVariable('y_range','d', ('side',))
	y_range.units = 'y'
	y_range[:] = [rlat[0],rlat[-1]]
	x_range = nc.createVariable('x_range','d', ('side',))
	x_range.units = 'x'
	x_range[:] = [rlon[0],rlon[-1]]
	z_range = nc.createVariable('z_range','d', ('side',))
	z_range.units = 'z'
	z_range[:] = [z.min(),z.max()]
	spacing = nc.createVariable('spacing','d',('side',))
	spacing[:] = [rlon[1]-rlon[0],rlat[1]-rlat[0]]
	dimension = nc.createVariable('dimension','i',('side',))
	dimension[:] = [len(rlon),len(rlat)]
	grid_data = nc.createVariable('z','f', ('xysize',))
	grid_data.scale_factor = np.array([1.])
	grid_data.add_offset = np.array([0.])
	grid_data.node_offset = np.array([0])
	q = np.flipud(z)
	q = q.flatten()
	grid_data[:] = q.astype('float32')
	nc.close()
	
    def plot_slip(self, slip, filename=None,slipmax=None):
	    
	sbfs=self.sbfs
	sbf_ind=self.sbf_ind
	num_sbr=self.num_sbr
	plt.clf()
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_axes([0.1, 0.15, .8, 0.6])

	m=self.drawmap()

	cmap=plt.cm.seismic#jet 

	#if slipmax==None:
	slipmax=slip.max()
	slipmin=slip.min()
	slipmag=slipmax-slipmin

	if slipmag==0:
	    slipmag=1.0
	plot_ind=True
	for k in range(num_sbr):
	    sbf=sbfs[k][1]
	    sfply = []
	    for xyz in sbf:
		sfply.append(m(xyz[0],xyz[1]))

	    cmval=(slip[k]+slipmax)/float(2*slipmax)
	    #print k, slip[k], cmval
	    color = cmap(cmval)
	    poly =  Polygon(sfply,facecolor=rgb2hex(color),edgecolor='blue')
	    ax.add_patch(poly)

	    #plt.text(ln,lt,'%s'%int(slip[k]))
	    if plot_ind:
		if k%13==0:
		    cntx,cnty=sbfs[k][0]
		    (ln,lt)=m(cntx,cnty)	
		    plt.text(ln,lt,'%d'%(sbf_ind[k]))

	# draw trench line
	x = []
	y = []
	for line in open('trench.xy','r'):
	    if line.startswith('>'):
		if len(x) > 0:
		    (x0,y0) = m(x,y)
		    plt.plot(x0,y0,'-r',linewidth=2)
		x = []
		y = []
		continue
	    x.append(float(line.split()[0]))
	    y.append(float(line.split()[1]))
	    
	hlon, hlat=self.hypo    
	yc,xc=m(hlon, hlat)     
	plt.plot(yc, xc,'r*',markersize=18)
	norm1 = MidpointNormalize(vmin=-slipmax, vmax=slipmax, midpoint=0)
	#norm1 = mpl.colors.Normalize(vmin=slipmin, vmax=slipmax)
	ax2 = fig.add_axes([0.72, .15, 0.045, 0.6])
	cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
			   norm=norm1,
			   orientation='vertical')
	cb1.set_label('amplitude(m)')		   
	extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())  		       
	plt.savefig(filename,bbox_inches='tight',pad_inches = 0.02)   
	 
    def create_full_source(self, slip, filename=None):
	from math import cos, pi, radians
	sbfs=self.sbfs
	lllon, urlon,lllat, urlat=self.source_dim
	res=60
	ds=res/3600.0
	nx=int((urlon-lllon)/ds +1)
	ny=int((urlat-lllat)/ds +1)
	print 'dimension: (%s,%s)'%(nx,ny)
	rlon=np.linspace(lllon, urlon,nx) #at every 30 sec
	rlat=np.linspace(lllat, urlat,ny)  
	x,y = np.meshgrid(rlon,rlat)
	zf=np.zeros((len(rlat),len(rlon)))
	for k in range(len(slip)):
	    cntr=sbfs[k][0]
	    
	    #print cntr
	    amp=slip[k]
	    z = self.rbf_cosine(cntr,x,y,tpr=15) 
	    zf=zf+amp*z
	if filename is not None:
	    #self.write_dz(filename+'.tt1',x,y,zf)
	    self.write_netcdf(filename+'.grd',rlon,rlat,zf, title='source derived by GFRI')

    def rbf_cosine(self,cntr,x,y,tpr=None):
	'''
	Calculate cosine basis function
	'''
	length=self.sbf_dx
	width=self.sbf_dy
	strike = self.strike
	cntrx,cntry=cntr
	#print length, cntrx,cntry
	g = Geod(ellps='WGS84')
	cbf = np.zeros(x.shape)
	lt2k = 6371.*mt.pi/180.
	ln2k = lt2k*mt.cos(mt.radians(cntry))
	#print ln2k
	rmax = mt.hypot(length,width)
	tmp = np.where(np.hypot((x-cntrx)*ln2k,(y-cntry)*lt2k) < rmax)
	cst = 1#mt.cos(mt.radians(sbf.strike))
	sst = 1#mt.sin(mt.radians(sbf.strike))
	# Take average of lat and lon grid interval
	dgrd = 0.5*((x[0,1]-x[0,0])*ln2k + (y[1,0]-y[0,0])*lt2k) 
	if tpr is None:
	    tpr = 7.5
	for iy,ix in zip(tmp[0],tmp[1]):
	    if True:
		az,baz,dist = g.inv(cntrx,cntry,x[iy,ix],y[iy,ix])
		xp = abs(0.001*dist*mt.cos(mt.radians(az-strike)))
		yp = abs(0.001*dist*mt.sin(mt.radians(az-strike)))
	    else:
		xp = abs( (x[iy,ix]-cntrx)*ln2k*sst + (y[iy,ix]-cntry)*lt2k*cst)
		yp = abs(-(x[iy,ix]-cntrx)*ln2k*cst + (y[iy,ix]-cntry)*lt2k*sst)
	    xbf = 1.
	    if   xp > 0.5*length + tpr:
		xbf = 0.
	    elif xp > 0.5*length - tpr:
		xbf = xbf*0.5*(1.+mt.cos(mt.pi*((xp-0.5*length+tpr)/(2.*tpr))))
	    if   yp > 0.5*width  + tpr:
		xbf = 0.
	    elif yp > 0.5*width  - tpr:
		xbf = xbf*0.5*(1.+mt.cos(mt.pi*((yp-0.5*width+tpr)/(2.*tpr))))
	    #if xbf != 0.:
	    #    print '##%7.3f %7.3f %7.2f %7.2f %g' % (x[iy,ix],y[iy,ix],xp,yp,xbf)
	    cbf[iy,ix] = xbf
	return cbf	
