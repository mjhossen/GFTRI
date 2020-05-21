#!/usr/bin/env python

import os
import sys
from mpl_toolkits.basemap import Basemap,interp
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches 
import numpy as np
from scipy.interpolate import interp2d
#from anuga.shallow_water.shallow_water_domain import Domain
#from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
#from anuga.abstract_2d_finite_volumes.quantity import Quantity
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
#sys.path.append("/usr/lib64/python2.6/dist-packages/")
from Scientific.IO.NetCDF import NetCDFFile
from scipy.io import netcdf
from pyproj import Proj
from matplotlib.colors import Normalize
import matplotlib as mpl

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
	     
def read_dem(filename):
    ncf = netcdf.netcdf_file(filename,'r')
    #ncf = NetCDFFile(filename,'r')
    # Check if not ANUGA dem
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
        lnll = xrng[0]
        lnur = xrng[1]
        ltll = yrng[0]
        ltur = yrng[1]
	nx,ny = ncf.variables['dimension']
	print nx,ny
        zdat = ncf.variables['z'][:].reshape(ny,nx)
        zdat = np.flipud(zdat)
        x = np.linspace(lnll,lnur,nx)
        y = np.linspace(ltll,ltur,ny)        
    elif 'xx_range' in ncf.variables:
        xrng = ncf.variables['x_range']
        yrng = ncf.variables['y_range']
        zrng = ncf.variables['z_range']
        dxdy = ncf.variables['spacing']
        lnll = xrng[0]
        lnur = xrng[1]
        ltll = yrng[0]
        ltur = yrng[1]
        # Pixel registration
        if ncf.variables['z'].node_offset[0] == 1:
            lnll += 0.5*dxdy[0]
            lnur -= 0.5*dxdy[0]
            ltll += 0.5*dxdy[1]
            ltur -= 0.5*dxdy[1]
        nx = 1+int(0.1+(lnur-lnll)/dxdy[0])
        ny = 1+int(0.1+(ltur-ltll)/dxdy[1])
        zdat = ncf.variables['z'][:].reshape(ny,nx)
        zdat = np.flipud(zdat)
        x = np.linspace(lnll,lnur,nx)
        y = np.linspace(ltll,ltur,ny)
    elif ncf.Conventions == 'CF-1.7':
        if 'x' in ncf.variables:
            x = ncf.variables['x'][:]
            y = ncf.variables['y'][:]
        elif 'longitude' in ncf.variables:
            x = ncf.variables['longitude'][:]
            y = ncf.variables['latitude'][:]
        else:
            print 'No x or lon found in COARDS netcdf file %s' %  filename
            reurn(-1)
        zdat=ncf.variables['z'][:]	
    elif ncf.Conventions == 'COARDS, CF-1.5':
        if 'x' in ncf.variables:
            x = ncf.variables['x'][:]
            y = ncf.variables['y'][:]
        elif 'longitude' in ncf.variables:
            x = ncf.variables['longitude'][:]
            y = ncf.variables['latitude'][:]
        else:
            print 'No x or lon found in COARDS netcdf file %s' %  filename
            reurn(-1)
        zdat=ncf.variables['z'][:]
	
    elif ncf.Conventions == 'COARDS/CF-1.0':
        if 'x' in ncf.variables:
            x = ncf.variables['x'][:]
            y = ncf.variables['y'][:]
        elif 'longitude' in ncf.variables:
            x = ncf.variables['longitude'][:]
            y = ncf.variables['latitude'][:]
        else:
            print 'No x or lon found in COARDS netcdf file %s' %  filename
            reurn(-1)
        zdat=ncf.variables['z'][:]	
    else:
        print 'File %s is not  ANUGA or GMT' % filename
        return(-1)
    ncf.close()
    return (x,y,zdat)

def drawmap(lllon,lllat,urlon,urlat):
   m = Basemap(lllon,lllat,urlon,urlat,projection='merc',resolution='l')
   m.drawmapboundary(fill_color='aqua') 
   # fill continents, set lake color same as ocean color. 
   m.drawcoastlines()
   #m.fillcontinents(color='coral',lake_color='aqua')
   if urlon - lllon < 2.5:
      dlon = 0.5
   elif urlon - lllon < 5.:
      dlon =1.0
   else:
      dlon = 1.
   m.drawparallels(np.arange(-90.,90.,1.0),labels=[1,0,0,0])
   m.drawmeridians(np.arange(-180.,180.,2.0),labels=[0,0,0,1])
   return m
   
   

def plot_data(fname, x, y, z):
    from matplotlib.mlab import griddata as mgrid
    from scipy.interpolate import griddata as sgrid  
    
    lnll=min(x)
    ltll=min(y)
    lnur=max(x)
    ltur=max(y)
    zone=54
    xi = np.linspace(lnll,lnur,100)
    yi = np.linspace(ltll,ltur,100)
    xig, yig=np.meshgrid(xi,yi)
    xg, yg=np.meshgrid(x,y)   
    points=zip(xg.flatten(),yg.flatten())
    #zi = mgrid(xg.flatten(),yg.flatten(),z.flatten(),xi,yi,interp='linear')
    zi = sgrid(points,z.flatten(),(xig,yig),method='nearest')
    plt.clf()
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_axes([0.1, 0.15, .8, 0.6])
    m = drawmap(lnll,ltll,lnur,ltur)
    xx, yy=m(xi,yi)
    xg, yg=np.meshgrid(xx,yy)
    #zi=zi*100 # convert to cm
    hm=zi.max()

    cmap1 = mpl.cm.seismic
    norm1 = MidpointNormalize(vmin=zi.min(), vmax=zi.max(), midpoint=0)
    #norm1 = MidpointNormalize(vmin=-hm, vmax=hm, midpoint=0)
    #clevs=np.arange(-hm,hm,2)
    clevs=np.linspace(-hm,hm,10)
    m.contourf(xg,yg,zi, len(clevs), cmap=cmap1, norm=norm1)
    clevs=[-.1, -.05, -.025, 0, .05, .1,.2,.3]
    CS=m.contour(xg,yg,zi,clevs,linewidths=0.5,colors='w')
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
	xx=float(line.split()[0])
	yy=float(line.split()[1])
	if xx <0:
	    xx=360+xx
	x.append(xx)
	y.append(yy)
    hypo_lon=210.9
    hypo_lat=56	
    yc,xc=m(hypo_lon,hypo_lat)  
    plt.plot(yc, xc,'g*',markersize=18)    
    ax2 = fig.add_axes([0.76, .15, 0.03, 0.6])
    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap1,
		       norm=norm1,
		       orientation='vertical')
    cb1.set_label('Water height (m)', fontsize=13)
    CS.levels = [nf(val) for val in CS.levels ]
    # Label levels with specially formatted floats
    if plt.rcParams["text.usetex"]:
	 fmt = r'%r'
    else:
	 fmt = '%r'
    plt.clabel(CS, CS.levels, fontsize=9,fmt=fmt, inline=True, colors='k')
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())  		       
    plt.savefig(fname.split('.')[0]+'.png',bbox_inches='tight', pad_inches = 0.02)
    
# main part
if __name__=="__main__":
    #time=sys.argv[1]
    #grdfile='disp_outer_padding_t%s.grd'%time
    grdfile=sys.argv[1]
    (x,y,zdat) = read_dem(grdfile)
    plot_data(grdfile, x,y,zdat)

