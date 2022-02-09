import numpy as np
import matplotlib.pyplot as plt
import my_toolbox as my
from netCDF4 import Dataset
import sys

savefig = False
path_fig = '../pics_nac/'
nnf=0

run = 'run_03'
path_data = '/Users/nbruggemann/work/proj_vmix_mw/onedmix/'+run+'/'
fnamenc = 'OneDmix.nc'

# read namelist
f = open(path_data+'OneDmix.nl')
for line in f.readlines():
  if '=' in line:
    line = line[1:-1]
    if line[-1]==',':
      line = line[:-1]
    try:                                                                             
      exec(line)                                                                     
    except:                                                                          
      vnam = line.split('=')[0]                                                      
      arg = line.split('=')[1]                                                       
      exec('%s = \'%s\'' % (vnam, arg))
f.close()

f = Dataset(path_data+fnamenc, 'r')
zt = f.variables['zt'][:]
zu = f.variables['zu'][:]
dzt = f.variables['dzt'][:]
dzw = f.variables['dzw'][:]
time = f.variables['t'][:]
lvars = f.variables.keys()
lvars = lvars[7:] # exclude grid
for var in lvars:
  if f.variables[var].ndim==4:
    vard = f.variables[var][:,:,0,0]
  if f.variables[var].ndim==3:
    vard = f.variables[var][:,0,0]
  exec('%s = vard'%var)
f.close()

H0 = dzw.sum()

tp = time/86400. 
ts = 'time [days]'

dens += -1000.
#N2 *= 1e6
#S2 *= 1e6
#Av *= 1e6
#kv *= 1e6

#N22 = np.zeros((ntt+1, nz))
#N22[:,1:] = -grav/rho0*(dens[:,:-1]-dens[:,1:])/dzt[np.newaxis,1:]
#N22[:,0] = N22[:,1]
#N22 *= 1e6

## list of variable names 
#fobj = open(path_data+'OneDmix_2D_varlist.txt', 'r')                                   
#varnames2D = fobj.readlines()                                                        
#nvars = len(varnames2D)                                                              
#for nn in range(nvars):                                                              
#  varnames2D[nn] = varnames2D[nn][1:-1]                                              
#  fobj.close() 
#
#fobj = open(path_data+'OneDmix_2Dp1_varlist.txt', 'r')                                   
#varnames2Dp1 = fobj.readlines()                                                        
#nvars = len(varnames2Dp1)                                                              
#for nn in range(nvars):                                                              
#  varnames2Dp1[nn] = varnames2Dp1[nn][1:-1]                                              
#  fobj.close() 
#
#fobj = open(path_data+'OneDmix_1D_varlist.txt', 'r')                                   
#varnames1D = fobj.readlines()                                                        
#nvars = len(varnames1D)                                                              
#for nn in range(nvars):                                                              
#  varnames1D[nn] = varnames1D[nn][1:-1]                                              
#  fobj.close() 

# ================================================================================ 
# Here starts plotting
# ================================================================================ 
plt.close("all")

nclev = len(time)
cmap = plt.cm.get_cmap('jet')
cols = [0]*nclev
for nn in range(nclev):
  cols[nn] = cmap(float(nn)/nclev)

# --- time series of state variables
hca, hcb = my.arrange_axes(3,2, plot_cb=True, sasp=0.66, fig_size_fac=1.5, \
                            sharex=False, sharey=False, xlabel="", ylabel="")
ii=-1

#pvars = ['uvel', 'vvel', 'buoy']
#pvars = varnames2D 
#pvars = ['uvel']
pvars = ['uvel', 'vvel', 'temp', 'salt', 'dens', 'temp-temp[0,:]']
clims = ['auto']*len(pvars)
clims[2] = [15,18]
clims[-1] = 'sym'
for nn, var in enumerate(pvars):
  exec( 'data = 1.0*(%s)' % var)
  ii+=1; ax=hca[ii]; cax=hcb[ii]
  #for l in range(len(time)):
  #  ax.plot(data[l,:], zt, color=cols[l])
  #ax.locator_params(nbins=4)
  my.shade2y(tp, zt, data.transpose(), ys=[-3000., -200., 0.], ax=ax, cax=cax, clim=clims[nn], ylabel="z [m]", xlabel=ts)
  ax.set_title(var)

nnf+=1
if savefig:
  print 'save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf)
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

# --- time series of mixing parameters
hca, hcb = my.arrange_axes(3,2, plot_cb=True, sasp=0.66, fig_size_fac=2., \
                            sharex=False, sharey=False)
ii=-1

#pvars = ['uvel', 'vvel', 'buoy']
#pvars = varnames2D 
#pvars = ['uvel']
#pvars = ['uvel', 'vvel', 'temp', 'salt', 'dens', 'ptra']
pvars = ['Av', 'kv', 'N2', 'S2', 'Ri']
for var in pvars:
  exec( 'data = 1.0*(%s)' % var)
  data[data<=0.0] = np.ma.masked
  data = np.log10(data)
  ii+=1; ax=hca[ii]; cax=hcb[ii]
  #for l in range(len(time)):
  #  ax.plot(data[l,:], zt, color=cols[l])
  #ax.locator_params(nbins=4)
  my.shade2y(tp, zu, data.transpose(), ys=[-3000., -200., 0.], ax=ax, cax=cax, clim=[-6,0], ylabel="z [m]", xlabel=ts)
  ax.set_title('log$_{10}$(%s)'%var)

nnf+=1
if savefig:
  print 'save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf)
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

# --- 2Dp1 variables
hca, hcb = my.arrange_axes(7,2, plot_cb=False, sasp=1.5, fig_size_fac=2., \
                            sharex=False, sharey=True, xlabel="", ylabel="z [m]")
ii=-1

pvars  = ['uvel', 'vvel', 'temp', 'salt', 'dens', 'ptra']
pvars += ['Av', 'kv', 'N2', 'S2', 'Ri', 'tke']
for var in pvars:
  exec( 'data = (%s)' % var)
  ii+=1; ax=hca[ii]; cax=hcb[ii]
  if data.shape[1]==nz:
    zp = 1.*zt
  else:
    zp = 1.*zu
  for l in range(len(time)):
    ax.plot(data[l,:], zp, color=cols[l])
  ax.set_title(var)
  ax.locator_params(nbins=4)

tke_Tres = tke_Ttot - (tke_Tbpr + tke_Tspr + tke_Tdif + tke_Tdis + tke_Twin + tke_Tiwf + tke_Tbck)
#pvars = ['tke_Tdif', 'tke_Tdis', 'tke_Tbpr', 'tke_Tspr', 'tke_Twin', 'tke_Ttot', 'tke_Tres']
pvars = ['tke_Tbpr', 'tke_Tspr', 'tke_Tdif', 'tke_Tdis', 'tke_Twin', 'tke_Tiwf', 'tke_Tbck', 'tke_Ttot', 'tke_Tres']
ii+=1; ax=hca[ii]; cax=hcb[ii]
for var in pvars:
  exec( 'data = (%s)' % var)
  ax.plot(data[-1,:]/1e-6, zu, label=var)
ax.set_title("TKE budget / 1e-6")
ax.locator_params(nbins=4)

ii+=1; ax=hca[ii]; cax=hcb[ii]
for var in pvars:
  exec( 'data = (%s)' % var)
  ax.plot(data[-1,:]/1e-6, zu, label=var)
ax.set_title("TKE budget / 1e-6")
ax.legend()
ax.locator_params(nbins=4)

nnf+=1
if savefig:
  print 'save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf)
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

# ---
hca, hcb = my.arrange_axes(3,2, plot_cb=False, sasp=0.66, fig_size_fac=2., \
                            sharex=True, sharey=False, xlabel=ts, ylabel="")
ii=-1

pvars = [u'taux_act', u'tauy_act', u'emp_act', u'q0_act']
for var in pvars:
  exec( 'data = (%s)' % var)
  ii+=1; ax=hca[ii]; cax=hcb[ii]
  ax.plot(tp, data, marker='x')
  ax.set_title(var)
  #ax.plot(u0, zt, color='orange')
  ax.locator_params(nbins=4)

#pvars = ['tke_Tdif', 'tke_Tdis', 'tke_Tbpr', 'tke_Tspr', 'tke_Twin', 'tke_Ttot', 'tke_Tres']
pvars = ['tke_Tbpr', 'tke_Tspr', 'tke_Tdif', 'tke_Tdis', 'tke_Twin', 'tke_Tiwf', 'tke_Tbck', 'tke_Ttot', 'tke_Tres']
ii+=1; ax=hca[ii]; cax=hcb[ii]
for var in pvars:
  exec( 'data = (%s)' % var)
  ax.plot(tp, data[:,0]/1e-6, label=var)
ax.set_title("TKE budget / 1e-6")
ax.locator_params(nbins=4)

ii+=1; ax=hca[ii]; cax=hcb[ii]
for var in pvars:
  exec( 'data = (%s)' % var)
  #ax.plot(tp, data[:,4]/1e-6, marker='.', label=var)
  ax.plot(tp, (data*dzt).sum(axis=1)/H0/1e-6, marker='x', label=var)
ax.set_title("TKE budget / 1e-6")
ax.locator_params(nbins=4)

nnf+=1
if savefig:
  print 'save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf)
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

# ---
hca, hcb = my.arrange_axes(1,1, plot_cb=False, sasp=0.66, fig_size_fac=2., \
                            sharex=False, sharey=False, xlabel=ts, ylabel="")
ii=-1

ke = (0.5*(uvel**2+vvel**2)*dzw).sum(axis=1) / H0
ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(tp, ke)
ax.set_title('ke [m$^2$/s$^2$]')

nnf+=1
if savefig:
  print 'save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf)
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

plt.show()
