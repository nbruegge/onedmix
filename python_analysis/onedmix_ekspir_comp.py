import numpy as np
import matplotlib.pyplot as plt
import my_toolbox as my
from netCDF4 import Dataset
import sys

savefig = True
path_fig = '../pics/'
nnf=0

class Run(object):
  def __init__(self):
    return

Rs = []

R=Run(); R.run = 'run_Ekman_const'; R.name = 'const $A_v$'
Rs.append(R)
R=Run(); R.run = 'run_Ekman_tke'; R.name = 'tke'
Rs.append(R)


for nn, R in enumerate(Rs):
  path_data = '/Users/nbruggemann/work/proj_vmix_mw/onedmix/'+R.run+'/'

  # --- read namelist
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
  
  # --- load data
  f = Dataset(path_data+'OneDmix.nc', 'r')
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
  dens += -1000.
  
  # --- time averages
  #it = 20
  #it = [it,it+1]
  #it = np.arange(5,20)
  it = np.arange(100,100+4*10)
  uvel_mean = uvel[it,:].mean(axis=0) 
  vvel_mean = vvel[it,:].mean(axis=0) 
  taux_act_mean = taux_act[it].mean(axis=0) 
  tauy_act_mean = tauy_act[it].mean(axis=0) 
  Av_mean = Av[it,:].mean(axis=0)
  N2_mean = N2[it,:].mean(axis=0)
  temp_mean = temp[it,:].mean(axis=0)
  
  # --- Ekman depth
  vabs = np.sqrt(uvel_mean**2+vvel_mean**2)
  kek = np.argmin(np.abs(vabs - vabs[0]/np.exp(1)))#, axis=1)
  dek = -zt[kek]
  kek100 = np.argmin(np.abs(vabs - vabs[0]/np.exp(1)/100.))#, axis=1)
  dek100 = -zt[kek100]
  
  # --- Ekman transport
  # integrating up to Ekman depth
  Uek_ekd = (uvel_mean*dzw)[:kek].sum()
  Vek_ekd = (vvel_mean*dzw)[:kek].sum()
  # integrating over total water column
  Uek = (uvel_mean*dzw).sum()#axis=1)
  Vek = (vvel_mean*dzw).sum()#axis=1)
  
  # --- theory
  # theoretical values
  Uektheo =  tauy_act_mean/fCor
  Vektheo = -taux_act_mean/fCor
  dektheo = np.sqrt(2*Avb/np.abs(fCor))
  
  # --- calculate shear components
  nt = uvel.shape[0]
  uz = np.zeros((nt, nz+1))
  uz[:,1:-1] = (uvel[:,:-1]-uvel[:,1:])/(zt[:-1]-zt[1:])
  vz = np.zeros((nt, nz+1))
  vz[:,1:-1] = (vvel[:,:-1]-vvel[:,1:])/(zt[:-1]-zt[1:])
  
  S2_mean = uz**2+vz**2
  
  Rvars = ['uvel_mean', 'vvel_mean', 'taux_act_mean', 'tauy_act_mean']
  for var in Rvars:
    exec('R.%s = 1.*%s' % (var, var))

# ================================================================================ 
# Here starts plotting
# ================================================================================ 
plt.close('all')
# --- time series of state variables
projection = ['rectilinear']*3; projection[0]='3d'
hca, hcb = my.arrange_axes(3,1, plot_cb=False, sasp=1.0, fig_size_fac=2., \
                            sharex=False, sharey=False, xlabel="", ylabel="", 
                            projection=projection, axlab_kw=None,
                          )
ii=-1

cols = ['C0', 'C1']
ii+=1; ax=hca[ii]; cax=hcb[ii]
for nn, R in enumerate(Rs):
  ax.plot(R.uvel_mean, R.vvel_mean, zt, color=cols[nn], label=R.name)
ax.set_title('velocity profile')
ax.set_xlabel('$u$ [m/s]')
ax.set_ylabel('$v$ [m/s]')
ax.set_zlabel('z [m]')
ax.view_init(azim=-9, elev=18)

ii+=1; ax=hca[ii]; cax=hcb[ii]
for nn, R in enumerate(Rs):
  for k, z in enumerate(zt):
    ax.plot([0,R.uvel_mean[k]], [0,R.vvel_mean[k]], color=cols[nn])
  #ax.plot(uvel_mean, vvel_mean)
  uscale = 0.5 * np.sqrt(R.uvel_mean**2+R.vvel_mean**2).max()
  ekscale = np.sqrt(Uek**2+Vek**2)
  tauscale = np.sqrt(R.taux_act_mean**2+R.tauy_act_mean**2)
  if nn==0:
    ax.plot([0,Uek*uscale/ekscale], [0,Vek*uscale/ekscale], color='b', label='Ekman transport')
    #ax.plot([0,Uek_ekd*uscale/ekscale], [0,Vek_ekd*uscale/ekscale], color='c')
    ax.plot([0,taux_act_mean*uscale/tauscale], [0,tauy_act_mean*uscale/tauscale], color='r', label='wind stress')

ax.set_xlabel('$u$ [m/s]')
ax.set_ylabel('$v$ [m/s]')
ax.set_title('velocity direction')
first_legend = ax.legend(loc='lower left')
ax.add_artist(first_legend)
l1, = ax.plot([],[], color=cols[0], linestyle='-', label=Rs[0].name)
l2, = ax.plot([],[], color=cols[1], linestyle='-', label=Rs[1].name)
ax.legend(handles=[l1,l2], loc='upper right')

ii+=1; ax=hca[ii]; cax=hcb[ii]
for nn, R in enumerate(Rs):
    ax.plot(R.uvel_mean, zt, color=cols[nn])
    ax.plot(R.vvel_mean, zt, color=cols[nn], linestyle='--')

# legend
l1, = ax.plot([],[], color=cols[0], linestyle='-', label=Rs[0].name)
l2, = ax.plot([],[], color=cols[1], linestyle='-', label=Rs[1].name)
l3, = ax.plot([],[], color='k', linestyle='-', label='$u$')
l4, = ax.plot([],[], color='k', linestyle='--', label='$v$')
ax.set_title('velocity [m/s]')
ax.set_ylabel('z [m]')
ax.legend()

nnf+=1
if savefig:
  print 'save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3],  nnf)
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

plt.show()
