import numpy as np
import matplotlib.pyplot as plt
import my_toolbox as my
from netCDF4 import Dataset
import sys

savefig = True
path_fig = '../pics/'
nnf=0

run = 'run_Ekman_const'; name = 'const'
#run = 'run_Ekman_tke'; name = 'tke'
#run = 'run_Ekman_tke_const_f'; name = 'const_f'
path_data = '/Users/nbruggemann/work/proj_vmix_mw/onedmix/'+run+'/'
fnamenc = 'OneDmix.nc'

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
try:
  tke_mean = tke[it,:].mean(axis=0)
  tke_Tbck_mean = tke_Tbck[it,:].mean(axis=0)
  tke_Tbpr_mean = tke_Tbpr[it,:].mean(axis=0)
  tke_Tdif_mean = tke_Tdif[it,:].mean(axis=0)
  tke_Tdis_mean = tke_Tdis[it,:].mean(axis=0)
  tke_Tiwf_mean = tke_Tiwf[it,:].mean(axis=0)
  tke_Tspr_mean = tke_Tspr[it,:].mean(axis=0)
  tke_Ttot_mean = tke_Ttot[it,:].mean(axis=0)
  tke_Twin_mean = tke_Twin[it,:].mean(axis=0)
  plot_tke_budget = True
except:
  tke_mean = np.zeros((nz+1))
  tke_Tbck_mean = np.zeros((nz+1))
  tke_Tbpr_mean = np.zeros((nz+1))
  tke_Tdif_mean = np.zeros((nz+1))
  tke_Tdis_mean = np.zeros((nz+1))
  tke_Tiwf_mean = np.zeros((nz+1))
  tke_Tspr_mean = np.zeros((nz+1))
  tke_Ttot_mean = np.zeros((nz+1))
  tke_Twin_mean = np.zeros((nz+1))
  plot_tke_budget = False

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

# ================================================================================ 
# Here starts plotting
# ================================================================================ 
plt.close("all")

def get_color(nclev, cmap='viridis'):
  cmap_mpl = plt.cm.get_cmap(cmap)
  cols = [0]*nclev
  for nn in range(nclev):
    cols[nn] = cmap_mpl(float(nn)/nclev)
  return cols

import matplotlib as mpl
def get_cols(vals, clim='auto', cmap='viridis'):
  vals = np.array(vals)
  if isinstance(clim, str) and clim=='auto':
    clim = [vals.min(), vals.max()]
  cmap_mpl = mpl.cm.get_cmap(cmap)
  norm = mpl.colors.Normalize(vmin=clim[0], vmax=clim[1])
  cols = cmap_mpl( norm(vals) )
  return cols, norm

def make_cbar(cax, norm, cmap='viridis'):
  cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
  return

tp = time/86400. 
ts = 'time [days]'

cols = get_color(len(time), cmap='viridis')

# --- time series of state variables
hca, hcb = my.arrange_axes(2,1, plot_cb=[0,1], sasp=0.66, fig_size_fac=2., \
                            sharex=True, sharey=True, xlabel=ts, ylabel="z [m]")
ii=-1

pvars = ['uvel', 'vvel']
tits = ['$u$ [m/s]', '$v$ [m/s]']
clims = ['sym']*len(pvars)
for nn, var in enumerate(pvars):
  exec( 'data = 1.0*(%s)' % var)
  ii+=1; ax=hca[ii]; cax=hcb[ii]
  #my.shade2y(tp, zt, data.transpose(), ys=[-3000., -200., 0.], ax=ax, cax=cax, clim=clims[nn], ylabel="z [m]", xlabel=ts)
  my.shade(tp, zt, data.transpose(), ax=ax, cax=cax, clim=clims[nn])
  ax.set_title(tits[nn])

for ax in hca:
  ax.axvline(tp[it[0]], color='0.7')
  ax.axvline(tp[it[-1]], color='0.7')

nnf+=1
if savefig:
  print 'save figure: %s_%s_%02d.pdf' % (__file__.split('/')[-1][:-3], name, nnf)
  plt.savefig("%s_%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], name, nnf))

# --- time series of state variables
hca, hcb = my.arrange_axes(3,1, plot_cb=True, sasp=1.0, fig_size_fac=2., \
                            sharex=False, sharey=False, xlabel="", ylabel="", 
                            projection=['rectilinear', '3d', 'rectilinear'], axlab_kw=None,
                          )
ii=-1

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(uvel[:,0], vvel[:,0], zorder=1)
hp = ax.scatter(uvel[:,0], vvel[:,0], c=tp, zorder=2)
plt.colorbar(mappable=hp, cax=cax)
ht = cax.set_title('days', fontsize=10)
ax.set_xlabel('$u$ [m/s]')
ax.set_ylabel('$v$ [m/s]')
ax.set_title('surface velocity')

#ii+=1; ax=hca[ii]; cax=hcb[ii]
#iz = np.argmin(np.abs(zt+3000)) 
#ax.plot(uvel_mean[:iz], vvel_mean[:iz])
#hp = ax.scatter(uvel_mean[:iz], vvel_mean[:iz], c=zt[:iz])
#plt.colorbar(mappable=hp, cax=cax)
#ax.set_xlabel('uvel [m/s]')
#ax.set_ylabel('vvel [m/s]')

cols, norm = get_cols(tp[it], clim='auto', cmap='viridis')
ii+=1; ax=hca[ii]; cax=hcb[ii]
for nn, itl in enumerate(it):
  ax.plot(uvel[itl,:], vvel[itl,:], zt, color=cols[nn,:])
ax.plot(uvel[it,:].mean(axis=0), vvel[it,:].mean(axis=0), zt, color='k')
ax.scatter(uvel[it,kek].mean(axis=0), vvel[it,kek].mean(axis=0), zt[kek])
make_cbar(cax, norm, cmap='viridis')
ht = cax.set_title('days', fontsize=10)
ax.set_title('velocity profile')
ax.set_xlabel('$u$ [m/s]')
ax.set_ylabel('$v$ [m/s]')
ax.set_zlabel('z [m]')
ax.view_init(azim=-9, elev=18)

cols, norm = get_cols(zt, clim='auto', cmap='viridis')
ii+=1; ax=hca[ii]; cax=hcb[ii]
#for k, z in enumerate(zt[:iz]):
for k, z in enumerate(zt):
  #ax.scatter([0,uvel[it,k]], [0,vvel[it,k]], c=[zt[k]]*2, linestyle='-', vmin=zt[iz], vmax=zt[0])
  ax.plot([0,uvel_mean[k]], [0,vvel_mean[k]], color=cols[k,:])
#ax.plot(uvel_mean, vvel_mean)
uscale = 0.5 * np.sqrt(uvel_mean**2+vvel_mean**2).max()
ekscale = np.sqrt(Uek**2+Vek**2)
tauscale = np.sqrt(taux_act_mean**2+tauy_act_mean**2)
ax.plot([0,Uek*uscale/ekscale], [0,Vek*uscale/ekscale], color='b', label='Ekman transport')
ax.plot([0,Uek_ekd*uscale/ekscale], [0,Vek_ekd*uscale/ekscale], color='c')
ax.plot([0,taux_act_mean*uscale/tauscale], [0,tauy_act_mean*uscale/tauscale], color='r', label='wind stress')
#ax.set_xlim(-1,1)
#ax.set_ylim(-1,1)
ax.set_xlabel('$u$ [m/s]')
ax.set_ylabel('$v$ [m/s]')
make_cbar(cax, norm, cmap='viridis')
ht = cax.set_title('z [m]', fontsize=10)
ax.set_title('velocity direction')
ax.legend(loc='lower left')

nnf+=1
if savefig:
  print 'save figure: %s_%s_%02d.pdf' % (__file__.split('/')[-1][:-3], name, nnf)
  plt.savefig("%s_%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], name, nnf))

# --- mean profiles
hca, hcb = my.arrange_axes(5,1, plot_cb=False, sasp=1.5, fig_size_fac=2., \
                            sharex=True, sharey=True, xlabel="", ylabel="z [m]")
ii=-1

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(uvel_mean, zt, color='r', label='$u$')
ax.plot(vvel_mean, zt, color='b', label='$v$')
ax.plot(vabs, zt, color='k', label=r'$|\vec{u}|$')
ax.legend(loc='best')
ax.set_title('velocity [m/s]')

ii+=1; ax=hca[ii]; cax=hcb[ii]
Av_vmean = (Av_mean*dzt).sum() / dzt.sum()
Av_vmean_dek100 = (Av_mean*dzt)[:kek100].sum() / dzt[:kek100].sum()
ax.plot(Av_mean, zu, color='r', label='$A_v$')
ax.axvline(Av_vmean, color='k')
ax.axvline(Av_vmean_dek100, color='0.5')
ax.set_title('viscosity $A_v$ [m$^2$/s]')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(N2_mean/1e-6, zu, color='r', label='$<N^2>$')
ax.plot(N2[0,:]/1e-6, zu, color='b', label='$N^2(t=0)$')
ax.set_title('$N^2$ [10$^{-6}$/s]')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(tke_mean/1e-4, zu, color='r', label='TKE')
ax.set_title('TKE [cm$^2$/s$^2$]')

ii+=1; ax=hca[ii]; cax=hcb[ii]
ax.plot(temp_mean, zt, color='r')
ax.plot(temp[0,:], zt, color='b')
ax.set_title('temp')

for ax in hca:
  ax.axhline(-dek, color='0.7')
  ax.axhline(-dek100, color='0.7')

nnf+=1
if savefig:
  print 'save figure: %s_%s_%02d.pdf' % (__file__.split('/')[-1][:-3], name, nnf)
  plt.savefig("%s_%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], name, nnf))

# --- tke budget
if plot_tke_budget:
  hca, hcb = my.arrange_axes(3,1, plot_cb=False, sasp=1.0, fig_size_fac=2., \
                              sharex=False, sharey=False, xlabel="", ylabel="",
                              oyt = 0.1, 
                              )
  ii=-1

  pvars = [
    'tke_Tbpr', 'tke_Tspr', 'tke_Tdif', 'tke_Tdis', 
    'tke_Twin', 'tke_Tiwf', 'tke_Tbck', 'tke_Ttot', 
    #'tke_Tres',
          ]

  ii+=1; ax=hca[ii]; cax=hcb[ii]
  for var in pvars:
    exec( 'data = 1.*(%s)' % var)
    data = (data*dzt).sum(axis=1) / dzt.sum()
    ax.plot(tp, data/1e-6, label=var)
  ax.set_title("vertically integrated\nTKE budget / 1e-6 [m$^2$/s$^3$]")
  ax.locator_params(nbins=4)
  ax.set_xlabel(ts)
  ax.set_ylim(-0.01,0.01)
  ax.axvline(tp[it[0]], color='0.7')
  ax.axvline(tp[it[-1]], color='0.7')

  ii+=1; ax=hca[ii]; cax=hcb[ii]
  for var in pvars:
    exec( 'data = 1.*(%s_mean)' % var)
    ax.plot(data/1e-6, zu, label=var)
  ax.set_title("TKE budget / 1e-6 [m$^2$/s$^3$]")
  ax.locator_params(nbins=4)
  ax.set_ylabel('z [m]')
  ax.set_ylim(zu[-1],zu[0])
  ax.legend(loc='best')

  ii+=1; ax=hca[ii]; cax=hcb[ii]
  for var in pvars:
    exec( 'data = 1.*(%s_mean)' % var)
    ax.plot(data/1e-6, zu, label=var)
  ax.set_title("TKE budget / 1e-6 [m$^2$/s$^3$]")
  ax.locator_params(nbins=4)
  ax.set_ylabel('z [m]')
  ax.set_ylim(0.5*zu[-1],zu[0])
  #ax.set_xlim(np.array([-1,1])*np.abs(tke_Tdif/1e-6).max())
  ax.set_xlim(-0.1,0.1)

nnf+=1
if savefig:
  print 'save figure: %s_%s_%02d.pdf' % (__file__.split('/')[-1][:-3], name, nnf)
  plt.savefig("%s_%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], name, nnf))

plt.show()
sys.exit()

# --- time series of mixing parameters
hca, hcb = my.arrange_axes(3,1, plot_cb=True, sasp=0.66, fig_size_fac=2., \
                            sharex=ts, sharey='z [m]')
ii=-1

ii+=1; ax=hca[ii]; cax=hcb[ii]
my.shade(tp, zu, uz.transpose(), ax=ax, cax=cax, clim='auto')
ax.set_title('uz')

ii+=1; ax=hca[ii]; cax=hcb[ii]
my.shade(tp, zu, vz.transpose(), ax=ax, cax=cax, clim='auto')
ax.set_title('vz')

ii+=1; ax=hca[ii]; cax=hcb[ii]
my.shade(tp, zu, S2_mean.transpose(), ax=ax, cax=cax, clim='auto')
ax.set_title('S2_mean')

for ax in hca:
  ax.axvline(tp[it[0]], color='0.7')
  ax.axvline(tp[it[-1]], color='0.7')

# --- time series of mixing parameters
hca, hcb = my.arrange_axes(3,2, plot_cb=True, sasp=0.66, fig_size_fac=2., \
                            sharex=True, sharey=True, xlabel=ts, ylabel='z [m]')
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
  #my.shade2y(tp, zu, data.transpose(), ys=[-3000., -200., 0.], ax=ax, cax=cax, clim=[-6,0], ylabel="z [m]", xlabel=ts)
  my.shade(tp, zu, data.transpose(), ax=ax, cax=cax, clim=[-6,0])
  ax.set_title('log$_{10}$(%s)'%var)

ii+=1; ax=hca[ii]; cax=hcb[ii]
my.shade(tp, zu, (temp-temp_mean).transpose(), ax=ax, cax=cax, clim='auto')
ax.set_title('temp - temp_mean')

for ax in hca:
  ax.axvline(tp[it[0]], color='0.7')
  ax.axvline(tp[it[-1]], color='0.7')

# ---
hca, hcb = my.arrange_axes(3,2, plot_cb=False, sasp=0.66, fig_size_fac=2., \
                            sharex=True, sharey=False, xlabel=ts, ylabel="")
ii=-1

ke = (0.5*(uvel**2+vvel**2)*dzw).sum(axis=1)
pvars = [u'taux_act', u'tauy_act', u'emp_act', u'q0_act', 'ke', 'uvel[:,-1]']
for var in pvars:
  exec( 'data = (%s)' % var)
  ii+=1; ax=hca[ii]; cax=hcb[ii]
  ax.plot(tp, data, marker='x')
  ax.set_title(var)
  #ax.plot(u0, zt, color='orange')
  ax.locator_params(nbins=4)

plt.show()
