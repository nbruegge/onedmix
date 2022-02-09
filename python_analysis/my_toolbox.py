import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
try:
  import mycmap
except:
  pass
import datetime

"""
functions contained in this toolbox:

# ================================================================================ 
# plotting functions (need to import matplotlib.pyplot and numpy)
# ================================================================================
- shade
- shade2y
- colorbar
- getticks
- plotsettings
- arrange_axes
- axlab
- autoclim
- figinfo
- dcmap
- saveshow
- make_inlay_axes
- sort_contour
- plt_stdcols
- tbox
- nlcmap

# ================================================================================ 
# map plotting functions (need to import mpl_toolkits.basemap, matplotlib.pyplot and numpy)
# ================================================================================ 
- plot_map

# ================================================================================ 
# math functions (need to import numpy)
# ================================================================================ 
- vectogrid
- stats
- nanmean
- bootstrap
- repmat
- roundsign
- indfind
- bilininterp
- sort_contour
- get_island_mask

# ================================================================================ 
# mpi4py object and functions
# ================================================================================ 

"""

# ================================================================================ 
# plotting functions (need to import matplotlib.pyplot and numpy)
# ================================================================================ 
def shade(x, y, datai,
            ax='auto', cax=0,
            cmap='auto',
            rasterized=True,
            clim=[None, None],
            extend='both',
            conts=None,
            nclev='auto',
            cint='auto',
            contcolor='k',
            contthick=0.,
            contfs=None,
            contlw=1.,
            use_pcol=True,
            cbticks='auto',
            adjust_axlims=True,
            bmp=None,
            transform=None,
            logplot=False,
            cbtitle='',
            edgecolor='none',
         ):
  """ Makes a nice pcolor(mesh) plot.

last change:
----------
2016-08-23
  """
  # mask 0 and negative values in case of log plot
  data = 1.*datai
  if logplot and isinstance(data, np.ma.MaskedArray):
    data[data<=0.0] = np.ma.masked
    data = np.ma.log10(data) 
  elif logplot and not isinstance(data, np.ma.MaskedArray):
    data[data<=0.0] = np.nan
    data = np.log10(data) 

  #clims
  if isinstance(clim, str) and clim=='auto':
    clim = [None, None]
  elif isinstance(clim, str) and clim=='sym':
    clim = np.abs(data).max()
  clim=np.array(clim)
  if clim.size==1:
    clim = np.array([-1, 1])*clim
  if clim[0] is None:
    clim[0] = data.min()
  if clim[1] is None:
    clim[1] = data.max()

  if (clim[0]==-clim[1]) and cmap=='auto':
    cmap = 'RdBu_r'
  elif cmap=='auto':
    #cmap = 'viridis'
    cmap = 'RdYlBu_r'

  # calculate contour x/y and contour levels if needed
  if conts is None:
    use_cont = False
  elif isinstance(conts,str) and conts=='auto':
    use_cont = True
    if isinstance(nclev,str) and nclev=='auto':
      conts = np.linspace(clim[0], clim[1], 11)
    else:
      conts = np.linspace(clim[0], clim[1], nclev)
    if not (isinstance(cint,str) and cint=='auto'):
      conts = np.arange(clim[0], clim[1]+cint, cint)
  else:
    use_cont = True
    conts = np.array(conts)

  if contfs is None:
    use_contf=False
  elif isinstance(contfs, str) and contfs=='auto':
    use_contf=True
    use_pcol=False
    if isinstance(nclev,str) and nclev=='auto':
      contfs = np.linspace(clim[0], clim[1], 11)
    else:
      contfs = np.linspace(clim[0], clim[1], nclev)
    if not (isinstance(cint,str) and cint=='auto'):
      contfs = np.arange(clim[0], clim[1]+cint, cint)
  elif isinstance(contfs, str) and contfs!='auto':
    use_contf=True
    use_pcol=False
    contfs = np.linspace(clim[0], clim[1], int(contfs))
  else:
    use_contf=True
    use_pcol=False
    contfs = np.array(contfs)

  ccrsdict = dict()
  if transform is not None:
    ccrsdict = dict(transform=transform)
    #adjust_axlims = False
    adjust_axlims = True
  
  # make axes if necessary
  if ax == 'auto':
    ax = plt.gca()

  # make x and y 2D
  if x.ndim==1:
    x, y = np.meshgrid(x, y)

  # convert to Basemap maps coordinates
  if bmp is not None:
    x, y = bmp(x, y)
    
  # bring x and y to correct shape for contour
  if (use_cont) or (use_contf):
    if x.shape[1] != data.shape[1]:
      xc = 0.25*(x[1:,1:]+x[:-1,1:]+x[1:,:-1]+x[:-1,:-1])
      yc = 0.25*(y[1:,1:]+y[:-1,1:]+y[1:,:-1]+y[:-1,:-1])
    else:
      xc = 1.*x
      yc = 1.*y
    
  hs = []
  # pcolor plot
  if use_pcol:
    hm = ax.pcolormesh(x, y, data, 
                        vmin=clim[0], vmax=clim[1],
                        cmap=cmap, 
                        rasterized=rasterized,
                        edgecolor=edgecolor,
                        **ccrsdict
                      )
    hs.append(hm)
  # contourf plot
  elif use_contf:
    hm = ax.contourf(xc, yc, data, contfs,
                        vmin=clim[0], vmax=clim[1],
                        cmap=cmap, 
                        extend='both',
                        **ccrsdict
                      )
  #  # this prevents white lines if fig is saved as pdf
  #  for cl in hm.collections: 
  #    cl.set_edgecolor("face")
  #    cl.set_rasterized(True)
  #  # add handle to hanlde list
  #  hs.append(hm)
  #  # rasterize
  #  if rasterized:
  #    #zorder = -5
  #    #ax.set_rasterization_zorder(zorder)
  #    for cl in hm.collections:
  #      #cl.set_zorder(zorder - 1)
  #      cl.set_rasterized(True)
  else:
    hm = None

  # extra contours
  if use_cont:
    hc = ax.contour(xc, yc, data, conts, colors=contcolor, linewidths=contlw, **ccrsdict)
    try:
      i0 = np.where(hc.levels==contthick)[0][0]
      #hc.collections[i0].set_linewidth(1.5)
      hc.collections[i0].set_linewidth(2.5*contlw)
    except:
      #print("::: Warning: Could not make contour contthick=%g thick. :::" % (contthick))
      pass
    hs.append(hc)

  # colorbar
  if ((cax is not None) and (cax!=0)) and (hm is not None): 
    if cax == 1:
      from mpl_toolkits.axes_grid1 import make_axes_locatable
      div = make_axes_locatable(ax)
      cax = div.append_axes("right", size="10%", pad=0.1)
    cb = plt.colorbar(mappable=hm, cax=cax, extend=extend)
    # this prevents white lines if fig is saved as pdf
    cb.solids.set_edgecolor("face")
    hs.append(cb)

    # colobar ticks
    if (not use_contf) and isinstance(cbticks, str) and cbticks=='auto':
      cb.formatter.set_powerlimits((-3, 2))
      tick_locator = ticker.MaxNLocator(nbins=8)
      cb.locator = tick_locator
      cb.update_ticks()
    elif use_contf:
      pass
    else:
      cb.formatter.set_powerlimits((-1, 1))
      #cb.formatter.set_scientific(False)
      cb.update_ticks()

    cax.set_title(cbtitle)

  # labels and ticks
  if adjust_axlims:
    ax.locator_params(nbins=5)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
  return hs 

#def shade(x, y, data, 
#  clims="auto", xlims="auto", ylims="auto", 
#  show_colorbar="yes", 
#  hca = "new", hcb = "new",
#  title = False, cb_title = False,
#  ys = False, frac_lower=0.5, dist_axes=0.0,
#  pcolormesh_kw=dict(cmap=plt.cm.RdYlBu_r, rasterized=True),
#  cbar_kw=dict()):
#  """ Makes a nice pcolor(mesh) plot.
#
#last change:
#----------
#2015-02-27
#  """
#  #=== some parameters
#
#  #-- ticks
#  ntickscb = 6
#
#  #-- create new axes if needed
#  if hca == "new":
#    hca = plt.axes()
#  else:
#    plt.sca(hca)
#
##  data[where(data==0.)]=nan
# 
#  #-- squeeze data
#  x = x.squeeze()
#  y = y.squeeze()
#  data = data.squeeze()
#
#  if isinstance(clims, basestring):
#  #if clims=="auto":
#    #clims=[0,0]
#    #clims[0]=np.nanmin(data)
#    #clims[1]=np.nanmax(data)
#    clims = autoclim(data, mode=clims, enab_rounded=True, enab_centered=False)
#
#  if xlims == "auto":
#    xlims=[0,0]
#    xlims[0]=np.amin(x)
#    xlims[1]=np.amax(x)
#
#  if ylims == "auto":
#    ylims=[0,0]
#    ylims[0]=np.amin(y)
#    ylims[1]=np.amax(y)
#  
#  if title != False:
#    hca.set_title(title)
#
#  if not ys:    # standard behaviour
#    hcp = plt.pcolormesh(x,y,data.transpose(),**pcolormesh_kw)
#
#    plt.clim(clims[0],clims[1])
#    plt.xlim(xlims[0],xlims[1])
#    plt.ylim(ylims[0],ylims[1])
#
#    # set ticks 
#    if True: #plt.getp(hca,'xticks').size > 0:
#      xticks = getticks(x) 
#      plt.setp(hca,'xticks',xticks)
#    if True: #plt.getp(hca,'yticks').size > 0:
#      yticks = getticks(y) 
#      plt.setp(hca,'yticks',yticks)
#
#  else:
#    # make y-vector ascending 
#    if y[1]<y[0]:
#      y = y[range(y.size-1,-1,-1)]
#      data = data[:,range(data.shape[1]-1,-1,-1)]
#
#    # position of original axes
#    pos_ori = hca.get_position().get_points()
#    xoori = pos_ori[0,0]
#    yoori = pos_ori[0,1]
#    xwori = pos_ori[1,0]-pos_ori[0,0] 
#    ywori = pos_ori[1,1]-pos_ori[0,1]
#
#    # position of lower and upper axes
#    yolo = yoori
#    ywlo = ywori*(1-dist_axes)*frac_lower
#    youp = yolo + ywlo + ywori*dist_axes
#    ywup = ywori*(1-dist_axes)*(1-frac_lower)
#
#    # change old axes to upper axes
#    hca.set_position([xoori,youp,xwori,ywup])
#    #hca_up = plt.axes([xoori,youp,xwori,ywup])
#    # create new lower axes
#    hca_lo = plt.axes([xoori,yolo,xwori,ywlo])
#
#    visible_xticks = plt.getp(hca.get_xticklabels()[0],"visible")
#    visible_yticks = plt.getp(hca.get_yticklabels()[0],"visible")
#
#    # divide data in upper and lower data
#    ind = np.argmin(np.abs(y-ys))
#    ylo = y[:ind+1]
#    yup = y[ind:]
#    dlo = data[:,:ind]
#    dup = data[:,ind:]
#
#    # plot upper data 
#    plt.axes(hca)
#    plt.pcolormesh(x,yup,dup.transpose(),**pcolormesh_kw)
#    plt.xlim(xlims[0],xlims[1])
#    plt.ylim([np.amin(yup),np.amax(yup)])
#    plt.clim(clims[0],clims[1])
#    plt.yticks(getticks(yup,nt=3)[1:])
#    plt.setp(hca.get_xticklabels(),visible=False)
#    plt.setp(hca_lo.get_yticklabels(),visible=visible_yticks)
#
#    # plot lower data
#    plt.axes(hca_lo)
#    hcp = plt.pcolormesh(x,ylo,dlo.transpose(),**pcolormesh_kw)
#    plt.xlim(xlims[0],xlims[1])
#    plt.ylim([np.amin(ylo),np.amax(ylo)])
#    plt.clim(clims[0],clims[1])
#    plt.yticks( np.concatenate((getticks(ylo,nt=3)[:-1],[ys])) )
#    plt.setp(hca_lo.get_xticklabels(),visible=visible_xticks)
#    plt.setp(hca_lo.get_yticklabels(),visible=visible_yticks)
#
#    # switch xlabel to lower axes
#    xlabel = hca.get_xlabel()
#    hca.set_xlabel("")
#    hca_lo.set_xlabel(xlabel)
#
#    plt.sca(hca)
#    
#  if hcb == "new":
#    if 'ticks' not in cbar_kw or cbar_kw['ticks'] == 'auto':
#      cbar_kw['ticks'] = np.linspace(clims[0],clims[1],ntickscb)
#    hcb = colorbar(hcp=hcp,**cbar_kw)
#  elif hcb != 0:
#    #plt.colorbar(hcp, cax=hcb, extend='both')
#    if cb_title:
#      hcb.set_title(cb_title)
#    cticks = getticks(clims)
#    plt.colorbar(hcp, cax=hcb, ticks=cticks, extend="both")
#
#  return(hcp)
#  #return(hca, hcb, hcp)
#
#def colorbar(hcp, cbtitle="", ticks="auto", hca="auto"):
#  from mpl_toolkits.axes_grid1 import make_axes_locatable
#  """
#last change:
#----------
#2014-04-28
# """ 
#
#  if hca == "auto":
#    hca = plt.gca()
#    
#  if ticks == "auto":
#    print("Error ticks=\"auto\" has not been coded yet.")
#  divider = make_axes_locatable(hca)
#  cax = divider.append_axes("right", "5%", pad="3%")
#  hcb = plt.colorbar(mappable=hcp, cax=cax, extend='both', ticks=ticks)
#  plt.title(cbtitle)
#  plt.sca(hca)
#  
#  return(hcb)

def shade2y(xp, yp, data, ys, yfrac=0.7, dyfrac=0.01, ylabel='', xlabel='', ax='auto', cax=0, clim=[None, None], cmap='auto', **kwargs):
  """ Makes a my.shade plot in with a splitted y axes. 

  ys:     y-Value where to split axes
  yfrac:  fraction of lower axes
  dyfrac: fraction of vertical distance between the two axes
  ylabel: ylabel
  xlabel: xlabel
  ax:     axes on which splitting is based on (frame for the two axes)
  cax:    axes for colorbar
  kwargs: keyword arguments for my.shade

  returns: ax1 and ax2; handles of the two axes 

last change:
----------
2017-09-08
  """
  #clims
  if isinstance(clim, str) and clim=='auto':
    clim = [None, None]
  elif isinstance(clim, str) and clim=='sym':
    clim = np.abs(data).max()
  clim=np.array(clim)
  if clim.size==1:
    clim = np.array([-1, 1])*clim
  if clim[0] is None:
    clim[0] = data.min()
  if clim[1] is None:
    clim[1] = data.max()

  if (clim[0]==-clim[1]) and cmap=='auto':
    cmap = 'RdBu_r'
  elif cmap=='auto':
    #cmap = 'viridis'
    cmap = 'RdYlBu_r'

  if isinstance(ax, str) and ax=='auto':
    ax = plt.gca()
  x0 = ax.get_position().x0
  y0 = ax.get_position().y0
  dx = ax.get_position().width
  dy = ax.get_position().height
  
  fig = plt.gcf()
  ax2 = ax
  ax1 = fig.add_axes([np.random.random(), 0.1,0.1,0.1])
  
  dy1 = dy*yfrac
  y01 = y0
  ax1.set_position([x0, y01, dx, dy1])
  
  dy2 = dy*(1-yfrac-dyfrac)
  y02 = y0+dy1 + dyfrac*dy
  ax2.set_position([x0, y02, dx, dy2])
  ax2.set_xticklabels([])
  
  ax1.set_ylim(ys[0], ys[1])
  ax2.set_ylim(ys[1], ys[2])
  
  ax1.set_yticks(np.linspace(ys[0], ys[1], 5))
  ax2.set_yticks(np.linspace(ys[1], ys[2], 3)[1:])
  ax1.set_xlim(xp.min(), xp.max())
  ax2.set_xlim(xp.min(), xp.max())

  hm = shade(xp, yp, data, ax=ax1, cax=cax, adjust_axlims=False, clim=clim, cmap=cmap, **kwargs) 
  hm = shade(xp, yp, data, ax=ax2, cax=cax, adjust_axlims=False, clim=clim, cmap=cmap, **kwargs)
  
  hylab = ax1.text(-0.12, 0.5/yfrac, ylabel, transform=ax1.transAxes, ha='center', va='center', rotation=90)
  hxlab = ax1.set_xlabel(xlabel)
  return ax1, ax2

def trishade(Tri, datai,
            ax='auto', cax=0,
            cmap='auto',
            rasterized=True,
            clim=[None, None],
            extend='both',
            edgecolor='none',
            conts=None,
            nclev='auto',
            cint='auto',
            contcolor='k',
            contthick=0.,
            contfs=None,
            contlw=1.,
            use_pcol=True,
            adjust_axlims=True,
            bmp=None,
            transform=None,
            logplot=False,
         ):
  """ Makes a nice tripcolor plot.

last change:
----------
2018-03-08
  """

  # mask 0 and negative values in case of log plot
  data = 1.*datai
  if logplot and isinstance(data, np.ma.MaskedArray):
    data[data<=0.0] = np.ma.masked
    data = np.ma.log10(data) 
  elif logplot and not isinstance(data, np.ma.MaskedArray):
    data[data<=0.0] = np.nan
    data = np.log10(data) 

  #clims
  if isinstance(clim, str) and clim=='auto':
    clim = [None, None]
  elif isinstance(clim, str) and clim=='sym':
    clim = np.abs(data).max()
  clim=np.array(clim)
  if clim.size==1:
    clim = np.array([-1, 1])*clim
  if clim[0] is None:
    clim[0] = data.min()
  if clim[1] is None:
    clim[1] = data.max()

  if (clim[0]==-clim[1]) and cmap=='auto':
    cmap = 'RdBu_r'
  elif cmap=='auto':
    #cmap = 'viridis'
    cmap = 'RdYlBu_r'

  # calculate contour x/y and contour levels if needed
  if conts is None:
    use_cont = False
  elif isinstance(conts,str) and conts=='auto':
    use_cont = True
    if isinstance(nclev,str) and nclev=='auto':
      conts = np.linspace(clim[0], clim[1], 11)
    else:
      conts = np.linspace(clim[0], clim[1], nclev)
    if not (isinstance(cint,str) and cint=='auto'):
      conts = np.arange(clim[0], clim[1]+cint, cint)
  else:
    use_cont = True
    conts = np.array(conts)

  if contfs is None:
    use_contf=False
  elif isinstance(contfs, str) and contfs=='auto':
    use_contf=True
    use_pcol=False
    if isinstance(nclev,str) and nclev=='auto':
      contfs = np.linspace(clim[0], clim[1], 11)
    else:
      contfs = np.linspace(clim[0], clim[1], nclev)
    if not (isinstance(cint,str) and cint=='auto'):
      contfs = np.arange(clim[0], clim[1]+cint, cint)
  elif isinstance(contfs, str) and contfs!='auto':
    use_contf=True
    use_pcol=False
    contfs = np.linspace(clim[0], clim[1], int(contfs))
  else:
    use_contf=True
    use_pcol=False
    contfs = np.array(contfs)

  ccrsdict = dict()
  if transform is not None:
    ccrsdict = dict(transform=transform)
    #adjust_axlims = False
    adjust_axlims = True
  
  # make axes if necessary
  if ax == 'auto':
    ax = plt.gca()

  #### make x and y 2D
  ###if x.ndim==1:
  ###  x, y = np.meshgrid(x, y)

  #### convert to Basemap maps coordinates
  ###if bmp is not None:
  ###  x, y = bmp(x, y)
  ###  
  #### bring x and y to correct shape for contour
  ###if (use_cont) or (use_contf):
  ###  if x.shape[1] != data.shape[1]:
  ###    xc = 0.25*(x[1:,1:]+x[:-1,1:]+x[1:,:-1]+x[:-1,:-1])
  ###    yc = 0.25*(y[1:,1:]+y[:-1,1:]+y[1:,:-1]+y[:-1,:-1])
  ###  else:
  ###    xc = 1.*x
  ###    yc = 1.*y
    
  hs = []
  # pcolor plot
  if use_pcol:

    hm = ax.tripcolor(Tri, data, 
                        edgecolor=edgecolor,
                        vmin=clim[0], vmax=clim[1],
                        cmap=cmap, 
                        rasterized=rasterized,
                        **ccrsdict
                      )
    hs.append(hm)
  # contourf plot
  elif use_contf:
    hm = ax.contourf(xc, yc, data, contfs,
                        vmin=clim[0], vmax=clim[1],
                        cmap=cmap, 
                        extend='both',
                        **ccrsdict
                      )
    # this prevents white lines if fig is saved as pdf
    for cl in hm.collections: 
      cl.set_edgecolor("face")
    # add handle to hanlde list
    hs.append(hm)
  else:
    hm = None

  # extra contours
  if use_cont:
    hc = ax.contour(xc, yc, data, conts, colors=contcolor, linewidths=contlw, **ccrsdict)
    try:
      i0 = np.where(conts==contthick)[0][0]
      #hc.collections[i0].set_linewidth(1.5)
      hc.collections[i0].set_linewidth(2.5*contlw)
    except:
      #print("::: Warning: Could not make contour contthick=%g thick. :::" % (contthick))
      pass
    hs.append(hc)

  # colorbar
  if ((cax is not None) and (cax!=0)) and (hm is not None): 
    if cax == 1:
      from mpl_toolkits.axes_grid1 import make_axes_locatable
      div = make_axes_locatable(ax)
      cax = div.append_axes("right", size="10%", pad=0.1)
    cb = plt.colorbar(mappable=hm, cax=cax, extend=extend)
    # this prevents white lines if fig is saved as pdf
    cb.solids.set_edgecolor("face")
    hs.append(cb)

    # colobar ticks
    cb.formatter.set_powerlimits((-3, 2))
    tick_locator = ticker.MaxNLocator(nbins=8)
    cb.locator = tick_locator
    cb.update_ticks()

  # labels and ticks
  if adjust_axlims:
    ax.locator_params(nbins=5)
    ax.set_xlim(Tri.x.min(), Tri.x.max())
    ax.set_ylim(Tri.y.min(), Tri.y.max())
  return hs 


def getticks(x, nt=5):
  """
last change:
----------
2014-05-16
 """ 
  # find min and max from x
  x = np.array(x)
  x = x[np.isnan(x)==0]
  x1 = np.amin(x)
  x2 = np.amax(x)
  
  # round intervall and starting value
  dx = x2-x1
  exp = np.ceil(np.log10(dx))
  intx = np.fix( dx/10**(exp-2) )*10**(exp-2)/(nt-1)
  xs   = np.fix( x1/10**(exp-2) )*10**(exp-2)

  # calc ticks
  ticks = np.arange(xs,x2+0.5*intx,intx)
  return ticks

def plotsettings(hca, hcb=[], fig_size_fac=1):
  """
  fsaxtic     font size of axes ticks
  fsaxlat     font size of axes label
  fsaxtic     font size of axes title
  fsaxlab     font size of axes lable [ (a), (b), (c), etc. ]
  fscbtit     font size of colorbar title
  fscbtic     font size of colorbar ticks
last change:
----------
2015-07-22
 """ 
  #-- fontsize axis
  fsaxtit = fsaxlab = 9.
  fsaxtic = 7.
  #-- fontsizes colorbar
  fscbtit = fscbtic = 7.
  #-- factor to increase distance between ticks and label
  ylabelpad_fac = 1.0


  # apply fig_size_fac
  fsaxtit = fsaxtit*fig_size_fac
  fsaxlab = fsaxlab*fig_size_fac
  fsaxtic = fsaxtic*fig_size_fac
  fscbtit = fscbtit*fig_size_fac
  fscbtic = fscbtic*fig_size_fac
  ylabelpad_fac = ylabelpad_fac*fig_size_fac
  
  # apply changes axes
  for ax in hca:
    # title
    ax.title.set_fontsize(fsaxtit)
    # axis labels
    ax.xaxis.label.set_fontsize(fsaxlab)
    ax.yaxis.label.set_fontsize(fsaxlab)
    # ticks
    for tt in ax.get_xticklabels()+ax.get_yticklabels():
        tt.set_fontsize(fsaxtic)
    # offset or multiplyer at axes (e.g. 1e9)
    ax.xaxis.get_offset_text().set_fontsize(fsaxtic)
    ax.yaxis.get_offset_text().set_fontsize(fsaxtic)
    # distance between ticks and labels of y axis
    ax.yaxis.labelpad = ylabelpad_fac * ax.yaxis.labelpad
    # label of axes (e.g. [(a), (b), etc.])
    if hasattr(ax, 'axlab'):
      ax.axlab.set_fontsize(fsaxlab)
  
  # apply changes colorbar
  for cax in hcb:
    if cax != 0:
      cax.title.set_fontsize(fscbtit)
      cax.yaxis.get_offset_text().set_fontsize(fsaxtic)
      for tt in cax.get_yticklabels():
          tt.set_fontsize(fscbtic)

def arrange_axes( nx,ny,
                  # height of and aspect ratio of subplot
                  asy  = 3.5,
                  sasp = 1.0,
                  # plot colorbar
                  plot_cb = False,
                  # have common x or y axes
                  sharex = True, sharey = True,
                  xlabel = "",   ylabel = "",
                  # additional space left right and above and below axes
                  oxl = 0.1, oxr = 0.0,
                  oyb = 0.0, oyt = 0.0,
                  # factor that increases distance between axes
                  axfac_x = 1., axfac_y = 1.,
                  # kw for axes labels [(a), (b), etc.]
                  axlab_kw = dict(),
                  # figure size and aspect ratio
                  fig_size     = 'auto',
                  fig_asp      = 'auto',
                  fig_size_fac = 1.,
                  # figure title
                  fig_title = None,
                  projection = None,
                  ):
  """
last change:
----------
2015-07-22
 """ 

  # all lengths are in cm
  cm2inch = 0.3937        # to convert cm into inch

  # horizontal standard spaces
  alx = 1.0
  asx = asy / sasp
  adx = 0.5    
  cdx = 0.2
  clx = 0.8
  csx = 0.32 
  
  # vertical standard spaces
  aly = 0.8
  asy = asy
  ady = 0.2  
  aty = 0.6
  fty = 1.               # extra space for figure title (set to zero if fig_title = None)

  # apply manual changes to spaces
  adx = adx * axfac_x 
  ady = ady * axfac_y 
  #cdx = cdx * axfac_x   # this is a fix I do not understand why adxv is set to cdx if icbspace==True
  clx = clx * axfac_x

  if fig_title==None:
    fty = 0.

  # make vector of plot_cb if it has been true or false before
  # plot_cb can have values [{1}, 0] 
  # with meanings:
  #   1: plot cb; 
  #   0: do not plot cb
  if isinstance(plot_cb, bool) and (plot_cb==True):
    plot_cb = np.ones((nx,ny))  
    nohcb = False
  elif isinstance(plot_cb, bool) and (plot_cb==False):
    plot_cb = np.zeros((nx,ny))
    nohcb = True
  else:
    plot_cb = np.array(plot_cb)
    if plot_cb.size!=nx*ny:
      raise ValueError('Vector plot_cb has wrong length!')
    if plot_cb.shape[0]==nx*ny:
      plot_cb = plot_cb.reshape(ny,nx).transpose()
    elif plot_cb.shape[0]==ny:
      plot_cb = plot_cb.transpose()
    nohcb = False

  if not isinstance(projection, list):
    projection = [projection]*nx*ny

  # make spaces vectors
  # horizontal
  alxv = np.array([alx]*(nx))
  asxv = np.array([asx]*(nx))
  adxv = np.array([adx]*(nx))
  clxv = np.array([clx]*(nx))
  csxv = np.array([csx]*(nx))

  icbspace = plot_cb.sum(axis=1)>0
  csxv[icbspace==False] = 0.0
  clxv[icbspace==False] = 0.0
  adxv[icbspace==True ] = cdx
  if sharey:
    alxv[1:] = 0.0  

  # vertical
  alyv = np.array([aly]*(ny))
  asyv = np.array([asy]*(ny))
  adyv = np.array([ady]*(ny))
  atyv = np.array([aty]*(ny))

  if sharex:
    alyv[:-1] = 0.0

  # calculate figure size
  fw_auto = ( oxl + (alxv+asxv+adxv+csxv+clxv).sum() + oxr       )
  fh_auto = ( oyb + (alyv+asyv+adyv+atyv).sum()      + oyt + fty )
  if fig_size == 'auto':
    fw = fw_auto 
    fh = fh_auto 
  elif fig_size == 'dina4pt':
    fw = 21.0
    fh = 29.7
  elif fig_size == 'dina4ls':
    fw = 29.7
    fh = 21.0
  elif fig_size == 'jpo':
    fw = 15.5
    if fig_asp == 'auto':
      fh = fh_auto
    else:
      fh = fw*fig_asp
  elif isinstance( fig_size, (int,float) ):
    fw = fig_size
    if fig_asp == 'auto':
      fh = fh_auto
    else:
      fh = fw*fig_asp

  # make figure
  fasp = fh/fw
  hcf = plt.figure(figsize=(fw*cm2inch*fig_size_fac, fh*cm2inch*fig_size_fac))

  if not fig_title == None:
    hcf.suptitle(fig_title)

  # handle for axes
  hca = [0]*(nx*ny) 
  hcb = [0]*(nx*ny)

  kk = -1
  for jj in range(ny):
    for ii in range(nx):
      kk += 1

      # set axes x offspring
      if ii == 0:
        oxa = oxl + alxv[ii]
      else:
        oxa = oxa + alxv[ii] + (asxv+adxv+csxv+clxv)[ii-1]

      # set axes y offsping
      #if jj == 0 and ii == 0:
      #  oya = oyb + alyv[jj]
      #elif jj != 0 and ii == 0:
      #  oya = oya + alyv[jj] + (asyv+adyv+atyv)[jj-1]

      if jj == 0 and ii == 0:
        oya = fh - oyt - fty - (atyv+asyv)[jj]
      elif jj != 0 and ii == 0:
        oya =      oya - alyv[jj-1] - (adyv+atyv+asyv)[jj]

      # set colorbar x offspring
      oxc = oxa + (asxv+adxv)[ii]

      # calculated rectangles for axes and colorbar
      rect   = np.array([oxa, oya/fasp, asxv[ii], asyv[jj]/fasp])/fw
      rectcb = np.array([oxc, oya/fasp, csxv[ii], asyv[jj]/fasp])/fw
      
      # plot axes
      if projection[kk] is None:
        hca[kk] = plt.axes(rect, xlabel=xlabel, ylabel=ylabel)
      else:
        hca[kk] = plt.axes(rect, xlabel=xlabel, ylabel=ylabel, projection=projection[kk])

      # delet labels for shared axes
      if sharex and jj!=ny-1:
        hca[kk].ticklabel_format(axis='x',style='plain',useOffset=False)
        hca[kk].tick_params(labelbottom='off')
        hca[kk].set_xlabel('')

      if sharey and ii!=0:
        hca[kk].ticklabel_format(axis='y',style='plain',useOffset=False)
        hca[kk].tick_params(labelleft='off')
        hca[kk].set_ylabel('')

      # plot colorbars
      if plot_cb[ii,jj] == 1:
        hcb[kk] = plt.axes(rectcb, xticks=[])
        hcb[kk].yaxis.tick_right()

  # add letters for subplots
  if axlab_kw is not None:
    hca = axlab(hca, fontdict=axlab_kw)
  
  # return axes handles
  if nohcb:
    #plotsettings(hca)
    return hca, hcb
  else:
    #plotsettings(hca,hcb)
    return hca, hcb

def axlab(hca, figstr=[], posx=[-0.0], posy=[1.08], fontdict=None):
  """
input:
----------
  hca:      list with axes handles
  figstr:   list with strings that label the subplots
  posx:     list with length 1 or len(hca) that gives the x-coordinate in ax-space
  posy:     list with length 1 or len(hca) that gives the y-coordinate in ax-space
last change:
----------
2015-07-21
  """

  # make list that looks like [ '(a)', '(b)', '(c)', ... ]
  if len(figstr)==0:
    lett = "abcdefghijklmnopqrstuvwxyz"
    lett = lett[0:len(hca)]
    figstr = ["z"]*len(hca)
    for nn, ax in enumerate(hca):
      figstr[nn] = "(%s)" % (lett[nn])
  
  if len(posx)==1:
    posx = posx*len(hca)
  if len(posy)==1:
    posy = posy*len(hca)
  
  # draw text
  for nn, ax in enumerate(hca):
    ht = hca[nn].text(posx[nn], posy[nn], figstr[nn], 
                      transform = hca[nn].transAxes, 
                      horizontalalignment = 'right',
                      fontdict=fontdict)
    # add text handle to axes to give possibility of changing text properties later
    # e.g. by hca[nn].axlab.set_fontsize(8)
    hca[nn].axlab = ht
  return hca

def autoclim(data, mode=None, enab_rounded=True, enab_centered=False):
  """ determine color limits for pcolor-like plots 

input:
----------
*data* : array     
  Data array for which clims should be calculated.
*mode* : string
  Mode decides how color levels are calculated:
  minmax:     simple minimum and maximum (default)
  cminmax:    zero centered minimum and maximum
*enab_round* : bool
  If True, clims are rounded to three significant digits

output:
----------
*clims* : np.array 
  calculated color levels

last change:
----------
2015-02-27
  """ 
  if type(data) is np.ma.core.MaskedArray:
    data = data.flatten()
    #data = np.ma.filled( data[data.mask==False] )
    data = np.ma.filled( data, fill_value=0.0 )
  else:
    data = data.flat
    data = data[np.isnan(data)==False]
  data = data[np.isinf(np.abs(data))==False]
  data = data[data!=0]
  clim = np.zeros((2,))

  if mode is None or mode == "auto":
    mode = "minmax"

  if   mode == "minmax":
    clim[0] = data.min()
    clim[1] = data.max()
  elif mode == "2std":
    #data = np.sort(data)
    #data = data[0:np.size(data)*9/10]
    mean_data = np.mean(data)
    std_data = np.std(data)
    clim[0] = mean_data - 2*std_data
    clim[1] = mean_data + 2*std_data
    
  if enab_centered:
    absmax = np.max( [np.abs(clim[0]), np.abs(clim[1])] )
    clim[0] = - absmax
    clim[1] =   absmax

  if enab_rounded:
    clim[0] = roundsign(clim[0],2)
    clim[1] = roundsign(clim[1],2)

  return clim

def figinfo(scriptname='', add_text='', tp=(0.98,0.02), fontsize=6, del_stdtext=False, hcf='auto'):
  """ write info line containing scriptname and date to figure
last change:
----------
2015-08-26
  """
  if hcf == 'auto':
    hcf = plt.gcf()

  date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
  text = 'Created with %s at %s' % (scriptname, date)

  if del_stdtext:
    text = ''

  if add_text != '':
    text = text + '\n' + add_text
  ht = hcf.text(tp[0],tp[1], text, fontsize=fontsize, ha='right', va='bottom')
  return ht

def dcmap(N=20,cmap=None):
  """ return colormap cmap with N discret levels and an under and over color-level set

example usage:
----------
dcmap = my.dcmap(N=20,cmap='jet')
hm    = ax.pcolormesh(X,Y,Z, cmap=dcmap, vmin=vmin, vmax=vmax)            
hcb   = plt.colorbar(ticks=np.linspace(vmin,vmax,N/2+1), mappable=hm, extend='both')
 
last change:
----------
2015-08-31
  """
  base = plt.cm.get_cmap(cmap)                                                      
  color_list = base(np.linspace(0, 1, N+2))                                         
  cmap_name  = base.name + str(N)                                                   
  new_cmap   = base.from_list(cmap_name, color_list[1:-1], N)                       
  new_cmap.set_over(color_list[-1])                                                 
  new_cmap.set_under(color_list[0])
  return new_cmap


def saveshow(add_figinfo=True, nfig=False, fig_name='auto', fig_format='png', path_fig='auto', savefig='auto'):
  """ either save or show a figure 

example usage:
----------
saveshow()
 
last change:
----------
2015-10-14
  """
  # load extra module
  if fig_name == 'auto' or add_figinfo:
    from inspect import stack
    caller_frame = stack()[1]
    calling_script = caller_frame[0].f_globals.get('__file__', None)
  
  # find auto fig_name
  if fig_name == 'auto':
    #fig_name = __file__.split('/')[-1][:-3]
    fig_name = calling_script[:-3]

  # find figure number
  if nfig==False:
    nfig=''
  elif isinstance(nfig, int):
    nfig = '_%02d' % (nfig)
  
  # check if figure should be saved and where it should be saved
  if savefig == 'auto':
    f = open('do_save_fig')
    lines = f.readlines()
    savefig = int(lines[0])
    path_fig = lines[1][:-1]
    f.close()
    if savefig==1:
      savefig=True
    else:
      savefig=False
  
  # find path to save figure
  if path_fig == 'auto':
    f = open('do_save_fig')
    lines = f.readlines()
    path_fig = lines[1][:-1]
    f.close()
  
  # add fig info
  if add_figinfo:
    #figinfo(__file__.split('/')[-1])
    figinfo(calling_script)
  
  # save or show figure
  if savefig:
    fname = "%s/%s%s.%s" % (path_fig, fig_name, nfig, fig_format)
    print('Saving figure %s' % (fname))
    plt.savefig(fname, dpi=400)
  else:
    plt.show()
  return

def patch(xp, yp, ax, fc='none', ec='k', lw=1, alpha=1):
  xp = np.array(xp)
  yp = np.array(yp)
  verts = np.concatenate((xp.reshape(-1,1),yp.reshape(-1,1)), axis=1)
  verts = np.concatenate( (verts, verts[0:1,:]), axis=0 )
  codes = [Path.MOVETO,
           Path.LINETO,
           Path.LINETO,
           Path.LINETO,
           Path.CLOSEPOLY,
           ]
  path = Path(verts, codes)
  patch = patches.PathPatch(path, fc=fc, lw=lw, ec=ec, alpha=alpha)
  ax.add_patch(patch)
  return patch

def make_inlay_axes(xy='auto', width='auto', height='auto', asp='auto', 
                    points='auto',
                    transform='auto', fig='auto', ax='auto'):
  """ Transforms coordinates of a rectangle from either data or axes coordinates into figure coordinates.
  
  ---
  Input:
  ---
  Either use:
  xy          [left, bottom] of axes rectangle
  width       width of axes rectangle
  height      height of axes rectangle
  asp         aspect ratio, if asp is given height or width has to be 'auto' and are calculated according to asp
  or:
  points      coordinates for rectangle in axes or data coordinates
              if points.ndim==1, points is assumed to be [left, bottom, width, height] of the rectangle 
              (same as specifying xy, width and height; kept for backward compatibility)
              if points.ndim==2, points is assumed to be lower left and upper right points of the rectangle

  transform   either ax.transAxes if points is in axes coordinates or ax.transData if it is in Data coordinates
              other possibilities might also be possible (default is ax.transAxes)
              (transform is a matplotlib.transforms.BboxTransformTo object)
  fig         figure for which the transformation is done (default is plt.gcf())
  ax          axes for which the transformation is done (default is plt.gca())

  ---
  Output:
  ---
  rect_fig    recatanlge in figure coordinates
              (to be used e.g. for creating inlay axes by ai = hcf.add_axes(rect_new))

  ---
  Example:
  ---
  rect = [0.1, 0.2, 0.8, 0.2]       # rectangle in axes coordinates
  rect_fig = trans_inlay_rect(rect)  # rectangle in figure coordinates
  ai = hcf.add_axes(rect_fig)       # create axes inlay with position needed in figure coordinates
  """
  # calculate height or width if asp is given
  if (xy!='auto') and (asp!='auto'):
    if (width=='auto') and isinstance(height, float):
      width  = height/asp
    elif (height=='auto') and isinstance(width, float):
      height = width*asp
    else:
      raise ValueError('::: Error: width or height need to be specified together with asp!:::')
    
  # convert xy into points if xy etc. are given
  if (xy!='auto') and ((width=='auto') or (height=='auto')):
    raise ValueError('::: Error: width and height has to be specified together with xy!:::')
  if xy!='auto':
    if points!='auto':
      print("::: Warning, points in make_inlay_axes is ignorred since xy etc. is specified!")
    #points = [xy[0], xy[1], xy[0]+width, xy[1]+height]
    points = [xy[0], xy[1], width, height]

  # convert points to lower left and upper right points if input was [left, bottom, width, height] vector
  points = np.array(points)
  if points.ndim==1:
    rect=points
    points = np.array([[rect[0], rect[1]],
                    [rect[0]+rect[2], rect[1]+rect[3]]])

  # set defaults
  if fig=='auto':
    fig = plt.gcf()
  if ax=='auto':
    ax = plt.gca()
  if transform=='auto':
    transform = ax.transAxes

  # transform from data or axes to display coordinates
  points_dis = transform.transform(points)
  # transform from display coordinates to figure coordinates
  points_fig = fig.transFigure.inverted().transform(points_dis)

  # calculate new rectangle with [left, bottom, width, height] from lower left and upper right points
  rect_fig = [points_fig[0,0], points_fig[0,1], points_fig[1,0]-points_fig[0,0], points_fig[1,1]-points_fig[0,1]]
  
  ax = fig.add_axes(rect_fig)
  return ax

# ================================================================================ 
# map plotting functions (need to import mpl_toolkits.basemap, matplotlib.pyplot and numpy)
# ================================================================================ 
def plot_map(lon_reg, lat_reg,
              dlon='auto', dlat='auto',
              resolution='auto',
              color_land='#EEE8AA',
              projection='mill',
              drawcoastlines=True,
              drawgrid=True,
              coastlinecolor='k',
              coastlinelinewidth=1.0,
              grid_lw=1.,  # good values 0. (no lines) or 1. (thin lines)
              labelstyle='', # '' with E and W '+/-'
              ax='auto',
              bmpkw=dict(),
              fix_aspect=True,
              ):
  try:
    from mpl_toolkits.basemap import Basemap
  except:
    raise ValueError('::: Error: Could not load Basemap! :::')

  if ax=='auto':
    ax = plt.gca()

  lon_reg = np.array(lon_reg)
  lat_reg = np.array(lat_reg)

  if dlon=='auto':
    dl = np.abs(lon_reg[1]-lon_reg[0])
    if dl<12.:
      dll = 2.
    elif dl<30.:
      dll = 5.
    elif dl<60:
      dll = 10.
    elif dl<180.:
      dll = 30.
    else:
      dll = 60
    dlon = dll

  if dlat=='auto':
    dl = np.abs(lat_reg[1]-lat_reg[0])
    if dl<12.:
      dll = 2.
    elif dl<30.:
      dll = 5.
    elif dl<60:
      dll = 10.
    elif dl<180.:
      dll = 30.
    else:
      dll = 60.
    dlat = dll

  if resolution=='auto':
    dl = np.amin([dlon, dlat])
    if dl<=5.:
      resolution='l'
    elif dl<=10 :
      resolution='l'
    else:
      resolution='c'

  #print('dl=%d, resolution = %s' % (dl, resolution))
  bmp = Basemap(projection=projection,
                llcrnrlon=lon_reg[0], urcrnrlon=lon_reg[1],
                llcrnrlat=lat_reg[0], urcrnrlat=lat_reg[1],
                lat_0=lat_reg.mean(), lon_0=lon_reg.mean(),
                resolution=resolution,
                ax=ax,
                fix_aspect=fix_aspect,
                **bmpkw
               )
  if drawcoastlines:
    bmp.drawcoastlines(color=coastlinecolor, linewidth=coastlinelinewidth)
  if color_land is not None:
    bmp.fillcontinents(color=color_land)
  if drawgrid:
    if grid_lw == 0.0:
      grid_lw = 1e-3
      print("Bad fix for grid_lw=0.0. Hope this is not necessary anymore in the future.")
    par = bmp.drawparallels(np.arange(-90.,91.,dlat),   labels=[1,0,0,0], linewidth=grid_lw, labelstyle=labelstyle)
    mer = bmp.drawmeridians(np.arange(-180.,181.,dlon), labels=[0,0,0,1], linewidth=grid_lw, labelstyle=labelstyle)
  #plot_bathymetry(bmp)
  return bmp

# ================================================================================ 
# math functions (need to import numpy)
# ================================================================================ 

def vectogrid(a,dims):
  """
last change:
----------
2014-05-05
 """ 
  
  # find index where a should not be reproduced 
  sd = dims.index(1)

  # reshape a
  A = a.squeeze()
  if len(dims) == 3:
    if sd == 0:
      tp = [2,0,1]
    elif sd == 1:
      tp = [1,2,0]
    elif sd == 2:
      tp = [0,1,2]
  if len(dims) == 2:
    if sd == 0:
      tp = [1,0]
    elif sd == 1:
      tp = [0,1]
  A = np.tile(np.kron(a,np.ones((1,)*len(dims))).transpose(tp[:]),dims)

  return(A)

def stats(a, formatstr = "%5.2g"):
  """
last change:
----------
2014-05-07
 """ 
  
  type_a = type(a)
  shape_a = "("+"%d,"*len(a.shape) % (a.shape)+")"
  numnans_a = np.isnan(a).sum()
  if numnans_a == a.size:
    print("a contains only nans")
    print("a.shape: %s" % (shape_a))
    return

  mask_a = np.isnan(a) == False
  max_a = np.max(a[mask_a])
  min_a = np.min(a[mask_a])
  mean_a = np.mean(a[mask_a])
  absmin_a = np.min(np.abs(a[mask_a]))
  absmax_a = np.max(np.abs(a[mask_a]))

  strr = "a.shape: %s; num of NaNs: %s; a.max: %s; a.min: %s; a.abs.max %s; a.abs.min: %s" \
    % ("%s", "%d", formatstr, formatstr, formatstr, formatstr)
  print(strr % (shape_a, numnans_a, max_a, min_a, absmax_a, absmin_a))

def nanmean(a, **kw):
  """
last change:
----------
2014-05-07
 """ 
  mask = np.isnan(a)==False
  an = a
  an[np.isnan(a)] = 0.0 
  nansum_a = np.sum(an,**kw)
  nansum_mask = np.sum(np.float32(mask),**kw)
  nansum_mask[nansum_mask==0.0] = np.nan
  
  nanmean_a = np.zeros(a.shape)
  nanmean_a = nansum_a/nansum_mask
  nanmean_a[np.isnan(nanmean_a)] = 0.0

  return(nanmean_a)

def bootstrap(x, n=None, nb=10000, alpha=0.05, lowhigh=False):
  """
input:
----------
x:        data
n:        length of bootstrap sample (standard: x.size)
nb:       number of bootstrap iterations 
alpha:    percentile to use for confidence interval
lowhigh:  if True min and max of bootstrap mean confidencence level is calculated
          if False np.max(np.abs(np.mean(x)-[min, max])) is calculated

output:
----------
bootstrap "alpha"-confidence percentile

last change:
----------
2010-06-10

  """
  if n == None:
    n = x.size
  if nb == None:
    nb = x.size
  
  xm_samp = np.zeros(nb)
  xs_samp = np.zeros(nb)
  for i in range(nb):
    ind = np.floor(n*np.random.rand(n)).astype(int)
    xm_samp[i] = np.mean(x[ind])

  xm_samp.sort()

  alval = np.array([alpha/2,1-alpha/2])
  indal = np.round(alval*nb).astype(int)
  xm_min = xm_samp[indal[0]]
  xm_max = xm_samp[indal[-1]]

  if lowhigh:
    bs = np.array([xm_min,xm_max])
  else:
    bs = np.max(np.abs(np.mean(x)-[xm_min,xm_max]))

  return bs

def repmat(x,rdim):
  """ matlab repmat-like function

input:
----------
x:      matrix that should be copied
rdim:   vector that states the number of copies along the respective dimension

output:
----------
X:      extended matrix

last change:
----------
2014-06-17
  """
  rdim = np.array(rdim)
  if rdim.size == x.ndim:
    X = np.tile(x,rdim)
  else:
    X = np.tile(x,rdim[::-1]).transpose(range(rdim.size)[::-1])
  return X

def roundsign(a, signum):
  """ round a to sign significant digits

intput:
----------
*a* : float
  value that should be rounded
*signum* : integer
  number of significant digits

output:
----------
*asig* : float
  rounded value

last change:
----------
2015-02-26
  """
  signum = int(signum)
  power = np.ceil(np.log10(np.abs(a)))
  asig = np.around(a / 10**(power-signum)) * 10**(power-signum)

  return asig

def indfind(elements, vector):                                                      
  """ return indices of elements that closest match elements in vector
last change:
----------
2015-08-31
  """
  # convert elements to np array                                                    
  if type(elements) is int or type(elements) is float:                              
    elements = np.array([elements])
  elif type(elements) is list:
    elements = np.array(elements)                                                   
  # find indices
  inds = [0]*elements.size                                                          
  for i in range(elements.size):                                                    
    inds[i] = np.argmin( np.abs(vector - elements[i]) )                             
  return inds

def bilininterp(x, y, z, xi, yi):
  """ bilinear interpolation from z[y,x] to points xi, yi, where x and y are assumed to build a regular grid
if x and y are not vectors they are reduced to x[0,:] and y[:,0]
x and y probably should not have any missing values
missing values in z, xi and yi are ok

last change:
----------
2015-11-18
  """
  # if x and y are grids take
  if x.ndim==2 and y.ndim==2:
    x = x[0,:]
    y = y[:,0]
  if (y.size, x.size) != z.shape:
    raise ValueError('::: Error: (y.size, x.size) has to be equal to z.shape! :::')

  # allocate zi
  zi = np.ma.array( 0.*xi, mask=0.*xi==0 )

  # detect "float"-indices such that xi-x0=dx*ie and yi-y0=dy*je
  ie = np.ma.array( (xi - x[0])/(x[1]-x[0]) )
  je = np.ma.array( (yi - y[0])/(y[1]-y[0]) )
  
  # remove points that are outside of grid
  indout = (ie>=x.size-1) | (ie<=0) | (je>=y.size-1) | (je<=0)
  ie = ie[~indout]
  je = je[~indout]
  #ie[np.ma.where( (ie>=x.size-1) | (ie<=0) )] = np.ma.masked
  #je[np.ma.where( (je>=y.size-1) | (je<=0) )] = np.ma.masked
  xi = xi[~indout]
  yi = yi[~indout]

  # allocate index matrices as integers
  i1 = np.empty_like(ie, dtype=np.int32)
  i2 = np.empty_like(ie, dtype=np.int32)
  j1 = np.empty_like(ie, dtype=np.int32)
  j2 = np.empty_like(ie, dtype=np.int32)

  # calculate edge points to xi and yi
  i1 = np.floor(ie, i1.astype(float))
  i2 = np.ceil( ie, i2.astype(float))
  j1 = np.floor(je, j1.astype(float))
  j2 = np.ceil( je, j2.astype(float))

  i1 = i1.astype(int)
  i2 = i2.astype(int)
  j1 = j1.astype(int)
  j2 = j2.astype(int)

  # interpolate first from   f(x1,y1) and f(x2,y1) to f(xi,y1)
  # interpolate than from    f(x1,y2) and f(x2,y2) to f(xi,y2)
  # finally interpolate from f(xi,y1) and f(xi,y2) to f(xi,yi) 
  zii = 1./((x[i2]-x[i1])*(y[j2]-y[j1])) * ( \
       z[j1,i1]*(x[i2]-xi)*(y[j2]-yi) + \
       z[j1,i2]*(xi-x[i1])*(y[j2]-yi) + \
       z[j2,i1]*(x[i2]-xi)*(yi-y[j1]) + \
       z[j2,i2]*(xi-x[i1])*(yi-y[j1]) )

  zi[~indout] = zii
  return zi

def find1st(cond, axis=0):
  """ Find 1st index of n-dimensional condition matrix 'cond' along axis 'axis' 

  Returns array 'ind1st' of indices as np.where, np.nonzero etc. but only 
  containing the first match along dimension

last change:
----------
2015-12-03
  """
  ndim = cond.ndim

  icond  = np.nonzero(cond)
  ind1st = [0]*ndim
  ind    = np.nonzero( np.diff( np.insert(icond[axis],0,-1) ) )[0]

  for nn in range(ndim):
    ind1st[nn] = icond[nn][ind]
  return ind1st

def sort_contour(xps, yps, ind=0, ds=1e-4):
  """ Sort contour points counter clockwise.

  Usage: 
  xps, yps = my.sort_contour(xps, yps, ind=0, ds=1e-4)

last change:
----------
2016-12-09
  """
  # get center of point
  xpc = xps.mean()
  ypc = yps.mean()
  # make input contour arrays masked arrays
  xps = np.ma.array(xps)
  yps = np.ma.array(yps)
  # allocate new contour arrays
  xpn = np.zeros(xps.shape)
  ypn = np.zeros(yps.shape)
  for nn in range(xps.size):
    # get contour points
    xp = xps[ind]
    yp = yps[ind]
    # save them in new arrays
    xpn[nn] = xp
    ypn[nn] = yp
    # set them to missing such that they are not found again
    xps[ind] = np.ma.masked
    yps[ind] = np.ma.masked
    # search for next point that is closest to infinitesimal tangent guess
    vpx = -ypc+yp
    vpy =  xpc-xp
    xg = xp+vpx*ds
    yg = yp+vpy*ds
    ind = ((xps-xg)**2+(yps-yg)**2).argmin()
  return xpn, ypn

def get_island_mask(mask, iis, jjs, 
                   maxit='auto', stencil='four_point', verbose=False):
  """ Get all points of an island (and not more) within a land mask.
  
  mask: mask in which island is identified.
  iis:  i-index of startpoint that is a mask point or has mask points around
  jjs:  j-index of startpoint that is a mask point or has mask points around

  Different stencils can be chosen:
    stencil='four_point':   yields probably the best results
    stencil='eight_point':  probably too large

  If maxi<mask.size='auto' it might be that not all points of island are found.w

  If things take forever you can see with verbose=True how long it will still last.

  Usage: 
  itrue, jtrue, mask_out = get_island_mask(mask, iis, jjs)

last change:
----------
2017-08-16
  """
  if isinstance(maxit, str) and maxit=='auto':
    maxit = mask.size
  if not isinstance(stencil, str): 
    isten = stencil[0]
    jsten = stencil[1]
  elif stencil == 'eight_point':
    isten = np.array([-1,  0,  1, -1,  1, -1, 0,  1])
    jsten = np.array([-1, -1, -1,  0,  0,  1, 1,  1])
  elif stencil == 'four_point':
    isten = np.array([ 0, -1,  1, 0])
    jsten = np.array([-1,  0,  0, 1])
  else:
    print("Wrong choice for stencil:")
    print(stencil)
    sys.exit()

  # copy mask
  mask = 1*mask

  if mask[jjs, iis] != 0:
    itrue = np.array([iis], dtype=int)
    jtrue = np.array([jjs], dtype=int)
  else:
    itrue = np.array([], dtype=int)
    jtrue = np.array([], dtype=int)
  #  print('Started with wrong index!')
  #  sys.exit()
  
  cont = True
  nn = 0
  itocheck = np.array([], dtype=int)
  jtocheck = np.array([], dtype=int)
  while cont and nn<maxit:
    nn +=1 
    if verbose and nn%100==0:
      print('nn = % 4d / %d' % (nn, maxit))
  
    # check stencil around next point
    icheck = iis-isten
    jcheck = jjs-jsten
    mm = mask[jcheck, icheck]==1
  
    itocheck = np.concatenate((itocheck, icheck[mm]))
    jtocheck = np.concatenate((jtocheck, jcheck[mm]))
  
    itrue = np.concatenate((itrue, icheck[mm]))  
    jtrue = np.concatenate((jtrue, jcheck[mm]))  
  
    # mask found values to not count them double
    mask[jcheck[mm],icheck[mm]] = 0
  
    if itocheck.size == nn:
      cont = False
    else:
      iis = itocheck[nn] 
      jjs = jtocheck[nn] 
  mask_out = np.zeros(mask.shape)
  mask_out[jtrue, itrue] = 1
  return itrue, jtrue, mask_out

def haversine_dist(lonr, latr, lonps, latps, degree=True):
  """ Calculates distance of points on a sphere. 

  lonr and latr:       reference point (or point 1)
  lonps and latps:     other points
  degree:              if True coordinates are in degrees else in radians
  for details see http://en.wikipedia.org/wiki/Haversine_formula
  """
  r = 6378.e3
  if degree:
    lonr  = lonr  * np.pi/180.
    lonps = lonps * np.pi/180.
    latr  = latr  * np.pi/180.
    latps = latps * np.pi/180.
  arg = np.sqrt( np.sin(0.5*(latps-latr))**2 + np.cos(latr)*np.cos(latps)*np.sin(0.5*(lonps-lonr))**2 )
  dist = 2*r * np.arcsin(arg)
  return dist

def plt_stdcols(nn):
  return plt.rcParams['axes.prop_cycle'].by_key()['color'][nn]

def tbox(text, loc, ax, facecolor='w', alpha=1.0):
  bbox=dict(facecolor=facecolor, alpha=alpha, edgecolor='none')
  if loc=='ul':
    x = 0.03; y=0.95
    ha='left'; va='top'
  elif loc=='ur':
    x = 0.98; y=0.95
    ha='right'; va='top'
  elif loc=='ll':
    x = 0.03; y=0.05
    ha='left'; va='bottom'
  elif loc=='lr':
    x = 0.98; y=0.05
    ha='right'; va='bottom'
  ht = ax.text(x, y, text, ha=ha, va=va, bbox=bbox, transform=ax.transAxes)
  return

class nlcmap(LinearSegmentedColormap):
  """A nonlinear colormap

Taken from here:
https://stackoverflow.com/questions/33873397/nonlinear-colormap-with-matplotlib/33875488
  """
  name = 'nlcmap'
  def __init__(self, cmap, levels):
    self.cmap = cmap
    self.monochrome = self.cmap.monochrome
    self.levels = np.asarray(levels, dtype='float64')
    self._x = self.levels-self.levels.min()
    self._x/= self._x.max()
    self._y = np.linspace(0, 1, len(self.levels))
  def __call__(self, xi, alpha=1.0, **kw):
    yi = np.interp(xi, self._x, self._y)
    return self.cmap(yi, alpha)

# ---- start mpi stuff ----
class MPIobj(object):
  def __init__(self):
    try:
      from mpi4py import MPI
      self.comm = MPI.COMM_WORLD
      self.rank = self.comm.Get_rank()
      self.npro = self.comm.Get_size()
      self.use_mpi = True
    except:
      print("::: Warning, could not load mpi4py. Continue without mpi. :::")
      self.comm = 0
      self.rank = 0
      self.npro = 1
      self.use_mpi = False
    return

  def mpi_wait(self):
    npro = self.npro
    rank = self.rank
    comm = self.comm
    if self.use_mpi:
      for nn in range(npro):
        comm.send(0, dest=nn, tag=211)
      for nn in range(npro):
        wait = comm.recv(source=nn, tag=211)
    #time.sleep(0.5)
    return
  
  def mpi_print(self, string):
    npro = self.npro
    rank = self.rank
    comm = self.comm
    if self.use_mpi:
      for nn in range(npro):
        comm.send(string, dest=nn, tag=111)
        comm.send(rank,   dest=nn, tag=112)
      for nn in range(npro):
        string_r = comm.recv(source=nn, tag=111)
        rank_r   = comm.recv(source=nn, tag=112)
        if rank==0:
          print('proc %d: %s' % (rank_r, string_r))
    else:
      print('proc %d: %s' % (rank, string))
    return

  def get_steps(self, steps_all):
    npro = self.npro
    rank = self.rank
    steps_all = np.array(steps_all)

    ntasks = steps_all.size
    ntasks_per_node = (ntasks - ntasks%npro) / npro + 1
    task_ids = np.arange(rank*ntasks_per_node, (rank+1)*ntasks_per_node)
    task_ids = task_ids[task_ids<ntasks]
    steps = steps_all[task_ids]
    nsteps = steps.size

    self.mpi_print('nsteps = %d/%d'%(nsteps, steps_all.size))
    return steps, nsteps
# ---- end mpi stuff ----

# ---- start netcdf stuff ----
class Ncvar_old(object):
  """ Attribute container for a variable. """
  def __init__(self, data, name, grid_data, grid_names):
    self.data = data
    self.name = name
    self.grid_data = grid_data
    self.grid_names = grid_names
    return

class Ncvar(object):
  """ Attribute container for a variable. """
  def __init__(self, data, name, grid):
    self.data = data
    self.name = name
    self.grid = grid
    return

from netCDF4 import Dataset
class myDataset(Dataset):
  """ Main object to write netcdf file.
  
Usage: 
nc = Netcdf('/path/to/ncfile/ncfile.nc')
ncv = nc.add_var(data2d, [yt, xt])
ncv = nc.add_var(data3d, [zt, yt, xt])
nc.write_ncfile()
  """
  #def __init__(self, fpath, globs):
  #  self.fpath = fpath
  #  self.vlist = []
  #  self.vlist_gvars = []
  #  self.globs = globs#.copy()
  #  # --- open netcdf file
  #  self.fo = Dataset(self.fpath, 'w', format='NETCDF4')
  #  return

  #def getobjname(self, obj):
  #  """ Find the name of variable (works at least for numpy objects)."""
  #  for name in self.globs:
  #    if self.globs[name] is obj:
  #      outname = name
  #  return outname

  #def add_var_old(self, vdata, gdata, units=''):
  #  """ Add a variable to list so that it can be saved to netcdf file later."""
  #  # --- get variable name  and grid variable names as string
  #  vname = self.getobjname(vdata)
  #  gnames = [0]*len(gdata)
  #  for nn, gvar in enumerate(gdata):
  #    gnames[nn] = self.getobjname(gvar)
  #  # --- check dimensions
  #  gdims = []
  #  for gvar in gdata:
  #    gdims.append(gvar.size)
  #  if not list(vdata.shape) == gdims:
  #    raise ValueError('::: Error: %s has dimension %s and its grid has dimension %s! :::' % (vname, vdata.shape.__str__(), gdims.__str__()))
  #  # --- create Variable
  #  ncv = Ncvar(vdata, vname, gdata, gnames)
  #  # --- add attributes
  #  ncv.units = units
  #  # --- add variable to vlist
  #  self.vlist.append( ncv )
  #  return ncv

  def __init__(self, fpath):
    Dataset.__init__(self, fpath, 'w', format='NETCDF4')
    #super(Dataset, self).__init__(fpath, 'w', format='NETCDF4')
    self.var_list = np.array([])
    return

  def add_var(self, name, data, grid, units=''):
    # --- add variable
    ncv = self.createVariable(name, 'f4', grid)
    ncv.units = units
    ncv[:] = data 
    print(name)
    self.var_list = np.concatenate([self.var_list, [name]])
    print(type(self.var_list))
    print(type(self))
    return ncv

  def add_gvar(self, name, data, units=''):
    grid = [name]
    self.createDimension(name, data.size)
    self.add_var(name, data, grid, units=units)
    return

  def add_gvars(self, names, datas, units=''):
    for nn in range(len(names)):
      self.add_gvar(names[nn], datas[nn], units[nn])
    return

  #def write_ncfile(self):
  #  """ Finally, write all variables added with self.add_var to netcdf file."""
  #  for ncv in self.vlist:
  #    print('Processing %s...'%ncv.name)
  #    # --- add grid variables
  #    for nn, gname in enumerate(ncv.grid_names):
  #      if not gname in self.fo.variables.keys() :
  #        nc = self.fo.createDimension(gname, ncv.grid_data[nn].size)
  #        nc = self.fo.createVariable(gname, 'f4', (gname,))
  #        nc[:] = ncv.grid_data[nn]
  #    # --- add variable
  #    nc = self.fo.createVariable(ncv.name, 'f4', ncv.grid_names)
  #    nc.units = ncv.units
  #    nc[:] = ncv.data 
  #  self.fo.close()
  #  return
# ---- end netcdf stuff ----
