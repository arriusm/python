#!/Users/arri/miniconda3/bin/python
##!/afs/ipp/.cs/anaconda/amd64_generic/3/2019.10/bin/python
##!/project/ls-gruen/ligodesi/python/venv/bin/python
##!/opt/miniconda3/bin/python
##!/Users/arri/miniconda3/bin/python
##!/project/ls-gruen/ligodesi/miniconda3/bin/python

VERSION='showfits.py V1.4.11 (2024-12-05) (c) USM, Arno Riffeser (arri@usm.lmu.de) based on showfits of Claus Goessl and fitsedit of Johannes Koppenhoefer'

#import tkinter as tk
import matplotlib as mpl
mpl.use('TkAgg')   # ATTENTION: for TkAgg on MacOSX the button_press_event is not working correctly
#mpl.use('MacOSX')
#mpl.use('Gtk3Agg')

#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import  argparse 
import  numpy                                                              as np
import  copy
import  warnings
from    os.path                  import basename
from    scipy.stats              import sigmaclip
from    scipy.optimize           import leastsq, curve_fit, least_squares
#from   scipy.spatial            import Delaunay
from    astropy.io               import fits                               as pyfits
from    astropy.wcs              import WCS
from    astropy.utils.exceptions import AstropyWarning
import  matplotlib.pyplot                                                  as plt
from    matplotlib.widgets       import Cursor, Button
#from   matplotlib.backend_bases                                           as bb
#import matplotlib.tri                                                     as mtri
#import matplotlib.cm                                                      as cm      # python 3.7
from    matplotlib               import colormaps                          as cm      # python 3.10
from    pathlib   import  Path

warnings.filterwarnings("ignore", message='Toggling axes navigation from the keyboard is deprecated since 3.3 and will be removed two minor releases later.')
warnings.simplefilter('ignore', category=AstropyWarning)


# https://matplotlib.org/stable/api/prev_api_changes/api_changes_3.3.0.html#toggling-axes-navigation-from-the-keyboard-using-a-and-digit-keys
# 
# Toggling axes navigation from the keyboard using "a" and digit keys
# 
# Axes navigation can still be toggled programmatically using Axes.set_navigate.
# 
# The following related APIs are also deprecated: backend_tools.ToolEnableAllNavigation,
# backend_tools.ToolEnableNavigation, and rcParams["keymap.all_axes"].




class Image_Analyzer:

  def __init__(self, fitradius=7, precaper=80, maxaper=80, zoom=0, cmap=None, fixcuts=False, cuts=[0,0], wcslim=False,
               interpolation=None, verbose=2,
               fulloutput=False, outfile=None, autosavemasks=False, autosavepng=False, autosavejpg=False) :
    self.verbose = verbose
    self.fig = None
    self.fig2 = None
    self.cut_zoom = 1.2
    self.click2 = False
    self.sky = 0.
    self.linec = None
    self.pointc = None
    self.line1 = None
    self.line2 = None
    self.lineh = None
    self.totflux = 0.
    self.sky = 0.
    self.rad0 = 0.
    self.rad1 = 0. 
    self.fitradius = fitradius
    self.precaper = precaper
    self.maxaper = maxaper
    self.r = None
    self.grow = None
    self.growsky = None
    self.maxgrow = 0.
    self.area = None
    self.sky0 = 0.
    self.xc = 0
    self.yc = 0
    self.extent = []
    self.imshow = None
    self.title = ''
    self.imagename = ''
    self.outfile = outfile
    self.imamat = None
    self.header = None
    self.imamat_plot = None
    self.imamat0 = None
    self.ylim0 = 0.
    self.ylim1 = 0.
    self.imagelist = []
    self.imagelistend = 0
    self.imagenr = 0
    self.cmap = None
    self.cmap_nr = cmap
    self.cmap_all = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']
    #self.cmap_all = ['heat','gist_gray','gist_heat_r','gist_gray_r','gist_ncar','gist_ncar_r',
    #                 'gist_rainbow','gist_rainbow_r','gist_stern','gist_stern_r','gist_yarg',
    #                 'gist_yarg_r','gist_earth','gist_earth_r']
    self.cmap_listend = len(self.cmap_all)
    self.fixcuts = fixcuts
    self.cuts = cuts
    self.fulloutput = fulloutput
    self.zoom = zoom
    self.xlim = None
    self.ylim = None
    self.wcslim = wcslim
    self.wcs = None
    self.ra = None
    self.de = None
    self.ra_rad = None
    self.de_rad = None
    self.id1 = None
    self.id2 = None
    self.id3 = None
    self.idxlim = None
    self.idxlim = None
    self.idprev = None
    self.idnext = None
    self.idcut0n = None
    self.idcut0p = None
    self.idcut1n = None
    self.idcut1p = None
    self.idcmapn = None
    self.idcmapp = None
    self.id2_button_press = None
    self.id2_close  = None
    self.ax = None
    self.ax2 = None
    self.axprev = None
    self.axnext = None
    self.axcut0n = None
    self.axcut0p = None
    self.axcut1n = None
    self.axcut1p = None
    self.axcmapn  = None
    self.axcmapp  = None
    self.axmask  = None
    self.axmskradn = None
    self.axmskradp = None
    #self.key4active = False
    self.autosavemasks = autosavemasks
    self.autosavepng = autosavepng
    self.autosavejpg = autosavejpg
    self.mask = None
    self.integrad = 30
    self.mask_radius = 2
    self.satx0 = None
    self.satx1 = None
    self.saty0 = None
    self.saty1 = None
    self.sat0 = None
    self.sat1 = None
    self.interpolation = None
    
  def start(self, fig) :
    self.fig = fig
    # self.ax muss als 2. stehen: button_press_event BUG bei click+'2' !!
    self.ax        = self.fig.add_axes([0.07, 0.08, 0.91, 0.91]) 
    self.axprev    = self.fig.add_axes([0.01, 0.01, 0.10, 0.03])
    self.axnext    = self.fig.add_axes([0.12, 0.01, 0.10, 0.03])
    self.axcut0n   = self.fig.add_axes([0.25, 0.01, 0.06, 0.03])
    self.axcut0p   = self.fig.add_axes([0.32, 0.01, 0.06, 0.03])
    self.axcut1n   = self.fig.add_axes([0.39, 0.01, 0.06, 0.03])
    self.axcut1p   = self.fig.add_axes([0.46, 0.01, 0.06, 0.03])
    self.axcmapn   = self.fig.add_axes([0.55, 0.01, 0.06, 0.03])
    self.axcmapp   = self.fig.add_axes([0.62, 0.01, 0.06, 0.03])
    self.axmskradn = self.fig.add_axes([0.71, 0.01, 0.06, 0.03])
    self.axmskradp = self.fig.add_axes([0.78, 0.01, 0.06, 0.03])
    self.axmask    = self.fig.add_axes([0.85, 0.01, 0.06, 0.03])

    
  def storelist(self,imagelist) :
    self.imagelist = imagelist
    self.imagelistend = len(self.imagelist)
    print('nr_images  ',self.imagelistend)
    for i,s in enumerate(imagelist) :
      print('{:<3} {}'.format(i+1,s))


  def save_mask(self,event) :
    if (self.imagename!='') : 
      image_path = Path(self.imagename)
      mask_path = image_path.with_name('mask-'+image_path.name)
      print('newmask    ', mask_path)
      bla=self.imamat
      ### bla[np.where((np.isnan(self.imamat)))]=0   # keep 0s
      ### bla[np.where((bla!=0))]=1                  # keep 0s
      bla[np.where(np.invert(np.isnan(self.imamat)))]=1
      bla[np.where(          np.isnan(self.imamat) )]=0
      hdu = pyfits.PrimaryHDU(np.transpose(bla))
      hdu.scale('uint8')
      hdul = pyfits.HDUList([hdu])
      #hdul.writeto(maskname, overwrite=True)
      hdul.writeto(mask_path, overwrite=True)

  def load_mask(self,event) :
    if (self.imagename!='') : 
      image_path = Path(self.imagename)
      mask_path = image_path.with_name('mask-'+image_path.name)
      print('newmask    ', mask_path)
      maskbuf=self.imamat
      ### bla[np.where((np.isnan(self.imamat)))]=0   # keep 0s
      ### bla[np.where((bla!=0))]=1                  # keep 0s
      maskbuf[np.where(np.invert(np.isnan(self.imamat)))]=1
      maskbuf[np.where(          np.isnan(self.imamat) )]=0

      hdu = pyfits.PrimaryHDU(np.transpose(bla))
      hdu.scale('uint8')
      hdul = pyfits.HDUList([hdu])
      #hdul.writeto(maskname, overwrite=True)
      hdul.writeto(mask_path, overwrite=True)

      
  def save_png(self,event) :
    if (self.imagename!='') : 
      image_path = Path(self.imagename)
      png_path = image_path.with_stem('png')
      print('png        ', png_path)
      self.fig.savefig(png_path)

  def save_jpg(self,event) :
    if (self.imagename!='') : 
      image_path = Path(self.imagename)
      jpg_path = image_path.with_stem('jpg')
      print('jpg        ', jpg_path)
      self.fig.savefig(jpg_path)

  def load(self,cut0=None,cut1=None) :
    if (self.autosavemasks) :
      self.save_mask(self)
    if (self.autosavepng) :
      self.save_png(self)
    if (self.autosavejpg) :
      self.save_jpg(self)
    if self.verbose>=2 :
        print('=============================================================================================================')
        
    self.imagename = self.imagelist[self.imagenr]
    self.title = '{}     {} of {}'.format(self.imagename,self.imagenr+1,self.imagelistend )
    self.header      = pyfits.getheader(self.imagename)
    self.imamat_plot = pyfits.getdata(self.imagename).astype(np.float32)
    self.imamat      = np.transpose(self.imamat_plot)
    (self.nx,self.ny) = np.shape(self.imamat)
    # small = self.imamat[0:int(self.nx/3),0:int(self.ny/3)]
    small = self.imamat[int(self.nx/3):int(2*self.nx/3),int(self.ny/3):int(2*self.ny/3)]
    # ima_med = np.nanmedian(self.imamat)
    # ima_med = np.nanmedian(small)
    # ima_std = np.nanstd(small) # NAN are more problematic....
    # ima_std = ima_med - np.nanpercentile(small,15.865)
    # #print('ima_std=',ima_std)
    # if -1e-8<ima_std<1e-8 :
    #   ima_std=(np.max(self.imamat)-np.min(self.imamat))/4.
    #ima_min = np.nanmin(small)
    #ima_max = np.nanmax(small)
    #ima_del = ima_max - ima_min
    if (cut0==None) :
      #cut0 = ima_min - 0.1 * ima_del
      #cut0 = ima_med - 1.5 * ima_std
      # good for flat images: cut0 = ima_med - 1.0 * ima_std 
      # bis 24-03-09  cut0 = ima_med - 5.0 * ima_std
      if self.cuts[0]==0 and self.cuts[1]==0 :
        cut0 = np.nanpercentile(small,1.)
      else :
        cut0 = self.cuts[0]
    if (cut1==None) :
      #cut1 = ima_max - 0.1 * ima_del
      #cut1 = ima_med + 4.0 * ima_std
      # good for flat images: cut1 = ima_med + 3.0 * ima_std
       # bis 24-03-09  cut1 = ima_med + 5.0 * ima_std
      if self.cuts[0]==0 and self.cuts[1]==0 :
        cut1 = np.nanpercentile(small,99.)*2.
      else :
        cut1 = self.cuts[1]
    if verbose>=3 :
      print('np.shape(self.imamat)=',np.shape(self.imamat))
    (self.nx,self.ny) = np.shape(self.imamat)
    self.cut0 = cut0
    self.cut1 = cut1
    self.dcut = cut1 - cut0
    if verbose>=1 :
      print('nr         ', self.imagenr+1,' of ',self.imagelistend)
      print('imagename  ', self.imagename)
      #print('exptime     {:<10.2f}'.format(self.header['EXPTIME']))
    if self.verbose>=2 :
      print('imagesize   {:<10.0f} {:<10.0f}'.format(self.nx,self.ny))
      print('cuts        {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    if verbose>=3 :
      print('#')
    self.extent = [0.5,self.nx+0.5,0.5,self.ny+0.5]
    # print('self.ax=',self.ax)
    self.ax.clear()

    self.cmap = copy.copy(cm.get_cmap(self.cmap_all[self.cmap_nr]))
    self.imshow = self.ax.imshow( self.imamat_plot, origin='lower',
                                  cmap=self.cmap, vmin=self.cut0,
                                  vmax=self.cut1, extent=self.extent, interpolation=self.interpolation)
    # self.fig.canvas.set_window_title(self.title)
    # self.imshow.figure.canvas.set_window_title(self.title)
    self.imshow.figure.canvas.manager.set_window_title(self.title)
    current_cmap = self.imshow.get_cmap()
    current_cmap.set_bad(color='blue')

    if self.wcslim :      
      self.wcs = WCS(self.header,relax=True,fix=False)
      if (self.ra!=None and self.de!=None and self.ra_rad!=None and self.de_rad!=None) :
        px = np.zeros(4)
        py = np.zeros(4)
        px[0],py[0] = self.wcs.world_to_pixel_values(self.ra-self.ra_rad, self.de-self.de_rad) + np.array([1.,1.])
        px[1],py[1] = self.wcs.world_to_pixel_values(self.ra+self.ra_rad, self.de-self.de_rad) + np.array([1.,1.])
        px[2],py[2] = self.wcs.world_to_pixel_values(self.ra+self.ra_rad, self.de+self.de_rad) + np.array([1.,1.])
        px[3],py[3] = self.wcs.world_to_pixel_values(self.ra-self.ra_rad, self.de+self.de_rad) + np.array([1.,1.])
        # dx = self.xlim[1]-self.xlim[0]
        # dy = self.ylim[1]-self.ylim[0]
        # self.xlim = (px-dx/2.,px+dx/2.)
        # self.ylim = (py-dy/2.,py+dy/2.)
        self.xlim = ( np.min(px), np.max(px) )
        self.ylim = ( np.min(py), np.max(py) )
      else :
        self.xlim = self.ax.get_xlim()
        self.ylim = self.ax.get_ylim()
        xcent = (self.xlim[0]+self.xlim[1])/2.
        ycent = (self.ylim[0]+self.ylim[1])/2.
        self.ra, self.de = self.wcs.pixel_to_world_values(xcent-1,ycent-1) # -1 because of an astropy python array inconsistency
        ra = np.zeros(4)
        de = np.zeros(4)
        ra[0], de[0] = self.wcs.pixel_to_world_values(self.xlim[0]-1,ycent-1)  # -1 because of an astropy python array inconsistency
        ra[1], de[1] = self.wcs.pixel_to_world_values(self.xlim[1]-1,ycent-1)  # -1 because of an astropy python array inconsistency
        ra[2], de[2] = self.wcs.pixel_to_world_values(xcent-1,self.ylim[0]-1)  # -1 because of an astropy python array inconsistency
        ra[3], de[3] = self.wcs.pixel_to_world_values(xcent-1,self.ylim[1]-1)  # -1 because of an astropy python array inconsistency
        self.ra_rad = np.abs( ra[1] - ra[0] ) / 2.
        self.de_rad = np.abs( de[3] - de[2] ) / 2.
      print('radec    {:10.4f} {:10.4f}'.format(self.ra, self.de))        
        
    if (self.xlim==None and self.ylim==None) :
      if self.zoom!=0 :
        self.ax.xaxis.zoom(self.zoom) 
        self.ax.yaxis.zoom(self.zoom)
    else :
      self.ax.set_xlim(self.xlim)
      self.ax.set_ylim(self.ylim)

    self.idxlim = self.ax.callbacks.connect('xlim_changed', self.on_xlims_change)
    self.idylim = self.ax.callbacks.connect('ylim_changed', self.on_ylims_change)

    # print('lim=',self.ax.get_xlim(),self.ax.get_xlim())

    # if self.zoom!=0 :
    #   if self.zoom==4 :
    #     dx = 100
    #     dy = 100
    #   if self.zoom==2 :
    #     dx = 200
    #     dy = 200
    #   if self.zoom==1 :
    #     dx = 400
    #     dy = 400
    #   if self.zoom==0.5 :
    #     dx = 800
    #     dy = 800
    #   self.ax.set_xlim([self.nx/2-dx,self.nx/2+dx])
    #   self.ax.set_ylim([self.ny/2-dy,self.ny/2+dy])
    self.fig.canvas.draw_idle()

    
  def connect(self):
    #print(self.imshow)
    self.id1 = self.fig.canvas.mpl_connect('button_press_event',self.onclick)
    # self.id1 = self.ax.figure.canvas.mpl_connect('button_press_event',self.onclick)
    # self.id1 = self.ax.on_clicked(self.onclick)
    self.id2 = self.fig.canvas.mpl_connect('key_press_event',self.keypress)
    # self.id3 = self.fig.canvas.mpl_connect('close_event',self.close)
    # cursor = Cursor(self.ax, useblit=False, color='g', linewidth=1 )
    self.idprev  = Button(self.axprev, 'previous')
    self.idnext  = Button(self.axnext, 'next')
    self.idcut0n = Button(self.axcut0n, 'Lcut -')
    self.idcut0p = Button(self.axcut0p, 'Lcut +')
    self.idcut1n = Button(self.axcut1n, 'Hcut -')
    self.idcut1p = Button(self.axcut1p, 'Hcut +')
    self.idcmapn = Button(self.axcmapn, 'cmap -')
    self.idcmapp = Button(self.axcmapp, 'cmap +')
    self.idmask = Button(self.axmask, 'mask')
    self.idmskradn = Button(self.axmskradn, 'rad -')
    self.idmskradp = Button(self.axmskradp, 'rad +')
    self.idprev.label.set_fontsize(8)
    self.idnext.label.set_fontsize(8)
    self.idcut0n.label.set_fontsize(8)
    self.idcut0p.label.set_fontsize(8)
    self.idcut1n.label.set_fontsize(8)
    self.idcut1p.label.set_fontsize(8)
    self.idcmapn.label.set_fontsize(8)
    self.idcmapp.label.set_fontsize(8)
    self.idmask.label.set_fontsize(8)
    self.idmskradn.label.set_fontsize(8)
    self.idmskradp.label.set_fontsize(8)
    self.idprev.on_clicked(self.prev)
    self.idnext.on_clicked(self.next)
    self.idcut0n.on_clicked(self.cut0n)
    self.idcut0p.on_clicked(self.cut0p)
    self.idcut1n.on_clicked(self.cut1n)
    self.idcut1p.on_clicked(self.cut1p)
    self.idcmapn.on_clicked(self.change_cmap_neg)
    self.idcmapp.on_clicked(self.change_cmap_plus)
    self.idmask.on_clicked(self.save_mask)
    self.idmskradn.on_clicked(self.set_mask_radius_minus)
    self.idmskradp.on_clicked(self.set_mask_radius_plus)

          
  def disconnect(self) :
    if (self.autosavemasks) :
      self.save_mask(self)
    if (self.autosavepng) :
      self.save_png(self)
    if (self.autosavejpg) :
      self.save_jpg(self)
    self.fig.canvas.mpl_disconnect(self.id1)
    self.fig.canvas.mpl_disconnect(self.id2)
    self.fig.canvas.mpl_disconnect(self.id3)
    self.ax.callbacks.disconnect(self.idxlim)
    self.ax.callbacks.disconnect(self.idylim)

      
  def close(self, event) :
    input('really close? <enter>')

  def set_mask_radius_minus(self, event) :
    self.mask_radius -= 1
    if self.mask_radius<0 :
      self.mask_radius = 0
    print('mask_radius {:<10.0f}'.format(self.mask_radius))

  def set_mask_radius_plus(self, event) :
    self.mask_radius += 1
    print('mask_radius {:<10.0f}'.format(self.mask_radius))

    
  def cut0n(self, event) :
    self.dcut = self.cut1 - self.cut0
    self.cut0 = self.cut1 - self.dcut*self.cut_zoom     
    print('cuts        {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    #self.imshow = self.ax.imshow( self.imamat_plot, origin='lower', cmap=self.cmap, vmin=self.cut0,
    #                              vmax=self.cut1, extent=self.extent)
    self.imshow.set_clim(vmin=self.cut0)
    self.fig.canvas.draw_idle()

    
  def cut0p(self, event) :
    self.dcut = self.cut1 - self.cut0
    self.cut0 = self.cut1 - self.dcut/self.cut_zoom 
    print('cuts        {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    #self.imshow = self.ax.imshow( self.imamat_plot, origin='lower', cmap=self.cmap, vmin=self.cut0,
    #                              vmax=self.cut1, extent=self.extent)
    self.imshow.set_clim(vmin=self.cut0)
    self.fig.canvas.draw_idle()

    
  def cut1n(self, event) :
    self.dcut = self.cut1 - self.cut0
    self.cut1 = self.cut0 + self.dcut/self.cut_zoom   
    print('cuts        {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    #self.imshow = self.ax.imshow( self.imamat_plot, origin='lower', cmap=self.cmap, vmin=self.cut0,
    #                              vmax=self.cut1, extent=self.extent)
    self.imshow.set_clim(vmax=self.cut1)
    self.fig.canvas.draw_idle()

    
  def cut1p(self, event):
    self.dcut = self.cut1 - self.cut0
    self.cut1 = self.cut0 + self.dcut*self.cut_zoom 
    print('cuts        {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    #self.imshow = self.ax.imshow( self.imamat_plot, origin='lower', cmap=self.cmap, vmin=self.cut0,
    #                              vmax=self.cut1, extent=self.extent)
    self.imshow.set_clim(vmax=self.cut1)
    self.fig.canvas.draw_idle()

    
  def change_cmap_plus(self, event):
    if (self.cmap_nr<self.cmap_listend-1) :
      self.cmap_nr = self.cmap_nr + 1
    else :
      self.cmap_nr = 0
    print('cmap        {:<10} {:<}'.format(self.cmap_nr,self.cmap_all[self.cmap_nr]))
    self.cmap = copy.copy(cm.get_cmap(self.cmap_all[self.cmap_nr]))   # python 3.7
    #self.imshow = self.ax.imshow( self.imamat_plot, origin='lower', cmap=self.cmap, vmin=self.cut0,
    #                              vmax=self.cut1, extent=self.extent)
    self.imshow.set_cmap(self.cmap)
    current_cmap = self.imshow.get_cmap()
    current_cmap.set_bad(color='blue')
    self.fig.canvas.draw_idle()    

  def change_cmap_neg(self, event):
    if (self.cmap_nr>0) :
      self.cmap_nr = self.cmap_nr - 1
    else :
      self.cmap_nr = self.cmap_listend-1
    print('cmap        {:<10} {:<}'.format(self.cmap_nr,self.cmap_all[self.cmap_nr]))
    self.cmap = copy.copy(       cm.get_cmap(self.cmap_all[self.cmap_nr]))   # python 3.7
    #self.imshow = self.ax.imshow( self.imamat_plot, origin='lower', cmap=self.cmap, vmin=self.cut0,
    #                              vmax=self.cut1, extent=self.extent)
    self.imshow.set_cmap(self.cmap)
    current_cmap = self.imshow.get_cmap()
    current_cmap.set_bad(color='blue')
    self.fig.canvas.draw_idle()    

    
  def test(self, nix=None) :
    print('id2 = ',self.ax)
    print('id2 = ',f'{id(self.ax)}')
    return 0


  def get_cursor(self, pos, verbose=0) :
    xc = int(pos[0]+0.5)
    yc = int(pos[1]+0.5)
    if self.verbose>=2 :
        print('{:<10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'value'))
    if self.verbose>=1 :         
        print('{:10} {:10.0f} {:10.0f} {:10.2f}'
              .format('pixel',xc,yc,self.imamat[xc-1,yc-1]))
    if self.fulloutput :
      self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10}\n'
                         .format('#','method','xc', 'yc', 'value'))
      self.outfile.write('{:<30} {:<10} {:10.0f} {:10.0f} {:10.2f}\n'
                         .format(self.imagename,'pixel',xc,yc,self.imamat[xc-1,yc-1]))
    else :
      self.outfile.write('{:<30} {:10.0f} {:10.0f} {:10.2f}\n'
                         .format(self.imagename,xc,yc,self.imamat[xc-1,yc-1]))
    #self.ax.plot(xc,yc,marker='2',c='k',ms=10)
    return xc,yc

  
  def gauss(self, xx, yy, xx0, yy0, sigx, sigy, phi, A, C=0. ):
    cw=np.cos(phi);
    sw=np.sin(phi);
    tx= (xx-xx0)*cw+(yy-yy0)*sw;
    ty=-(xx-xx0)*sw+(yy-yy0)*cw;
    # double ert = -( a(3) * tx * tx + a(4) * ty * ty );
    ert = 0.5*(tx/sigx)**2 +0.5*(ty/sigy)**2
    return C + A * np.exp( -ert )
  
    
  def residuals_gauss(self, a, XY, Z):
    yy, xx = XY
    #xx0, yy0, sigx, sigy, phi, A = a
    #C = 0.
    xx0, yy0, sigx, sigy, phi, A, C = a 
    return self.gauss(xx, yy, xx0, yy0, sigx, sigy, phi, A, C) - Z   
  
  
  def moffat(self, xx, yy, xx0, yy0, sigx, sigy, phi, A, C ):
    cw=np.cos(phi);
    sw=np.sin(phi);
    tx= (xx-xx0)*cw+(yy-yy0)*sw;
    ty=-(xx-xx0)*sw+(yy-yy0)*cw;
    # double ert = 1. + ( a(3) * tx ) * ( a(3) * tx ) + ( a(4) * ty ) * ( a(4) * ty );
    ert = 1.+0.5*(tx/sigx)**2 +0.5*(ty/sigy)**2
    return C + A * ert**(-3.0)
  
    
  def residuals_moffat(self, a, XY, Z):
    yy, xx = XY
    #xx0, yy0, sigx, sigy, phi, A = a
    #C = 0.
    xx0, yy0, sigx, sigy, phi, A, C = a 
    return self.moffat(xx, yy, xx0, yy0, sigx, sigy, phi, A, C) - Z


  def show_psf(self, pos, verbose=0) :
    xc = int(pos[0]+0.5)
    yc = int(pos[1]+0.5)
    rad = self.fitradius
    boxsize = 2*rad+1
    #print('x,y=',pos[0],pos[1])
    #AR x0 = xc-rad
    #AR x1 = xc+rad
    #AR y0 = yc-rad
    #AR y1 = yc+rad
    #AR DATA = self.imamat[x0-1:x1,y0-1:y1]    # minus 1: fits array -> python array
    #AR max = np.nanmax(DATA)
    #AR cen = DATA[rad,rad]
    #AR x = np.linspace(x0, x1, boxsize, dtype=int)
    #AR y = np.linspace(y0, y1, boxsize, dtype=int)
    #AR #XV, YV = np.meshgrid(y,x)
    #AR for ix in range(0,boxsize) :
    #AR     for iy in range(0,boxsize) :
    #AR         if DATA[ix,iy]==max :
    #AR             xc = x[ix]
    #AR             yc = y[iy]
    #AR print('x,y=',xc,yc)
    x0 = xc-rad
    x1 = xc+rad
    y0 = yc-rad
    y1 = yc+rad
    DATA = self.imamat[x0-1:x1,y0-1:y1]    # minus 1: fits array -> python array
    x = np.linspace(x0, x1, boxsize, dtype=float)
    y = np.linspace(y0, y1, boxsize, dtype=float)
    lensubs = 101
    subs=np.linspace(-5.0, 5.0, lensubs, dtype=float)
    #ARprint('subs=',subs)
    r = np.zeros((boxsize**2,lensubs**2),dtype=float)
    p = np.zeros((boxsize**2,lensubs**2),dtype=float)
    F = np.zeros(boxsize**2,dtype=float)
    ixy=0
    for ix in range(0,boxsize) :
      for iy in range(0,boxsize) :
        k=0
        #p[ixy,k]=np.arctan2(y[iy]-pos[1],x[ix]-pos[0])/np.pi*180.
        p[ixy,k]=np.arctan((y[iy]-pos[1])/(x[ix]-pos[0]))/np.pi*180.
        r[ixy,k]=np.sqrt( (x[ix]-pos[0])**2 + (y[iy]-pos[1])**2 )
        # for i,xsub in enumerate(subs) :
        #   for j,ysub in enumerate(subs) :
        #     r[ixy,k]=np.sqrt( (x[ix]-xc-xsub)**2 + (y[iy]-yc-ysub)**2 )
        #     k+=1
        F[ixy]=DATA[ix,iy]
        ixy+=1
    kbest=0
    #AR #print('indx=',indx)
    #AR scat = np.zeros(lensubs**2)
    #AR for k in range(lensubs**2) :
    #AR   indx = np.argsort(r[:,k])
    #AR   #rnew = r[i,indx]
    #AR   Fnew = F[indx]
    #AR   for z in range(1,boxsize**2) :
    #AR     scat[k] += abs(Fnew[z]-Fnew[z-1])
    #AR print(scat)
    #AR kbest = np.argmin(scat)
    #AR print(kbest)
    #AR k=0
    #AR for i,xsub in enumerate(subs) :
    #AR   for j,ysub in enumerate(subs) :
    #AR     if k==kbest :
    #AR       print('center=',xc+xsub,yc+ysub)
    #AR     k+=1
    indx = np.argsort(r[:,kbest])
    rnew = r[indx,kbest]
    pnew = p[indx,kbest]
    Fnew = F[indx]
    #print('r=',r)
    #print('F=',F)
    return rnew,pnew,Fnew
                
  
  def fit_psf(self, pos, use_moffat=False, verbose=0) :
    xc = int(pos[0]+0.5)
    yc = int(pos[1]+0.5)
    rad = self.fitradius
    boxsize = 2*rad+1
    x0 = xc-rad
    x1 = xc+rad
    y0 = yc-rad
    y1 = yc+rad
    sky = np.nanmedian(self.imamat)
    # print('sky=',sky)
    # DATA = self.imamat[x0-1:x1,y0-1:y1] - sky  # minus 1 in den indizes
    #                                              wegen image array -> python array
    DATA = self.imamat[x0-1:x1,y0-1:y1]    # minus 1 in den indizes wegen image
    #                                        array -> python array
    #bg  = np.nanmedian(DATA)
    bg  = np.nanpercentile(DATA,0.1)
    max = np.nanmax(DATA)
    amp = max-bg
    cen = DATA[rad,rad]
    # x, y = np.linspace(x0, x1, boxsize, dtype=float), np.linspace(y0, y1,boxsize, dtype=float)
    x = np.linspace(x0, x1, boxsize, dtype=float)
    y = np.linspace(y0, y1, boxsize, dtype=float)
    # print('x,y=',x,y)
    XV, YV = np.meshgrid(y,x)
    # for ix in range(0,boxsize) :
    #     for iy in range(0,boxsize) :
    #         if DATA[ix,iy]==MD :
    #             xc = x[ix]
    #             yc = y[iy]

    ix,iy = np.unravel_index(np.argmax(DATA, axis=None), DATA.shape)
    X_in = XV.ravel()
    Y_in = YV.ravel()
    Z_in = DATA.ravel()
    bad = np.where( np.isnan(Z_in) )
    X = np.delete(X_in,bad)
    Y = np.delete(Y_in,bad)
    Z = np.delete(Z_in,bad)
    XY = np.vstack((X,Y))

    afit=np.array([np.nan,np.nan])
    if not use_moffat : ########### GAUSS ###########
      #astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp )
      #astart = ( xc , yc , 4.1 , 4. , 45./180.*np.pi , amp , bg )  #GOOD
      #astart = ( xc , yc , 4.5 , 4.0 , 60./180.*np.pi , amp , bg )  # AR231215
      astart = ( x0+ix , y0+iy , 3.5 , 4.0 , 60./180.*np.pi , amp , bg )  # AR231215
      res = least_squares(self.residuals_gauss, astart, method='lm', args=(XY,Z), verbose=0)
      #afit = np.zeros(7)
      #afit[0:6] = res.x
      #afit[6] = 0.
      afit = res.x
      if ( not np.isnan(afit[0]) and not np.isnan(afit[1]) and
           afit[0]>x0 and afit[0]<x1 and afit[1]>y0 and afit[1]<y1 and afit[5]>5.*abs(np.sqrt(afit[6]))) :
        # double ert = -( a(3) * tx * tx + a(4) * ty * ty );  # usmlib.h
        # fwhm0  = 2. * sqrt( M_LN2 / a(3) );                 # starphot.cpp
        # refa(0)  =  refa(2)*M_PI/sqrt(refa(3)*refa(4));     # starphot.cpp
        fwhmx = abs(afit[2]) * np.sqrt(8.*np.log(2.))
        fwhmy = abs(afit[3]) * np.sqrt(8.*np.log(2.))
        fwhm = np.sqrt(fwhmx*fwhmy)
        a3 = 1./(2.*afit[2]**2)
        a4 = 1./(2.*afit[3]**2)
        totflux = afit[5] * np.pi / np.sqrt( a3 * a4);
        if verbose>=5 :
          xn = np.linspace(-self.integrad, self.integrad, 2*self.integrad+1, dtype=float)
          yn = np.linspace(-self.integrad, self.integrad, 2*self.integrad+1, dtype=float)
          Xn, Yn = np.meshgrid(yn,xn)
          yyn = Xn.ravel()
          xxn = Yn.ravel()
          xx0, yy0, sigx, sigy, phi, A, C = afit
          print('numerical totflux = ',
                np.sum(self.gauss(xxn, yyn, 0., 0., sigx, sigy, phi, A, 0.)))        
        if verbose>=2 :
          print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'
                .format('#','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
        if verbose>=1 :
          print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}'
                .format('gauss',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))

        if verbose>=0 :
          if self.fulloutput :
            self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'
                               .format('#','method','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
            self.outfile.write('{:<30} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'
                               .format(self.imagename,'gauss',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),
                                       abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
          else :
            self.outfile.write('{:<30} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'
                             .format(self.imagename,afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),
                                     abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))            
        #self.ax.plot(afit[0],afit[1],marker='+',c='g',ms=10)
      else :
        if np.isnan(afit[0]) or not afit[0]>x0 or not afit[0]<x1 :
          print("   fit error: xc=",xc,"-->",xc+afit[0])
        if np.isnan(afit[1]) or not afit[1]>x0 or not afit[1]<x1 :
          print("   fit error: yc=",yc,"-->",yc+afit[1])
        if afit[5]<0.1*afit[6] :
          print("   fit error:  amp (too small)",afit[5],5.*np.sqrt(abs(afit[6])),afit[6])
        afit[0]=np.nan
        afit[1]=np.nan
    else : #################### MOFFAT ####################
      #astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp )
      astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp , bg )
      res = least_squares(self.residuals_moffat, astart, method='lm', args=(XY,Z), verbose=0) 
      #afit = np.zeros(7)
      #afit[0:6] = res.x
      #afit[6] = 0.
      afit = res.x
      if (not np.isnan(afit[0]) and not np.isnan(afit[1]) and
          afit[0]>x0 and afit[0]<x1 and afit[1]>y0 and afit[1]<y1 and afit[5]>5.*np.sqrt(abs(afit[6]))) :
        # ert = 1.+0.5*(tx/sigx)**2 +0.5*(ty/sigy)**2
        # double ert = 1. + ( a(3) * tx ) * ( a(3) * tx ) + ( a(4) * ty ) * ( a(4) * ty );   # usmlib.h
        # fwhm0 = 2. * sqrt( pow( 2., 1. / a(8) ) - 1. ) / a(3);                             # starphot.cpp
        # a0(0)  =  a(2) * M_PI / ( (a(8)-1.) * a(3) * a(4) );                               # starphot.cpp
        a8 = 3.0
        a3 = np.sqrt( 1./(2.*afit[2]**2) )
        a4 = np.sqrt( 1./(2.*afit[3]**2) )
        fwhmx = 2. * np.sqrt( 2.**( 1. / a8 ) - 1. ) / a3
        fwhmy = 2. * np.sqrt( 2.**( 1. / a8 ) - 1. ) / a4
        fwhm = np.sqrt(fwhmx*fwhmy)
        totflux = afit[5] * np.pi / ( (a8-1.) * a3 * a4 )
        if verbose>=5 :
          print('analytical totflux = ',totflux)
          xn = np.linspace(-self.integrad, self.integrad, 2*self.integrad+1, dtype=float)
          yn = np.linspace(-self.integrad, self.integrad, 2*self.integrad+1, dtype=float)
          Xn, Yn = np.meshgrid(yn,xn)
          yyn = Xn.ravel()
          xxn = Yn.ravel()
          xx0, yy0, sigx, sigy, phi, A, C = afit
          print('numerical totflux = ',np.sum(self.moffat(xxn, yyn, 0., 0., sigx, sigy, phi, A, 0.)))        
        if verbose>=2 :
          print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'
                .format('#','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
        if verbose>=1 :
          print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}'
                .format('moffat',afit[0],afit[1],totflux,afit[6],afit[5],
                        abs(afit[2]),abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
        if self.fulloutput :
          if verbose>=0 :
            self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'
                               .format('#','method','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
            self.outfile.write('{:<30} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'
                               .format(self.imagename,'moffat',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),
                                       abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
        else :
            self.outfile.write('{:<30} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'
                               .format(self.imagename,afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),
                                       abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
            
          #self.ax.plot(afit[0],afit[1],marker='x',c='b',ms=10)
      else :
        if np.isnan(afit[0]) or not afit[0]>x0 or not afit[0]<x1 :
          print("   fit error: xc=",xc,"-->",xc+afit[0])
        if np.isnan(afit[1]) or not afit[1]>x0 or not afit[1]<x1 :
          print("   fit error: yc=",yc,"-->",yc+afit[1])
        if afit[5]<0.1*afit[6] :
          print("   fit error:  amp (too small)",afit[5],5.*np.sqrt(abs(afit[6])),afit[6])
        afit[0]=np.nan
        afit[1]=np.nan

    #if np.isnan(afit[0]) and np.isnan(afit[1]) :
    #  print('wrong click?')

    return afit[0],afit[1]


  def growth_curve(self, event) :
    if not self.click2 and event.xdata!=None and event.ydata!=None :
      r1 = event.xdata
      self.line1.set_xdata([r1,r1])
      self.line1.set_ydata([self.ylim0,self.ylim1])
      self.line1.figure.canvas.draw_idle()
      self.rad0 = r1
      self.rad1 = r1
      self.click2 = True
    elif event.xdata!=None and event.ydata!=None :
      r2 = event.xdata
      #print(r2,self.rad0 ,self.rad1)
      if   r2<self.rad0 :
        self.rad0 = r2
      elif r2>self.rad1 :
        self.rad1 = r2
      elif (r2-self.rad0)<(self.rad1-r2) :
        self.rad0 = r2
      else :
        self.rad1 = r2
      # extrapolate to sky by fitting parabola to intervall
      notfit = np.where( (self.r<self.rad0) | (self.rad1<self.r) )
      rfit = np.delete(self.r, notfit)
      areafit = np.delete(self.area, notfit)
      growfit = np.delete(self.grow, notfit)
      growskyfit = np.delete(self.growsky, notfit)
      self.rad0 = rfit[0]
      self.rad1 = rfit[-1]
      # https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
      # fitlist =  opt.curve_fit(sky, rfit,growfit, x0, sigma)
      # A = np.vstack([rfit**2*np.pi,np.ones(len(rfit))]).T
      A = np.vstack([areafit,np.ones(len(rfit))]).T
      lstsqfit = np.linalg.lstsq(A, growfit,rcond=None)
      #print('lstsqfit=',lstsqfit)
      a2,a0 = lstsqfit[0]
      rx = np.linspace(0.,100,100)
      self.linef.set_xdata(rx)
      self.linef.set_ydata(a2*rx**2*np.pi+a0)
      self.linef.figure.canvas.draw_idle()
      self.totflux = a0
      #print(growfit[0]-a2*areafit[0])
      self.sky = self.sky0 + a2
      new_grow = self.grow - a2*self.area
      new_growfit = np.delete(new_grow , notfit)
      totall = growskyfit[0]
      errtotall = np.sqrt(totall)
      #totall = self.totflux+self.sky*self.rad0**2*np.pi
      if self.verbose>=3 :
        print('rad0 =',self.rad0)
        print('rad1 =',self.rad1)
        print('totflux        = ',self.totflux)
        print('new_growfit[0] = ',new_growfit[0])
        print('totflux+sky*rad0**2*np.pi = ',self.totflux+self.sky*self.rad0**2*np.pi)
        #print('totflux+sky*areafit[0]    = ',new_growfit[0]+self.sky*areafit[0])
        #print('totflux+sky*pi*r0^2       = ',new_growfit[0]+self.sky*areafit[0])
        print('growskyfit[0]             = ',growskyfit[0])
        #print('totflux+sky*rfit[0]**2*pi = ',self.totflux+self.sky*rfit[0]**2*np.pi)
      self.linec.set_ydata(new_grow)
      self.pointc.set_ydata(new_grow)
      self.lineh.set_xdata([self.r[0],self.r[-1]])
      self.lineh.set_ydata([a0,a0])
      new_grow_min = np.nanmin(new_grow)
      new_grow_max = np.nanmax(new_grow)
      new_grow_del = new_grow_max - new_grow_min
      # replot
      self.ylim0 = new_grow_min-0.1*new_grow_del
      self.ylim1 = new_grow_max+0.1*new_grow_del
      self.line1.set_xdata([self.rad0,self.rad0])
      self.line1.set_ydata([self.ylim0,self.ylim1])
      self.line2.set_xdata([self.rad1,self.rad1])
      self.line2.set_ydata([self.ylim0,self.ylim1])
      self.line2.figure.canvas.draw_idle()
      self.ax2.set_ylim([self.ylim0,self.ylim1])
      ax2title = '{}  {:.1f}  {}'.format(self.imagename,self.totflux,'(growth curve)')
      self.ax2.set_title(ax2title)
      self.lineh.figure.canvas.draw_idle()
      #self.click2 = False
      if (event.dblclick) :
        print('---------------------------------------------------------------------------------------')
      if self.verbose>=1 :         
        print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f} {:10.1f}'
              .format('grow_curve',self.xc,self.yc,self.totflux,self.sky,self.rad0,
                      self.rad1,errtotall))
      if (event.dblclick) :
        print('---------------------------------------------------------------------------------------')
        if self.fulloutput :
          self.outfile.write('{:<30} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f} {:10.1f}\n'
                             .format(self.imagename,'grow_curve',self.xc,self.yc,
                                     self.totflux,self.sky,self.rad0,self.rad1,errtotall))
        else :
          self.outfile.write('{:<30} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f} {:10.1f}\n'
                             .format(self.imagename,self.xc,self.yc,self.totflux,
                                     self.sky,self.rad0,self.rad1,errtotall))

          
  def growth_curve_close(self, event_buf) :
    print('---------------------------------------------------------------------------------------')
    self.fig2.canvas.mpl_disconnect(self.id2_button_press)
    self.fig2.canvas.mpl_disconnect(self.id2_close)
    self.ax2.set_navigate(False)
    self.ax.set_navigate(True)


  def on_xlims_change(self,event):
    self.xlim = event.get_xlim()
    if self.verbose>=3 :
      print('xlimits     {:<10.3f} {:<10.3f}'.format(self.xlim[0],self.xlim[1]))
    if (self.wcslim and self.xlim!=None and self.ylim!=None) :
      xcent = (self.xlim[0]+self.xlim[1])/2.
      ycent = (self.ylim[0]+self.ylim[1])/2.
      self.ra, self.de = self.wcs.pixel_to_world_values(xcent-1,ycent-1) # -1 because of an astropy python array inconsistency
      ra = np.zeros(4)
      de = np.zeros(4)
      ra[0], de[0] = self.wcs.pixel_to_world_values(self.xlim[0]-1,ycent-1)  # -1 because of an astropy python array inconsistency
      ra[1], de[1] = self.wcs.pixel_to_world_values(self.xlim[1]-1,ycent-1)  # -1 because of an astropy python array inconsistency
      ra[2], de[2] = self.wcs.pixel_to_world_values(xcent-1,self.ylim[0]-1)  # -1 because of an astropy python array inconsistency
      ra[3], de[3] = self.wcs.pixel_to_world_values(xcent-1,self.ylim[1]-1)  # -1 because of an astropy python array inconsistency
      self.ra_rad = np.abs( ra[1] - ra[0] ) / 2.
      self.de_rad = np.abs( de[3] - de[2] ) / 2.
        
  def on_ylims_change(self,event):
    self.ylim =  event.get_ylim()
    if self.verbose>=3 :
        print('ylimits     {:<10.3f} {:<10.3f}'.format(self.ylim[0],self.ylim[1]))
    if (self.wcslim and self.xlim!=None and self.ylim!=None) :
      xcent = (self.xlim[0]+self.xlim[1])/2.
      ycent = (self.ylim[0]+self.ylim[1])/2.
      self.ra, self.de = self.wcs.pixel_to_world_values(xcent-1,ycent-1) # -1 because of an astropy python array inconsistency
      ra = np.zeros(4)
      de = np.zeros(4)
      ra[0], de[0] = self.wcs.pixel_to_world_values(self.xlim[0]-1,ycent-1)  # -1 because of an astropy python array inconsistency
      ra[1], de[1] = self.wcs.pixel_to_world_values(self.xlim[1]-1,ycent-1)  # -1 because of an astropy python array inconsistency
      ra[2], de[2] = self.wcs.pixel_to_world_values(xcent-1,self.ylim[0]-1)  # -1 because of an astropy python array inconsistency
      ra[3], de[3] = self.wcs.pixel_to_world_values(xcent-1,self.ylim[1]-1)  # -1 because of an astropy python array inconsistency
      self.ra_rad = np.abs( ra[1] - ra[0] ) / 2.
      self.de_rad = np.abs( de[3] - de[2] ) / 2.
    
  def prev(self, event):
    #print('prev')
    if (self.imagenr>0) :
      self.imagenr = self.imagenr - 1
    else :
      self.imagenr = self.imagelistend - 1
    if self.fixcuts :
      self.load(self.cut0,self.cut1)
    else :
      self.load()
    
  def next(self, event):
    #print('next')
    #print('self.imagenr=',self.imagenr)
    if (self.imagenr<self.imagelistend-1) :
      self.imagenr = self.imagenr + 1
    else :
      self.imagenr = 0          
    if self.fixcuts :
      self.load(self.cut0,self.cut1)
    else :
      self.load()

    
  def onclick(self, event):

    if (self.verbose>=5) :
      print('%s click: button=%s, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))

    #if (event.key=='4' and self.key4active) : # workaround for button_press_bug
    #  return
    #else :
    #  self.key4active=False
      
    #if not event.inaxes == self.ax:
    #  return

    #if (event.button==1) :
      # print('left mouse')


    #if (event.button==1) :
    #  print('left mouse')

    if (event.button==2)  : # Previous
      if (self.imagenr>0) :
        self.imagenr = self.imagenr - 1
      #else :
      #  self.imagenr = self.imagelistend - 1 
      if self.fixcuts :
        self.load(self.cut0,self.cut1)
      else :
        self.load()
    
    elif (event.button==3)  : # next
      if (self.imagenr<self.imagelistend-1) :
        self.imagenr = self.imagenr + 1
      #else :
      #  self.imagenr = 0          
      if self.fixcuts :
        self.load(self.cut0,self.cut1)
      else :
        self.load()


    else :
      pass  
      
    self.ax.set_navigate(True)

    
  def keypress(self, event):
    #print('event=',event)
    if (self.verbose>=5) :
        print('event.key = ',event.key)

        
    if (event.key=='1')  :
        #print(event)
        pos = [event.xdata, event.ydata]
        self.xc, self.yc = self.get_cursor(pos,verbose=self.verbose)
        #if self.zoom!=0 :
        #  self.ax.xaxis.zoom(self.zoom) 
        #  self.ax.yaxis.zoom(self.zoom)
        #self.fig.canvas.draw_idle()
        #self.fig.canvas.set_message('HELLO')


    elif (event.key=='2') :
        #print(event)
        pos = [event.xdata, event.ydata]
        self.xc, self.yc = self.fit_psf(pos,use_moffat=False,verbose=self.verbose)
        #if self.zoom!=0 :
        #  self.ax.xaxis.zoom(self.zoom) 
        #  self.ax.yaxis.zoom(self.zoom)
        #self.fig.canvas.draw_idle()

        
    elif (event.key=='3') :
        pos = [event.xdata, event.ydata]
        self.xc, self.yc = self.fit_psf(pos,use_moffat=True,verbose=self.verbose)
        #if self.zoom!=0 :
        #  self.ax.xaxis.zoom(self.zoom) 
        #  self.ax.yaxis.zoom(self.zoom)
        #self.fig.canvas.draw_idle()

        
    elif (event.key=='6') : 
      pos = [event.xdata, event.ydata]
      self.xc, self.yc = self.fit_psf(pos,use_moffat=False,verbose=-1)
      if (not np.isnan(self.xc) and not np.isnan(self.yc) ) :
        xc = int(self.xc+0.5)
        yc = int(self.yc+0.5)
        rad = self.fitradius
        x0 = xc-rad
        x1 = xc+rad
        y0 = yc-rad
        y1 = yc+rad
        DATA = self.imamat[x0-1:x1,y0-1:y1]    # minus 1: fits array -> python array
        DATA_plot = np.transpose(DATA)
        Fmin = np.nanpercentile(DATA,0.1)
        Fmax = np.nanmax(DATA)
        ix,iy = np.unravel_index(np.argmax(DATA, axis=None), DATA.shape)
        print("Fmin=",Fmin)
        print("Fmax=",Fmax)
        DF = (Fmax-Fmin)/10.
        contours = np.linspace(Fmin+DF,Fmax-DF,9)
        contoursall = np.linspace(Fmin+DF/5,Fmax-DF/5,49)
        #print("contours =",contours)
        #print("contoursall =",contoursall)
        X = np.linspace(x0,x1,2*rad+1)
        Y = np.linspace(y0,y1,2*rad+1)
        XX, YY = np.meshgrid(X,Y)
        self.fig2, self.ax2 = plt.subplots(figsize=(6,6),facecolor='w')
        self.ax2.scatter(x0+ix,y0+iy, marker='x', s=10, color='red',zorder=5)
        # algorithm{'mpl2005', 'mpl2014', 'serial', 'threaded'} default: 'mpl2014')
        #myalgorithm='mpl2014'
        #myalgorithm='mpl2005'
        #mycorner_mask=False
        self.ax2.contour(XX,YY,DATA_plot,contoursall,colors='lightgrey',linewidths=0.5)
        self.ax2.contour(XX,YY,DATA_plot,contours,colors='grey',linewidths=1)
        self.ax2.contour(XX,YY,DATA_plot,[contours[4]],colors='red',linewidths=2)
        # matplotlib.pyplot.tricontour
        ### positions = np.vstack(list(zip(XX.ravel(), YY.ravel())))
        ### ri = Delaunay(positions)
        ### plt.triplot(X, Y, tri.simplices)
        ### self.ax2.tricontour(DATA,contoursall,colors='lightgrey',linewidths=0.5)
        #self.ax2.tricontour(DATA,contoursall,colors='lightgrey',linewidths=0.5)
        #self.ax2.tricontour(DATA,contours,colors='grey',linewidths=1)
        #self.ax2.tricontour(DATA,[contours[4]],colors='red',linewidths=2)
        self.ax2.set_xlabel('x')
        self.ax2.set_ylabel('y')
        self.ax2.set_aspect('equal', 'box')
        self.axtitle = self.ax2.set_title('{} {}'
                                          .format(self.imagename,'(psf contour)'))
        
        self.ax.set_navigate(False)
        self.ax2.set_navigate(True)
        if (self.autosavepng) :
          self.fig2.savefig(self.imagename+'_contour.pdf')
          self.fig2.savefig(self.imagename+'_contour.png')
        self.fig2.show()

      
    elif (event.key=='5') :
      pos = [event.xdata, event.ydata]
      self.xc, self.yc = self.fit_psf(pos,use_moffat=False,verbose=-1)
      if (not np.isnan(self.xc) and not np.isnan(self.yc) ) :
        r,p,F = self.show_psf([self.xc, self.yc],verbose=self.verbose)
        colormap = plt.cm.twilight_shifted
        colorst = [ colormap(i) for i in ((p+90.)/180.) ]
        self.fig2, self.ax2 = plt.subplots(figsize=(6,4),facecolor='w')
        # self.pointc,= self.ax2.plot(r,F, marker='o', linestyle='-', ms=3,
        #                             linewidth=1,color='b',zorder=4)
        self.pointc = self.ax2.scatter(r, F, c=colorst, marker='o', s=30, alpha=0.8,
                                       edgecolors='none')
        Fmin = np.nanpercentile(F,0.25)
        Fmax = np.nanmax(F)
        Fhalf= Fmin+(Fmax-Fmin)/2.
        rmax = np.max(r)
        self.ax2.plot([0,rmax],[Fhalf,Fhalf], marker=None, linestyle='-', ms=2,
                                     linewidth=1,color='r',zorder=5)
        self.ax2.set_xlabel('radius')
        self.ax2.set_ylabel('flux')
        self.axtitle = self.ax2.set_title('{} {}'
                                          .format(self.imagename,'(psf profile)'))
        self.ax.set_navigate(False)
        self.ax2.set_navigate(True)
        self.fig2.show()

        
    # elif (event.key=='9') :
    #   #self.key4active = True
    #   pos = [event.xdata, event.ydata]
    #   self.xc, self.yc = self.fit_psf(pos,use_moffat=False,verbose=-1)        
    #   if (np.isnan(self.xc) or np.isnan(self.yc) ) :
    #     print('BAD gaussian fit! using clicked postions')
    #   if (not np.isnan(self.xc) and not np.isnan(self.yc) ) :
    #     ixc=int(self.xc+0.5)
    #     iyc=int(self.yc+0.5)
    #     self.wcs = WCS(self.header,relax=True,fix=False)
    #     ra1, de1 = self.wcs.pixel_to_world_values(self.xc-1,self.yc-1)  # -1 because of an astropy python array inconsistency
    #     print(' >>>>> ',self.xc,self.yc,ra1, de1)
    # 
    #     for iname in self.imagelist :
    #       iheader = pyfits.getheader(iname)
    #       iwcs = WCS(iheader,relax=True,fix=False)
    #       iplot = pyfits.getdata(iname).astype(np.float32)
    #       imat = np.transpose(iplot)
    #       fx,fy = iwcs.world_to_pixel_values(ra1,de1) + np.array([1.,1.])
    #       ix,iy = int(fx+0.5), int(fy+0.5)
    #       zp50 = iheader['ZP50']
    #       flux = np.sum((imat[ix-5:ix+5,iy-5:iy+5]))
    #       sky = np.median((imat[ix+15:ix+25,iy+15:iy+25]))*100.
    #       print(iname,ix,iy,flux,sky,-2.5*np.log10(flux-sky)+zp50)

        
    elif (event.key=='4') :
      #self.key4active = True
      pos = [event.xdata, event.ydata]
      self.xc, self.yc = self.fit_psf(pos,use_moffat=False,verbose=-1)        
      if (np.isnan(self.xc) or np.isnan(self.yc) ) :
        self.xc = event.xdata
        self.yc = event.ydata
        print('BAD gaussian fit! using clicked postions')
      if (not np.isnan(self.xc) and not np.isnan(self.yc) ) :
        ixc=int(self.xc+0.5)
        iyc=int(self.yc+0.5)
        
        ####### growth curve
        prec = self.precaper
        maxaper = self.maxaper
        # self.sky0 = np.nanmean(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50])
        # self.sky0 = np.nanmedian(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50])
        #c, low, upp = sigmaclip(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50],3.,2.,
        #                        nan_policy=‘omit’)
        #cenfunc='np.nanmean',stdfunc='np.nanstd)')
        #c, low, upp = sigmaclip(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50],3.,2.)
        #self.sky0 = np.average(c)
        #self.sky0 = np.nanmean(c)
        #self.sky0 = np.nanpercentile(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50],0.45)
        self.sky0 = np.nanmedian(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50])
        #print(' self.sky0=', self.sky0)
        #print('c, low, upp  = ',np.average(c), low, upp )
        #print(' self.sky0=', self.sky0)
        self.r = np.linspace(1.,maxaper,prec)
        #print(self.r)
        aper2 = self.r**2+1.
        self.grow = np.zeros(prec)
        self.area = np.zeros(prec)        
        self.growsky = np.zeros(prec)
        for i in range(prec) :
          s=0.
          n=0
          #print('{:10} {:10.0f} {:10.0f} {}'.format('test mask',ixc,iyc,
          #                                          self.imamat[ixc-1,iyc-1]))
          for ix in range(ixc-maxaper-1,ixc+maxaper+1) :
            for iy in range(iyc-maxaper-1,iyc+maxaper+1) :
              r2 = (ix-self.xc)**2+(iy-self.yc)**2
              # if r2<=aper2[i] :              
              if r2<aper2[i] :
                value = self.imamat[ix-1,iy-1]
                #print('value=',value)
                if not np.isnan(value) :
                  #print('value=',value)                  
                  s += value - self.sky0
                  n += 1
                #else :
                #  print('nan=',value)
          self.grow[i] = s
          self.area[i] = n
          self.growsky[i] = s+n*self.sky0
        #print('self.grow=',self.grow)
        self.maxgrow = np.nanmax(self.grow)
        grow_min = np.nanmin(self.grow)
        grow_max = np.nanmax(self.grow)
        grow_del = grow_max - grow_min
        #print('grow_min=',grow_min)
        #print('grow_max=',grow_max)
        self.ylim0 = grow_min-0.1*grow_del
        self.ylim1 = grow_max+0.1*grow_del
        
        self.fig2, self.ax2 = plt.subplots(figsize=(6,4),facecolor='w')
        self.lineg, = self.ax2.plot(self.r,self.grow,color='lightblue',zorder=1)
        self.linef, = self.ax2.plot([None,None],[None,None],color='mistyrose',zorder=0)
        self.linec, = self.ax2.plot(self.r,self.grow,color='b',zorder=3)
        self.pointc,= self.ax2.plot(self.r,self.grow, marker='o', linestyle=None, ms=3,
                                    linewidth=0,color='b',zorder=4)
        self.line1, = self.ax2.plot([None,None],[None,None],color='grey',zorder=2)
        self.line2, = self.ax2.plot([None,None],[None,None],color='grey',zorder=2)
        self.lineh, = self.ax2.plot([None,None],[None,None],color='r',   zorder=2)
        #linef, = plt.plot([p1.x_fwhm,p1.x_fwhm],[0.,1.2*p1.totflux],color='g',zorder=0)
        self.ax2.set_xlabel('radius')
        self.ax2.set_ylabel('flux')
        self.axtitle = self.ax2.set_title('{} {:10.1f} {}'
                                          .format(self.imagename,self.totflux,'(growth curve)'))
        self.ax2.set_ylim([self.ylim0,self.ylim1])
        
        self.click2 = False
        #self.id2_cursor = Cursor( self.ax2, useblit=False, color='g', linewidth=1 )
        self.id2_button_press = self.fig2.canvas.mpl_connect('button_press_event',self.growth_curve)
        self.id2_close        = self.fig2.canvas.mpl_connect('close_event',self.growth_curve_close)
        #cid4 = self.fig2.canvas.mpl_connect('key_press_event',self.keypress)
        if self.verbose>=2 :
          print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'
                .format('#','xc', 'yc', 'totflux', 'sky', 'r0', 'r1', 'errtotflux' ))
        if self.fulloutput :
          self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'
                             .format('#','method','xc', 'yc', 'totflux',
                                     'sky', 'r0', 'r1', 'errtotflux' ))
    
        #self.fig.canvas.mpl_disconnect(self.id1)  
        self.ax.set_navigate(False)
        self.ax2.set_navigate(True)
        self.fig2.show()
        ##self.fig.canvas.draw_idle()
    
      else :
        pass  
        
    elif (event.key=='0')  : # mask circle
      self.imamat0 = np.array(self.imamat)
      self.xc = int(event.xdata+0.5)
      self.yc = int(event.ydata+0.5)
      if self.verbose>=2 :
        print('{:<10} {:>10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'oldvalue', 'value'))
      if self.verbose>=1 :         
        print('{:10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}'
              .format('mask',self.xc,self.yc,self.imamat[self.xc-1,self.yc-1],np.nan))
      if self.fulloutput :
        self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10}\n'
                           .format('#','method','xc', 'yc', 'oldvalue', 'value'))
        self.outfile.write('{:<30} {:<10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}\n'
                           .format(self.imagename,'mask',self.xc,self.yc,
                                   self.imamat[self.xc-1,self.yc-1],np.nan))
      else :
        self.outfile.write('{:<30} {:10.0f} {:10.0f} {:10.2f} {:10.2f}\n'
                           .format(self.imagename,self.xc,self.yc,
                                   self.imamat[self.xc-1,self.yc-1],np.nan))
      mx0 = self.xc-self.mask_radius-1
      mx1 = self.xc+self.mask_radius+1
      my0 = self.yc-self.mask_radius-1
      my1 = self.yc+self.mask_radius+1
      if mx0<0       : mx0=0
      if mx1>self.nx : mx1=self.nx
      if my0<0       : my0=0
      if my1>self.ny : my1=self.ny
      r2max = (self.mask_radius+0.4)**2 # 0.4 for a convenient cross for small rad
      for ix in range(mx0,mx1+1) :
        for iy in range(my0,my1+1) :
            r2 = (ix-self.xc)**2+(iy-self.yc)**2
            if r2<=r2max :
              self.imamat[ix-1,iy-1] = np.nan 
      #print('imamat',self.imamat[mx0-1:mx1,my0-1:my1])
      #if self.verbose>=1 :         
      #  print('{:10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}'.format('mask',self.xc,self.yc,
      #        self.imamat[self.xc-1,self.yc-1],np.nan))
      #self.imshow = self.ax.imshow( self.imamat_plot, origin='lower', cmap=self.cmap,
      #                              vmin=self.cut0, vmax=self.cut1, extent=self.extent)
      self.imshow.set_data(self.imamat_plot)

      self.fig.canvas.draw_idle()
      #self.load()


    elif (event.key=='b')  : # mask circle
      self.imamat0 = np.array(self.imamat)
      self.xc = int(event.xdata+0.5)
      self.yc = int(event.ydata+0.5)
      if self.verbose>=2 :
        print('{:<10} {:>10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'oldvalue', 'value'))
      if self.verbose>=1 :         
        print('{:10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}'
              .format('mask',self.xc,self.yc,self.imamat[self.xc-1,self.yc-1],np.nan))
      if self.fulloutput :
        self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10}\n'
                           .format('#','method','xc', 'yc', 'oldvalue', 'value'))
        self.outfile.write('{:<30} {:<10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}\n'
                           .format(self.imagename,'mask',self.xc,self.yc,
                                   self.imamat[self.xc-1,self.yc-1],np.nan))
      else :
        self.outfile.write('{:<30} {:10.0f} {:10.0f} {:10.2f} {:10.2f}\n'
                           .format(self.imagename,self.xc,self.yc,
                                   self.imamat[self.xc-1,self.yc-1],np.nan))
      mx0 = self.xc-50-self.mask_radius-1
      mx1 = self.xc+50+self.mask_radius+1
      my0 = self.yc-50-self.mask_radius-1
      my1 = self.yc+50+self.mask_radius+1
      if mx0<0       : mx0=0
      if mx1>self.nx : mx1=self.nx
      if my0<0       : my0=0
      if my1>self.ny : my1=self.ny
      r2max = (50+self.mask_radius+0.4)**2 # 0.4 for a convenient cross for small rad
      for ix in range(mx0,mx1+1) :
        for iy in range(my0,my1+1) :
            r2 = (ix-self.xc)**2+(iy-self.yc)**2
            if r2<=r2max :
              self.imamat[ix-1,iy-1] = np.nan 
      #print('imamat',self.imamat[mx0-1:mx1,my0-1:my1])
      #if self.verbose>=1 :         
      #  print('{:10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}'.format('mask',self.xc,self.yc,
      #        self.imamat[self.xc-1,self.yc-1],np.nan))
      #self.imshow = self.ax.imshow( self.imamat_plot, origin='lower', cmap=self.cmap,
      #                              vmin=self.cut0, vmax=self.cut1, extent=self.extent)
      self.imshow.set_data(self.imamat_plot)

      self.fig.canvas.draw_idle()
      #self.load()


    elif (event.key=='8') or (event.key=='7')  : # satellite
      if ( self.satx0 == None and self.saty0 == None ): 
        self.imamat0 = np.array(self.imamat)
        if (event.key=='7') :
          self.sat0 = True
        else :
          self.sat0 = False
        self.satx0 = int(event.xdata+0.5)
        self.saty0 = int(event.ydata+0.5)
        if self.verbose>=1 :         
          print('{:10} {:10.0f} {:10.0f}'
                .format('sat0',self.satx0,self.saty0))
        ###### mark start point
        mx0 = self.satx0-self.mask_radius-1
        mx1 = self.satx0+self.mask_radius+1
        my0 = self.saty0-self.mask_radius-1
        my1 = self.saty0+self.mask_radius+1
        if mx0<0       : mx0=0
        if mx1>self.nx : mx0=self.nx
        if my0<0       : mx0=0
        if my1>self.ny : mx0=self.ny
        r2max = (self.mask_radius+0.4)**2 # 0.4 for a convenient cross for small rad
        for ix in range(mx0,mx1+1) :
          for iy in range(my0,my1+1) :
              r2 = (ix-self.satx0)**2+(iy-self.saty0)**2
              if r2<=r2max :
                self.imamat[ix-1,iy-1] = np.nan
        self.imshow.set_data(self.imamat_plot)
        self.fig.canvas.draw_idle()
        #################
      else :
        if (event.key=='7') :
          self.sat1 = True
        else :
          self.sat1 = False
        self.satx1 = int(event.xdata+0.5)
        self.saty1 = int(event.ydata+0.5)
        if (self.satx1-self.satx0) != 0 :
          ax = (self.saty1-self.saty0) / float(self.satx1-self.satx0)
          bx = self.saty1 - ax*self.satx1
        else :
          ax = np.nan
          bx = np.nan
        if (self.saty1-self.saty0) != 0 :
          ay = (self.satx1-self.satx0) / float(self.saty1-self.saty0)
          by = self.satx1 - ay*self.saty1
        else :
          ay = np.nan
          by = np.nan
        if self.verbose>=1 :         
          print('{:10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}'
                .format('sat1',self.satx1,self.saty1,ax,bx))
        ###### mark end point
        mx0 = self.satx1-self.mask_radius-1
        mx1 = self.satx1+self.mask_radius+1
        my0 = self.saty1-self.mask_radius-1
        my1 = self.saty1+self.mask_radius+1
        if mx0<0       : mx0=0
        if mx1>self.nx : mx0=self.nx
        if my0<0       : mx0=0
        if my1>self.ny : mx0=self.ny
        r2max = (self.mask_radius+0.4)**2 # 0.4 for a convenient cross for small rad
        for ix in range(mx0,mx1+1) :
          for iy in range(my0,my1+1) :
              r2 = (ix-self.satx1)**2+(iy-self.saty1)**2
              if r2<=r2max :
                self.imamat[ix-1,iy-1] = np.nan
        self.imshow.set_data(self.imamat_plot)
        self.fig.canvas.draw_idle()
        #################
        startx = 0
        starty = 0
        endx = self.nx
        endy = self.ny
        if self.sat0 : 
          if self.satx0<=self.satx1 :
            startx = self.satx0
          else :
            endx   = self.satx0
          if self.saty0<=self.saty1 :
            starty = self.saty0
          else :
            endy   = self.saty0
        if self.sat1 : 
          if self.satx0<=self.satx1 :
            endx   = self.satx1
          else :
            startx = self.satx1
          if self.saty0<=self.saty1 :
            endy   = self.saty1
          else :
            starty = self.saty1
        if ax<1 :
          for ix in range(startx,endx+1) :
            my0 = int(ax*ix + bx - self.mask_radius)
            my1 = int(ax*ix + bx + self.mask_radius)
            if my0<1       : my0=1
            if my1<1       : my1=1
            if my0>self.ny : my0=self.ny
            if my1>self.ny : my1=self.ny
            self.imamat[ix-1,my0-1:my1] = np.nan 
        if ay<1 :
          for iy in range(starty,endy+1) :
            mx0 = int(ay*iy + by - self.mask_radius)
            mx1 = int(ay*iy + by + self.mask_radius)
            if mx0<1       : mx0=1
            if mx1<1       : mx1=1
            if mx0>self.nx : mx0=self.nx
            if mx1>self.nx : mx1=self.nx
            self.imamat[mx0-1:mx1,iy-1] = np.nan 
        self.satx0 = None
        self.saty0 = None
        self.satx1 = None
        self.saty1 = None
        self.imshow.set_data(self.imamat_plot)
        self.fig.canvas.draw_idle()

    if (event.key=='@') :
        print('UNDO')
        self.imamat = np.array(self.imamat0)
        self.imamat_plot = np.transpose(self.imamat)
        self.imshow.set_data(self.imamat_plot)
        self.fig.canvas.draw_idle()

        
    if (event.key=='z' or event.key=='Z' or event.key=='x' or event.key=='X' or
        event.key=='m' or event.key=='n' or event.key=='j' ) :              
      self.dcut = self.cut1 - self.cut0
      if (event.key=='z') :
        self.cut0n(self)
        #self.cut0 = self.cut1 - self.dcut*self.cut_zoom        
      if (event.key=='x') :
        self.cut1n(self)
        #self.cut0 = self.cut1 - self.dcut/self.cut_zoom        
      if (event.key=='Z') :
        self.cut0p(self)
        #self.cut1 = self.cut0 + self.dcut/self.cut_zoom
      if (event.key=='X') :
        self.cut1p(self)
        #self.cut1 = self.cut0 + self.dcut*self.cut_zoom        
      if (event.key=='m') :
        self.save_mask(self)
      if (event.key=='n') :
        self.save_png(self)
      if (event.key=='j') :
        self.save_jpg(self)
      #if (event.key=='0') :
      #  self.cut0 = 0.
      #AR print('cuts        {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
      #self.imshow = self.ax.imshow( self.imamat_plot, origin='lower',
      #                              cmap=self.cmap, vmin=self.cut0, vmax=self.cut1)
      #AR self.imshow = self.ax.imshow( self.imamat_plot, origin='lower', cmap=self.cmap,
      #AR                              vmin=self.cut0, vmax=self.cut1, extent=self.extent)
      #self.fig.canvas.draw_idle()
    elif event.key==' ' :
      if verbose>=1 :
        print('#########################################################################################################')
        #if self.fulloutput :
        self.outfile.write('#########################################################################################################\n')
    # elif event.key=='o' :
    #   xn, yn = np.linspace(0.5,self.nx+0.5, self.nx, dtype=float), np.linspace(0.5,self.ny+0.5, self.ny, dtype=float)
    #   Xn, Yn = np.meshgrid(xn,yn)
    #   Z = self.imamat_plot
    #   ima_min = np.nanmin(self.imamat)
    #   ima_max = np.nanmax(self.imamat)
    #   ima_del = ima_max - ima_min
    #   #levels = np.linspace(ima_min+0.05*ima_del, ima_max-0.05*ima_del, 20)
    #   levels = np.linspace(ima_min, ima_max, 20)
    #   idc = self.ax.contour(Xn, Yn, Z, levels, colors='g')
    #   self.fig.canvas.draw_idle()
    else :
      pass  

    self.ax.set_navigate(True)
      
        
########################### MAIN ###########################

parser = argparse.ArgumentParser(description=VERSION+'\n')
parser.add_argument(              'imagelist',          nargs='*',                     help='image list' )
parser.add_argument('-figsize',   dest='figsize',       type=str,  default='10.0,8.5', help='[%(default)s] figsize' )
parser.add_argument('-v',         dest='verbose',       type=int,  default=2,          help='[%(default)s] verbose' )
parser.add_argument('-precaper',  dest='precaper',      type=int,  default=80,         help='[%(default)s] precaper' )
parser.add_argument('-maxaper',   dest='maxaper',       type=int,  default=80,         help='[%(default)s] maxaper' )
parser.add_argument('-fitradius', dest='fitradius',     type=int,  default='15',       help='[%(default)s] fitradius' )
parser.add_argument('-zoom',      dest='zoom',          type=int,  default='0',        help='[%(default)s] zoom' )
parser.add_argument('-interpol',  dest='interpol',      type=int,  default='1',        help='[%(default)s] interpolation' )
parser.add_argument('-cmap',      dest='cmap',          type=int,  default='102',      help='[%(default)s] cmap' )
parser.add_argument('-full',      dest='fulloutput',    action='store_true',           help='[%(default)s] full output' )
parser.add_argument('-automask',  dest='autosavemasks', action='store_true',           help='[%(default)s] autosave all masks' )
parser.add_argument('-autopng',   dest='autosavepng',   action='store_true',           help='[%(default)s] autosave all images as png' )
parser.add_argument('-autojpg',   dest='autosavejpg',   action='store_true',           help='[%(default)s] autosave all images as jpg' )
parser.add_argument('-wcslim',    dest='wcslim',        action='store_true',           help='[%(default)s] keep lim according WCS' )
parser.add_argument('-fixcuts',   dest='fixcuts',       action='store_true',           help='[%(default)s] keep cuts fixed from previous image' )
parser.add_argument('-cuts',      dest='cuts',          type=str,  default='0,0',      help='[%(default)s] cuts low and high' )

args = parser.parse_args()


figsize   = np.array(args.figsize.split(','),dtype=float)
cuts      = np.array(args.cuts.split(','),dtype=float)
verbose   = args.verbose
fitradius = args.fitradius
precaper  = args.precaper
maxaper   = args.maxaper
interpol  = None
methods = ['none', 'nearest', 'bilinear', 'bicubic', 'spline16',
           'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
           'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
if args.interpol != None :
  interpol = methods[args.interpol]



if args.imagelist==[] :
    parser.print_help()
    print('\nusage in plot window:')
    print('          1 - get cursor    at mouse position')
    print('          2 - center gauss  at mouse position')
    print('          3 - center moffat at mouse position')
    print('          4 - growth curve  at mouse position')
    print('          5 - psf profile   at mouse position')
    print('          6 - star contour  at mouse position')
    print('          7 - mask     satellite with two mouse positions')
    print('          8 - mask inf satellite with two mouse positions')

    print('          0 - mask circle   at mouse position')
    print('          b - big mask circle   at mouse position')
    print('          m - save mask')
    print('          n - save mask')
    print('          j - save mask')
    print('          @ - UNDO last mask command')    
    print('          z - decrease lower cut')
    print('          Z - increase lower cut')
    print('          x - decrease higher cut')
    print('          X - increase higher cut')
    print('        h r - reset (zoom home)')
    print('        < c - zoom back')
    print('        > v - zoom forward')
    print('          p - pan (move center)')
    print('          o - Zoom-to-rect')
    print('          s - save image')
    print('          f - full frame')
    print('          g - Toggle major grids')
    print('          G - Toggle minor grids')
    print('        k L - x axis log/linear')
    print('          l - y axis log/linear')
    print('          q - quit')
    exit(0)

# see    https://matplotlib.org/stable/users/explain/interactive.html#key-event-handling
# h r  < c  > v  p  o  s  f  g  G  k L  l  q



imagelist = args.imagelist
# outfile  = open('/project/ls-gruen/ligodesi/work/showfits.tab', 'w')
outfile  = open( 'showfits.tab' , 'w' )
fig = plt.figure(figsize=(figsize[0],figsize[1]),facecolor='w',dpi=100)

ia = Image_Analyzer(fitradius=fitradius,
                    precaper=precaper,
                    maxaper=maxaper,
                    zoom=args.zoom,
                    interpolation=args.interpol,                    
                    cmap=args.cmap,
                    fixcuts=args.fixcuts,
                    cuts=cuts,
                    wcslim=args.wcslim,
                    verbose=verbose,
                    fulloutput=args.fulloutput,
                    outfile=outfile,
                    autosavemasks=args.autosavemasks,
                    autosavepng=args.autosavepng,
                    autosavejpg=args.autosavejpg)


ia.start(fig)
ia.storelist(imagelist)
if cuts[0]!=0 or cuts[1]!=0 :
  ia.load(cuts[0],cuts[1])
else :
  ia.load()

#ia.load(-400,400)
ia.connect()

plt.show()

ia.disconnect()

outfile.close()
  

