#!/Users/arri/miniconda3/bin/python
## !/Users/arri/miniconda3/bin/python
## !/opt/miniconda3/bin/python
## !/afs/ipp/.cs/anaconda/amd64_generic/3/2019.10/bin/python
## !/afs/ipp-garching.mpg.de/home/a/arri/miniconda3/bin/python
## !/opt/miniconda3/bin/python
## !/opt/anaconda3.7/bin/python


# (2021-03-13) copyright by Arno Riffeser (arri@usm.lmu.de)


import matplotlib as mpl
#mpl.use('TkAgg')   # ATTENTION: TkAgg on MacOSX: button_press_event not working correctly!!
mpl.use('MacOSX')

from math import *
import numpy as np
import argparse
import copy
from   scipy.stats     import sigmaclip
from   astropy.io      import fits as pyfits
from   scipy.optimize  import leastsq, curve_fit,least_squares
import matplotlib.pyplot  as plt
import matplotlib.mlab    as mlab
from   matplotlib.widgets import Cursor, Slider, Button, RadioButtons


axcolor = 'lightblue'


class Image_Analyzer:

  def __init__(self, fitradius=7, precaper=80, maxaper=80, zoom=0, verbose=2, divvy=False, outfile=None):
    self.verbose = verbose
    #self.ax = axlist[0]
    #self.axprev = axlist[1]
    #self.axnext = axlist[2]
    self.fig = None
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
    self.ylim0 = 0.
    self.ylim1 = 0.
    self.imagelist = []
    self.imagelistend = 0
    self.imagenr = 0
    self.cmap = None
    self.cmap_nr = 0
    self.cmap_all = ["gist_heat","gist_heat_r","gist_gray","gist_gray_r","gist_ncar","gist_ncar_r","gist_rainbow","gist_rainbow_r","gist_stern","gist_stern_r","gist_yarg","gist_yarg_r","gist_earth","gist_earth_r"]
    self.cmap_listend = len(self.cmap_all)
    self.divvy = divvy
    self.zoom = zoom
    self.xlim = None
    self.ylim = None
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
    self.idcmap = None
    self.ax = None
    self.ax2 = None
    self.axprev = None
    self.axnext = None
    self.axcut0n = None
    self.axcut0p = None
    self.axcut1n = None
    self.axcut1p = None
    self.axcmap  = None
    
      
  def start(self, fig) :
    self.ax      = fig.add_axes([0.07, 0.07, 0.91, 0.91])
    self.axprev  = fig.add_axes([0.01, 0.01, 0.10, 0.03])
    self.axnext  = fig.add_axes([0.12, 0.01, 0.10, 0.03])
    self.axcut0n = fig.add_axes([0.23, 0.01, 0.06, 0.03])
    self.axcut0p = fig.add_axes([0.30, 0.01, 0.06, 0.03])
    self.axcut1n = fig.add_axes([0.37, 0.01, 0.06, 0.03])
    self.axcut1p = fig.add_axes([0.44, 0.01, 0.06, 0.03])
    self.axcmap  = fig.add_axes([0.51, 0.01, 0.06, 0.03])
  
  def connect(self):
    #print(self.imshow)
    self.id1 = self.imshow.figure.canvas.mpl_connect('button_press_event',self.onclick)
    # self.id1 = self.ax.figure.canvas.mpl_connect('button_press_event',self.onclick)
    # self.id1 = self.ax.on_clicked(self.onclick)
    self.id2 = self.imshow.figure.canvas.mpl_connect('key_press_event',self.keypress)
    # self.id3 = self.imshow.figure.canvas.mpl_connect('close_event',self.close)
    # cursor = Cursor(self.ax, useblit=False, color='g', linewidth=1 )
    self.idprev = Button(self.axprev, 'previous')
    self.idprev.on_clicked(self.prev)
    self.idnext = Button(self.axnext, 'next')
    self.idnext.on_clicked(self.next)
    self.idcut0n = Button(self.axcut0n, 'cut0 -')
    self.idcut0n.on_clicked(self.cut0n)
    self.idcut0p = Button(self.axcut0p, 'cut0 +')
    self.idcut0p.on_clicked(self.cut0p)
    self.idcut1n = Button(self.axcut1n, 'cut1 -')
    self.idcut1n.on_clicked(self.cut1n)
    self.idcut1p = Button(self.axcut1p, 'cut1 +')
    self.idcut1p.on_clicked(self.cut1p)
    self.idcmap = Button(self.axcmap, 'cmap')
    self.idcmap.on_clicked(self.change_cmap)


  def disconnect(self):
    self.imshow.figure.canvas.mpl_disconnect(self.id1)
    #self.imshow.figure.canvas.mpl_disconnect(self.id2)
    #self.imshow.figure.canvas.mpl_disconnect(self.id3)
    #self.ax.callbacks.disconnect(self.idxlim)
    #self.ax.callbacks.disconnect(self.idylim)

      
  def close(self, event) :
    input("really close? <enter>")

      
  def cut0n(self, event):
    self.dcut = self.cut1 - self.cut0
    self.cut0 = self.cut1 - self.dcut*self.cut_zoom     
    print('cuts           {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower',
                                  cmap=self.cmap, vmin=self.cut0, vmax=self.cut1, extent=self.extent)
    self.imshow.figure.canvas.draw_idle()

  def cut0p(self, event):
    self.dcut = self.cut1 - self.cut0
    self.cut0 = self.cut1 - self.dcut/self.cut_zoom 
    print('cuts           {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower',
                                  cmap=self.cmap, vmin=self.cut0, vmax=self.cut1, extent=self.extent)
    self.imshow.figure.canvas.draw_idle()

  def cut1n(self, event):
    self.dcut = self.cut1 - self.cut0
    self.cut1 = self.cut0 + self.dcut/self.cut_zoom   
    print('cuts           {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower',
                                  cmap=self.cmap, vmin=self.cut0, vmax=self.cut1, extent=self.extent)
    self.imshow.figure.canvas.draw_idle()

  def cut1p(self, event):
    self.dcut = self.cut1 - self.cut0
    self.cut1 = self.cut0 + self.dcut*self.cut_zoom 
    print('cuts           {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower',
                                  cmap=self.cmap, vmin=self.cut0, vmax=self.cut1, extent=self.extent)
    self.imshow.figure.canvas.draw_idle()

  def change_cmap(self, event):
    if (self.cmap_nr<self.cmap_listend-1) :
      self.cmap_nr = self.cmap_nr + 1
    else :
      self.cmap_nr = 0
    print('cmap           {:<10} {:<}'.format(self.cmap_nr,self.cmap_all[self.cmap_nr]))
    self.cmap = copy.copy(mpl.cm.get_cmap(self.cmap_all[self.cmap_nr]))
    self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower',
                                  cmap=self.cmap, vmin=self.cut0, vmax=self.cut1, extent=self.extent)
    self.imshow.figure.canvas.draw_idle()    
      
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
        print('{:10} {:10.0f} {:10.0f} {:10.2f}'.format('pixel',xc,yc,self.imamat[xc-1,yc-1]))
    if self.divvy :
      self.outfile.write('{:<30} {:10.0f} {:10.0f} {:10.2f}\n'.format(self.imagename,xc,yc,self.imamat[xc-1,yc-1]))
    else :
      self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10}\n'.format('#','method','xc', 'yc', 'value'))
      self.outfile.write('{:<30} {:<10} {:10.0f} {:10.0f} {:10.2f}\n'.format(self.imagename,'pixel',xc,yc,self.imamat[xc-1,yc-1]))
    self.ax.plot(xc,yc,marker="2",c='k',ms=10)
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
    # print("sky=",sky)
    # DATA = self.imamat[x0-1:x1,y0-1:y1] - sky    # minus 1 in den indizes wegen image array -> python array
    DATA = self.imamat[x0-1:x1,y0-1:y1]    # minus 1 in den indizes wegen image array -> python array
    bg = np.nanmedian(DATA)
    max = np.nanmax(DATA)
    amp = max-bg
    cen = DATA[rad,rad]
    # x, y = np.linspace(x0, x1, boxsize, dtype=float), np.linspace(y0, y1,boxsize, dtype=float)
    x = np.linspace(x0, x1, boxsize, dtype=float)
    y = np.linspace(y0, y1, boxsize, dtype=float)
    # print("x,y=",x,y)
    XV, YV = np.meshgrid(y,x)
    # for ix in range(0,boxsize) :
    #     for iy in range(0,boxsize) :
    #         if DATA[ix,iy]==MD :
    #             xc = x[ix]
    #             yc = y[iy]

    X_in = XV.ravel()
    Y_in = YV.ravel()
    Z_in = DATA.ravel()
    bad = np.where( np.isnan(Z_in) )
    X = np.delete(X_in,bad)
    Y = np.delete(Y_in,bad)
    Z = np.delete(Z_in,bad)
    XY = np.vstack((X,Y))
    
    if not use_moffat : # gauss
      #astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp )
      astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp , bg )
      res = least_squares(self.residuals_gauss, astart, method='lm', args=(XY,Z), verbose=0)
      #afit = np.zeros(7)
      #afit[0:6] = res.x
      #afit[6] = 0.
      afit = res.x
      # double ert = -( a(3) * tx * tx + a(4) * ty * ty );  # usmlib.h
      # fwhm0  = 2. * sqrt( M_LN2 / a(3) );                 # starphot.cpp
      # refa(0)  =  refa(2)*M_PI/sqrt(refa(3)*refa(4));     # starphot.cpp
      fwhmx = abs(afit[2]) * np.sqrt(8.*log(2.))
      fwhmy = abs(afit[3]) * np.sqrt(8.*log(2.))
      fwhm = np.sqrt(fwhmx*fwhmy)
      a3 = 1./(2.*afit[2]**2)
      a4 = 1./(2.*afit[3]**2)
      totflux = afit[5] * np.pi / np.sqrt( a3 * a4);
      if verbose>=5 :
        xn, yn = np.linspace(-30., 30., 61, dtype=float), np.linspace(-30., 30., 61, dtype=float)
        Xn, Yn = np.meshgrid(yn,xn)
        yyn = Xn.ravel()
        xxn = Yn.ravel()
        xx0, yy0, sigx, sigy, phi, A, C = afit
        print("numerical totflux = ",np.sum(self.gauss(xxn, yyn, 0., 0., sigx, sigy, phi, A, 0.)))        
      if verbose>=2 :
        print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'
              .format('#','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
      if verbose>=1 :
        print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}'
              .format('gauss',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
      if self.divvy :
        self.outfile.write('{:<30} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'
                           .format(self.imagename,afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),
                                   abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))            
      else :
        if verbose>=0 :
          self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'
                             .format('#','method','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
          self.outfile.write('{:<30} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'
                             .format(self.imagename,'gauss',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),
                                     abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
      self.ax.plot(afit[0],afit[1],marker="+",c='g',ms=10)
      
    else :
      #astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp )
      astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp , bg )
      res = least_squares(self.residuals_moffat, astart, method='lm', args=(XY,Z), verbose=0) 
      #afit = np.zeros(7)
      #afit[0:6] = res.x
      #afit[6] = 0.
      afit = res.x
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
        xn, yn = np.linspace(-30., 30., 61, dtype=float), np.linspace(-30., 30., 61, dtype=float)
        Xn, Yn = np.meshgrid(yn,xn)
        yyn = Xn.ravel()
        xxn = Yn.ravel()
        xx0, yy0, sigx, sigy, phi, A, C = afit
        print("numerical totflux = ",np.sum(self.moffat(xxn, yyn, 0., 0., sigx, sigy, phi, A, 0.)))        
      if verbose>=2 :
        print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
      if verbose>=1 :
        print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}'.format('moffat',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
      if self.divvy :
          self.outfile.write('{:<30} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'.format(self.imagename,afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
      else :
        if verbose>=0 :
          self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'.format('#','method','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
          self.outfile.write('{:<30} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'.format(self.imagename,'moffat',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
      self.ax.plot(afit[0],afit[1],marker="x",c='b',ms=10)
    return afit[0],afit[1]


  def growing_curve(self, event) :
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
      # fitlist =  opt.curve_fit(sky, rfit,growfit, x0, sigma)
      # A = np.vstack([rfit**2*np.pi,np.ones(len(rfit))]).T
      A = np.vstack([areafit,np.ones(len(rfit))]).T
      lstsqfit = np.linalg.lstsq(A, growfit,rcond=None)
      #print("lstsqfit=",lstsqfit)
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
        print("rad0 =",self.rad0)
        print("rad1 =",self.rad1)
        print("totflux        = ",self.totflux)
        print("new_growfit[0] = ",new_growfit[0])
        print("totflux+sky*rad0**2*np.pi = ",self.totflux+self.sky*self.rad0**2*np.pi)
        #print("totflux+sky*areafit[0]    = ",new_growfit[0]+self.sky*areafit[0])
        #print("totflux+sky*pi*r0^2       = ",new_growfit[0]+self.sky*areafit[0])
        print("growskyfit[0]             = ",growskyfit[0])
        #print("totflux+sky*rfit[0]**2*pi = ",self.totflux+self.sky*rfit[0]**2*np.pi)
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
      ax2title = '{}  {:.1f}  {}'.format(self.imagename,self.totflux,'(growing curve)')
      self.ax2.set_title(ax2title)
      self.lineh.figure.canvas.draw_idle()
      #self.click2 = False
      if (event.dblclick) :
        print('---------------------------------------------------------------------------------------')
      if self.verbose>=1 :         
        print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f} {:10.1f}'
              .format('grow.curve',self.xc,self.yc,self.totflux,self.sky,self.rad0,self.rad1,errtotall))
      if (event.dblclick) :
        print('---------------------------------------------------------------------------------------')
        if self.divvy :
          self.outfile.write('{:<30} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f} {:10.1f}\n'
                             .format(self.imagename,self.xc,self.yc,self.totflux,self.sky,self.rad0,self.rad1,errtotall))
        else :
          self.outfile.write('{:<30} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f} {:10.1f}\n'
                             .format(self.imagename,'grow.curve',self.xc,self.yc,self.totflux,self.sky,self.rad0,self.rad1,errtotall))
    #def growing_curve_close(self, event) :
    #print('---------------------------------------------------------------------------------------')
    #if self.rad0 < self.rad1 :
    #  if self.verbose>=2 :
    #    print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'
    #          .format('#','xc', 'yc', 'totflux', 'sky', 'r0', 'r1' ))
    #  if self.verbose>=1 :         
    #    print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f}'
    #          .format('grow.curve',self.xc,self.yc,self.totflux,self.sky,self.rad0,self.rad1))
    #  self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'
    #                     .format('#','method','xc', 'yc', 'totflux', 'sky', 'r0', 'r1' ))
    #  self.outfile.write('{:<30} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f}\n'
    #                      .format(self.imagename,'grow.curve',self.xc,self.yc,self.totflux,self.sky,self.rad0,self.rad1))


  def storelist(self,imagelist) :
    self.imagelist = imagelist
    self.imagelistend = len(self.imagelist)
    print('nr_of_images  ',self.imagelistend)
    for i,s in enumerate(imagelist) :
      print('{:<3} {}'.format(i+1,s))



  def on_xlims_change(self,event):
    self.xlim = event.get_xlim()
    if self.verbose>=3 :
        print('xlimits = {:.3f} {:.3f}'.format(self.xlim[0],self.xlim[1]))

  def on_ylims_change(self,event):
    self.ylim =  event.get_ylim()
    if self.verbose>=3 :
        print('ylimits = {:.3f} {:.3f}'.format(self.ylim[0],self.ylim[1]))
   
          
  def load(self,cut0=None,cut1=None) :
    if self.verbose>=2 :
        print('=============================================================================================================')
    self.imagename = self.imagelist[self.imagenr]
    self.title = '{}     {} of {}'.format(self.imagename,self.imagenr+1,self.imagelistend )
    self.imamat = np.transpose(pyfits.getdata(self.imagename).astype(np.float32))
    ima_med = np.nanmedian(self.imamat)
    ima_std = np.nanstd(self.imamat)
    ima_min = np.nanmin(self.imamat)
    ima_max = np.nanmax(self.imamat)
    ima_del = ima_max - ima_min
    if (cut0==None) :
      #cut0 = ima_min - 0.1 * ima_del
      cut0 = ima_med - 3. * ima_std
    if (cut1==None) :
      #cut1 = ima_max - 0.1 * ima_del
      cut1 = ima_med + 5. * ima_std
    if verbose>=3 :
      print('np.shape(self.imamat)=',np.shape(self.imamat))
    (self.nx,self.ny) = np.shape(self.imamat)
    self.cut0 = cut0
    self.cut1 = cut1
    self.dcut = cut1 - cut0
    if verbose>=1 :
      print('nr            ', self.imagenr+1,' of ',self.imagelistend)
      print('imagename     ', self.imagename)
    if self.verbose>=2 :
      print('imagesize      {:<10.0f} {:<10.0f}'.format(self.nx,self.ny))
      print('cuts           {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
    if verbose>=1 :
      print("#")
    self.extent = [0.5,self.nx+0.5,0.5,self.ny+0.5]
    # print('self.ax=',self.ax)
    self.ax.clear()
    self.cmap = copy.copy(mpl.cm.get_cmap(self.cmap_all[self.cmap_nr]))
    self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower', cmap=self.cmap, vmin=self.cut0,
                                  vmax=self.cut1, extent=self.extent)
    self.imshow.figure.canvas.set_window_title(self.title)
    current_cmap = self.imshow.get_cmap()
    current_cmap.set_bad(color='blue')
    if (self.xlim==None and self.ylim==None) :
      if self.zoom!=0 :
        self.ax.xaxis.zoom(self.zoom) 
        self.ax.yaxis.zoom(self.zoom)
    else :
      self.ax.set_xlim(self.xlim)
      self.ax.set_ylim(self.ylim)

    self.idxlim = self.ax.callbacks.connect('xlim_changed', self.on_xlims_change)
    self.idylim = self.ax.callbacks.connect('ylim_changed', self.on_ylims_change)

    # print("lim=",self.ax.get_xlim(),self.ax.get_xlim())

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
    self.imshow.figure.canvas.draw_idle()

  def prev(self, event):
    #print('prev')
    if (self.imagenr>0) :
      self.imagenr = self.imagenr - 1
    else :
      self.imagenr = self.imagelistend - 1 
    self.load()
        
  def next(self, event):
    #print('next')
    #print('self.imagenr=',self.imagenr)
    if (self.imagenr<self.imagelistend-1) :
      self.imagenr = self.imagenr + 1
    else :
      self.imagenr = 0          
    self.load()

    
  def onclick(self, event):
    # print(event.button,event.key)
    if (self.verbose>2) :
      print('%s click: button=%d, key=%s, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.key,event.x, event.y, event.xdata, event.ydata))

    #if (event.button==1)  : 
    #  print('left mouse')
      
    if (event.button==2)  : # Previous
      # print('middle mouse')
      if (self.imagenr>0) :
        self.imagenr = self.imagenr - 1
      else :
        self.imagenr = self.imagelistend - 1 
      self.load()

    if (event.button==3)  : # next
      # print('right mouse')
      # print('self.imagenr=',self.imagenr)
      if (self.imagenr<self.imagelistend-1) :
        self.imagenr = self.imagenr + 1
      else :
        self.imagenr = 0          
      self.load()

    if (event.button==1 and event.key=='0')  :
      self.xc = int(event.xdata+0.5)
      self.yc = int(event.ydata+0.5)
      if self.verbose>=2 :
        print('{:<10} {:>10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'oldvalue', 'value'))
      if self.verbose>=1 :         
        print('{:10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}'
              .format('mask',self.xc,self.yc,self.imamat[self.xc-1,self.yc-1],np.nan))
      if self.divvy :
        self.outfile.write('{:<30} {:10.0f} {:10.0f} {:10.2f} {:10.2f}\n'
                           .format(self.imagename,self.xc,self.yc,self.imamat[self.xc-1,self.yc-1],np.nan))
      else :
        self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10}\n'
                           .format('#','method','xc', 'yc', 'oldvalue', 'value'))
        self.outfile.write('{:<30} {:<10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}\n'
                           .format(self.imagename,'mask',self.xc,self.yc,self.imamat[self.xc-1,self.yc-1],np.nan))
      mx0 = self.xc-5
      mx1 = self.xc+5
      my0 = self.yc-5
      my1 = self.yc+5
      r2max = 5**2+1
      #self.imamat[mx0-1:mx1,my0-1:my1] = np.nan
      for ix in range(mx0,mx1+1) :
        for iy in range(my0,my1+1) :
            r2 = (ix-self.xc)**2+(iy-self.yc)**2
            if r2<=r2max :
              self.imamat[ix-1,iy-1] = np.nan 
      #print('imamat',self.imamat[mx0-1:mx1,my0-1:my1])
      #if self.verbose>=1 :         
      #  print('{:10} {:10.0f} {:10.0f} {:10.2f} {:10.2f}'.format('mask',self.xc,self.yc,
      #        self.imamat[self.xc-1,self.yc-1],np.nan))
      self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower', cmap=self.cmap,
                                    vmin=self.cut0, vmax=self.cut1, extent=self.extent)
      self.imshow.figure.canvas.draw_idle()
      #self.load()
      
    if (event.button==1 and event.key=='1')  :
      pos = [event.xdata, event.ydata]
      self.xc, self.yc = self.get_cursor(pos,verbose=self.verbose)
      if self.zoom!=0 :
        self.ax.xaxis.zoom(self.zoom) 
        self.ax.yaxis.zoom(self.zoom)
      self.imshow.figure.canvas.draw_idle()

    if (event.button==1 and event.key=='2') :
      pos = [event.xdata, event.ydata]
      self.xc, self.yc = self.fit_psf(pos,use_moffat=False,verbose=self.verbose)
      if self.zoom!=0 :
        self.ax.xaxis.zoom(self.zoom) 
        self.ax.yaxis.zoom(self.zoom)
      self.imshow.figure.canvas.draw_idle()
  
    if (event.button==1 and event.key=='3') :
      pos = [event.xdata, event.ydata]
      self.xc, self.yc = self.fit_psf(pos,use_moffat=True,verbose=self.verbose)
      if self.zoom!=0 :
        self.ax.xaxis.zoom(self.zoom) 
        self.ax.yaxis.zoom(self.zoom)
      self.imshow.figure.canvas.draw_idle()
  
    if (event.button==1 and event.key=='4') :
      pos = [event.xdata, event.ydata]
      self.xc, self.yc = self.fit_psf(pos,use_moffat=False,verbose=-1)
      ixc=int(self.xc+0.5)
      iyc=int(self.yc+0.5)
  
      ####### growing curve
      prec = self.precaper
      maxaper = self.maxaper
      # self.sky0 = np.nanmean(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50])
      # self.sky0 = np.nanmedian(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50])
      #c, low, upp = sigmaclip(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50],3.,2.,nan_policy=‘omit’)
      #cenfunc='np.nanmean',stdfunc='np.nanstd)')
      #c, low, upp = sigmaclip(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50],3.,2.)
      #self.sky0 = np.average(c)
      #self.sky0 = np.nanmean(c)
      self.sky0 = np.nanmedian(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50])
      #print("c, low, upp  = ",np.average(c), low, upp )
      #print(" self.sky0=", self.sky0)
      self.r = np.linspace(1.,maxaper,prec)
      #print(self.r)
      aper2 = self.r**2+1.
      self.grow = np.zeros(prec)
      self.area = np.zeros(prec)        
      self.growsky = np.zeros(prec)
      for i in range(prec) :
        s=0.
        n=0
        #print('{:10} {:10.0f} {:10.0f} {}'.format('test mask',ixc,iyc,self.imamat[ixc-1,iyc-1]))
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
      #print("self.grow=",self.grow)
      self.maxgrow = np.nanmax(self.grow)
      grow_min = np.nanmin(self.grow)
      grow_max = np.nanmax(self.grow)
      grow_del = grow_max - grow_min
      #print("grow_min=",grow_min)
      #print("grow_max=",grow_max)
      self.ylim0 = grow_min-0.1*grow_del
      self.ylim1 = grow_max+0.1*grow_del

      fig2, self.ax2 = plt.subplots(figsize=(6,4),facecolor='w')
      self.lineg, = self.ax2.plot(self.r,self.grow,color='lightblue',zorder=1)
      self.linef, = self.ax2.plot([None,None],[None,None],color='mistyrose',zorder=0)
      self.linec, = self.ax2.plot(self.r,self.grow,color='b',zorder=3)
      self.pointc,= self.ax2.plot(self.r,self.grow, marker='o', linestyle=None, ms=3, linewidth=0,color='b',zorder=4)
      self.line1, = self.ax2.plot([None,None],[None,None],color='grey',zorder=2)
      self.line2, = self.ax2.plot([None,None],[None,None],color='grey',zorder=2)
      self.lineh, = self.ax2.plot([None,None],[None,None],color='r',   zorder=2)
      #linef, = plt.plot([p1.x_fwhm,p1.x_fwhm],[0.,1.2*p1.totflux],color='g',zorder=0)
      self.ax2.set_xlabel("radius")
      self.ax2.set_ylabel("flux")
      self.axtitle = self.ax2.set_title('{} {:10.1f} {}'.format(self.imagename,self.totflux,'(growing curve)'))
      self.ax2.set_ylim([self.ylim0,self.ylim1])

      self.click2 = False
      cursor2 = Cursor( self.ax2, useblit=False, color='g', linewidth=1 )
      cid2 = fig2.canvas.mpl_connect('button_press_event',self.growing_curve)
      #cid3 = fig2.canvas.mpl_connect('close_event',self.growing_curve_close)
      cid4 = fig2.canvas.mpl_connect('key_press_event',self.keypress)
      if self.verbose>=2 :
        print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'totflux', 'sky', 'r0', 'r1', 'errtotflux' ))
      if not(self.divvy) :
        self.outfile.write('{:<30} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'.format('#','method','xc', 'yc', 'totflux', 'sky', 'r0', 'r1', 'errtotflux' ))

      fig2.show()

        
  def keypress(self, event):
      #print("event.key = ",event.key)
      #if (event.key=='z' or event.key=='Z' or event.key=='x' or event.key=='X' or event.key=='0') :              
      if (event.key=='z' or event.key=='x' or event.key=='Z' or event.key=='X') :              
        self.dcut = self.cut1 - self.cut0
        if (event.key=='z') :
          self.cut0 = self.cut1 - self.dcut*self.cut_zoom        
        if (event.key=='x') :
          self.cut0 = self.cut1 - self.dcut/self.cut_zoom        
        if (event.key=='Z') :
          self.cut1 = self.cut0 + self.dcut/self.cut_zoom
        if (event.key=='X') :
          self.cut1 = self.cut0 + self.dcut*self.cut_zoom        
        #if (event.key=='0') :
        #  self.cut0 = 0.
        print('cuts       {:<10.1f} {:<10.1f}'.format(self.cut0,self.cut1))
        #self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower',
        #                              cmap=self.cmap, vmin=self.cut0, vmax=self.cut1)
        self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower',
                                      cmap=self.cmap, vmin=self.cut0, vmax=self.cut1, extent=self.extent)
        self.imshow.figure.canvas.draw_idle()
      if event.key==' ' :
        if verbose>=1 :
          print('---------------------------------------------------------------------------------------------------------')
          if not(self.divvy) :
            self.outfile.write('---------------------------------------------------------------------------------------------------------\n')

        
      if event.key=='o' :
        xn, yn = np.linspace(0.5,self.nx+0.5, self.nx, dtype=float), np.linspace(0.5,self.ny+0.5, self.ny, dtype=float)
        Xn, Yn = np.meshgrid(xn,yn)
        Z = np.transpose(self.imamat)
        ima_min = np.nanmin(self.imamat)
        ima_max = np.nanmax(self.imamat)
        ima_del = ima_max - ima_min
        levels = np.linspace(ima_min+0.05*ima_del, ima_max-0.05*ima_del, 20)
        idc = self.ax.contour(Xn, Yn, Z, levels, colors='g')
        self.imshow.figure.canvas.draw_idle()


      
        
########################### MAIN

parser = argparse.ArgumentParser(description='showfits.py (2021-03-13) copyright Arno Riffeser (arri@usm.lmu.de)\n')
parser.add_argument(              'imagelist',       nargs='*',                   help='image list')
parser.add_argument('-figsize',   dest='figsize',    type=str,   default='8.75,8.75',   help='[%(default)s] figsize' )
parser.add_argument('-v',         dest='verbose',    type=int,   default=2,       help='[%(default)s] verbose' )
parser.add_argument('-precaper',  dest='precaper',   type=int,   default=80,      help='[%(default)s] precaper' )
parser.add_argument('-maxaper',   dest='maxaper',    type=int,   default=80,      help='[%(default)s] maxaper' )
parser.add_argument('-fitradius', dest='fitradius',  type=int,   default='7',     help='[%(default)s] fitradius' )
parser.add_argument('-divvy',     dest='divvy',      action='store_true',         help='[%(default)s] divvy' )
parser.add_argument('-zoom',      dest='zoom',       type=int,   default='0',     help='[%(default)s] zoom' )
args = parser.parse_args()

figsize   = np.array(args.figsize.split(','),dtype=float)
verbose   = args.verbose
fitradius = args.fitradius
precaper  = args.precaper
maxaper   = args.maxaper


if args.imagelist==[] :
    parser.print_help()
    print("\nusage in plot window:")
    print("    1+click - get cursor")
    print("    2+click - center gauss")
    print("    3+click - center moffat")
    print("    4+click - growing curve")
    print("    0+click - masking")
    print("    z - decrease lower cut")
    print("    x - increase lower cut")
    print("    Z - decrease higher cut")
    print("    X - increase higher cut")
    print("    p - pan (move center)")
    print("    r - reset (zoom home)")
    print("    f - full frame")
    print("    c - zoom back")
    print("    v - zoom forward")
    print("    left  - zoom back")
    print("    right - zoom forward")
    print("    k - x axis log")
    print("    l - y axis log")
    print("    s - save image")
    print("    q - quit")
    exit(0)


    
# https://matplotlib.org/3.3.0/users/interactive_guide.html

imagelist = args.imagelist
outfile  = open( 'showfits.tab' , 'w' )
fig = plt.figure(figsize=(figsize[0],figsize[1]),facecolor='w',dpi=100)


ia = Image_Analyzer(fitradius,precaper,maxaper,args.zoom,verbose,args.divvy,outfile)
ia.start(fig)
ia.storelist(imagelist)
ia.load()
ia.connect()
plt.show()
ia.disconnect()



outfile.close()
  
