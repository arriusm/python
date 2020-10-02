#!/Users/arri/miniconda3/bin/python
#/opt/anaconda3.7/bin/python
#

# (2020-10-02) by Arno Riffeser (arri@usm.lmu.de)


from math import *
import numpy as np
from   scipy.optimize import curve_fit,least_squares
import matplotlib.pyplot as plt
import matplotlib.mlab   as mlab
from   matplotlib.widgets import Cursor
from   scipy.optimize     import leastsq
from   scipy.stats        import sigmaclip
from   astropy.io         import fits as pyfits
import argparse

class Image_Analyzer:

  def __init__(self, ax, size=15, verbose=2, outfile=None):
      self.ax = ax
      self.ax2 = None
      self.cut_zoom = 1.2
      self.click2 = False
      self.sky = 0.
      self.linec = None
      self.line1 = None
      self.line2 = None
      self.lineh = None
      self.totflux = 0.
      self.sky = 0.
      self.rad0 = 0.
      self.rad1 = 0. 
      self.size = size
      self.verbose = verbose
      self.r = None
      self.grow = None
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
      rad = self.size//2
      x0 = xc-rad
      x1 = xc+rad
      y0 = yc-rad
      y1 = yc+rad
      sky = np.median(self.imamat)
      DATA = self.imamat[x0-1:x1,y0-1:y1]-sky
      bg = np.median(DATA)
      max = np.max(DATA)
      amp = max-bg
      cen = DATA[rad,rad]
      x, y = np.linspace(x0, x1, self.size, dtype=float), np.linspace(y0, y1, self.size, dtype=float)
      X, Y = np.meshgrid(y,x)
      # for ix in range(0,self.size) :
      #     for iy in range(0,self.size) :
      #         if DATA[ix,iy]==MD :
      #             xc = x[ix]
      #             yc = y[iy]
      XY = np.vstack((X.ravel(), Y.ravel()))
      Z = DATA.ravel()
      if not use_moffat : # gauss
        #astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp )
        astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp , bg )
        res = least_squares(self.residuals_gauss, astart, method='lm', args=(XY,DATA.ravel()), verbose=0)
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
          print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
        if verbose>=1 :
          print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}'.format('gauss',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
        if verbose>=0 :
          self.outfile.write('{:<20} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'.format('#','method','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
          self.outfile.write('{:<20} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'.format(self.imagename,'gauss',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
      else :
        #astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp )
        astart = ( xc , yc , 2.1 , 2. , 45./180.*np.pi , amp , bg )
        res = least_squares(self.residuals_moffat, astart, method='lm', args=(XY,DATA.ravel()), verbose=0) 
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
        if verbose>=0 :
          self.outfile.write('{:<20} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'.format('#','method','xc', 'yc', 'totflux' ,'sky', 'A', 'sigx', 'sigy', 'phi', 'fwhm' ))
          self.outfile.write('{:<20} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.2f} {:10.2f} {:10.2f} {:10.2f} {:10.2f}\n'.format(self.imagename,'moffat',afit[0],afit[1],totflux,afit[6],afit[5],abs(afit[2]),abs(afit[3]),afit[4]%np.pi/np.pi*180., fwhm))
      return afit[0],afit[1]

  def load(self,imamat,cut0=0.,cut1=1000.,imagename='') :
      self.title = imagename
      self.imagename = imagename
      self.imamat = imamat
      (self.nx,self.ny) = np.shape(self.imamat)
      self.cut0 = cut0
      self.cut1 = cut1
      self.dcut = cut1 - cut0
      if self.verbose>=2 :
        print('ima.size   {:10.0f} {:10.0f}'.format(self.nx,self.ny))
        print('cuts       {:10.1f} {:10.1f}'.format(self.cut0,self.cut1))
      self.extent = [0.5,self.nx+0.5,0.5,self.ny+0.5]
      self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower', cmap='gist_heat', vmin=self.cut0, vmax=self.cut1, extent=self.extent)
      #self.ax.set_xlim([1150,1170])
      #self.ax.set_ylim([870,890])
      self.imshow.figure.canvas.draw_idle()

  def connect(self):
      self.id1 = self.imshow.figure.canvas.mpl_connect('button_press_event',self.onclick)
      self.id2 = self.imshow.figure.canvas.mpl_connect('key_press_event',self.keypress)
      #self.id3 = self.imshow.figure.canvas.mpl_connect('close_event',self.close)
      self.imshow.figure.canvas.set_window_title(self.title)

  def disconnect(self):
      self.imshow.figure.canvas.mpl_disconnect(self.id1)
      self.imshow.figure.canvas.mpl_disconnect(self.id2)
      #self.imshow.figure.canvas.mpl_disconnect(self.id3)

  def close(self, event) :
      input("really close? <enter>")
      
  def keypress(self, event):
      #print("event.key = ",event.key)
      if (event.key=='z' or event.key=='Z' or event.key=='x' or event.key=='X' or event.key=='0') :              
        self.dcut = self.cut1 - self.cut0
        if (event.key=='z') :
          self.cut0 = self.cut1 - self.dcut*self.cut_zoom        
        if (event.key=='Z') :
          self.cut0 = self.cut1 - self.dcut/self.cut_zoom        
        if (event.key=='x') :
          self.cut1 = self.cut0 + self.dcut/self.cut_zoom
        if (event.key=='X') :
          self.cut1 = self.cut0 + self.dcut*self.cut_zoom        
        if (event.key=='0') :
          self.cut0 = 0.
        print('cuts       {:10.1f} {:10.1f}'.format(self.cut0,self.cut1))
        self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower', cmap='gist_heat', vmin=self.cut0, vmax=self.cut1)
        self.imshow = self.ax.imshow( np.transpose(self.imamat), origin='lower', cmap='gist_heat', vmin=self.cut0, vmax=self.cut1, extent=self.extent)
        self.imshow.figure.canvas.draw_idle()
      if event.key==' ' :
        if verbose>=1 :
          print('----------------------------------------------------------------------------------------------')
          self.outfile.write('----------------------------------------------------------------------------------------------\n')

        
      if event.key=='o' :
        xn, yn = np.linspace(0.5,self.nx+0.5, self.nx, dtype=float), np.linspace(0.5,self.ny+0.5, self.ny, dtype=float)
        Xn, Yn = np.meshgrid(xn,yn)
        Z = np.transpose(self.imamat)
        ima_min = np.min(self.imamat)
        ima_max = np.max(self.imamat)
        ima_del = ima_max - ima_min
        levels = np.linspace(ima_min+0.05*ima_del, ima_max-0.05*ima_del, 20)
        idc = self.ax.contour(Xn, Yn, Z, levels, colors='g')
        self.imshow.figure.canvas.draw_idle()

        
  def onclick(self, event):
      # print(event.button)
      # print('%s click: button=%d, key=%s, x=%d, y=%d, xdata=%f, ydata=%f' %
      #           ('double' if event.dblclick else 'single', event.button, event.key,
      #            event.x, event.y, event.xdata, event.ydata))  

      if (event.button==1 and event.key=='1')  :
        self.xc = int(event.xdata+0.5)
        self.yc = int(event.ydata+0.5)
        if self.verbose>=2 :
          print('{:<10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'value'))
        if self.verbose>=1 :         
          print('{:10} {:10.0f} {:10.0f} {:10.2f}'.format('pixel',self.xc,self.yc,self.imamat[self.xc-1,self.yc-1]))
        self.outfile.write('{:<20} {:<10} {:>10} {:>10} {:>10}\n'.format('#','method','xc', 'yc', 'value'))
        self.outfile.write('{:<20} {:<10} {:10.0f} {:10.0f} {:10.2f}\n'.format(self.imagename,'pixel',self.xc,self.yc,self.imamat[self.xc-1,self.yc-1]))
      # if (event.button==1 and event.key=='2') or event.dblclick :
      if (event.button==1 and event.key=='2') :
        pos = [event.xdata, event.ydata]
        self.xc, self.yc = self.fit_psf(pos,use_moffat=False,verbose=self.verbose)
  
      if (event.button==1 and event.key=='3') :
        pos = [event.xdata, event.ydata]
        self.xc, self.yc = self.fit_psf(pos,use_moffat=True,verbose=self.verbose)
  
      if (event.button==1 and event.key=='4') :
        pos = [event.xdata, event.ydata]
        self.xc, self.yc = self.fit_psf(pos,use_moffat=False,verbose=-1)
        ixc=int(self.xc+0.5)
        iyc=int(self.yc+0.5)
  
        ####### growing curve
        prec = 40
        maxaper = 80;
        # self.sky0 = np.average(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50])
        # self.sky0 = np.median(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50])
        c, low, upp = sigmaclip(self.imamat[ixc-1-50:ixc-1+50,iyc-1-50:iyc-1+50],3.,2.)
        self.sky0 = np.average(c)
        #print("c, low, upp  = ",np.average(c), low, upp )

        self.r = np.linspace(2.,maxaper,prec)
        #print(self.r)
        aper2 = self.r**2+1.
        self.grow = np.zeros(prec)
        self.area = np.zeros(prec)        
        for i in range(prec) :
          s=0.
          n=0
          for ix in range(ixc-maxaper-1,ixc+maxaper+1) :
            for iy in range(iyc-maxaper-1,iyc+maxaper+1) :
              r2 = (ix-self.xc)**2+(iy-self.yc)**2
              if r2<=aper2[i] :
                s += self.imamat[ix-1,iy-1] - self.sky0
                n += 1
          self.grow[i] = s
          self.area[i] = n
        self.maxgrow = np.max(self.grow)
        grow_min = np.min(self.grow)
        grow_max = np.max(self.grow)
        grow_del = grow_max - grow_min
        self.ylim0 = grow_min-0.1*grow_del
        self.ylim1 = grow_max+0.1*grow_del

        fig2, self.ax2 = plt.subplots(figsize=(6,4),facecolor='w')
        self.lineg, = self.ax2.plot(self.r,self.grow,color='lightblue',zorder=1)
        self.linef, = self.ax2.plot([None,None],[None,None],color='mistyrose',zorder=0)
        self.linec, = self.ax2.plot(self.r,self.grow,color='b',zorder=3)
        self.line1, = self.ax2.plot([None,None],[None,None],color='grey',zorder=2)
        self.line2, = self.ax2.plot([None,None],[None,None],color='grey',zorder=2)
        self.lineh, = self.ax2.plot([None,None],[None,None],color='r',   zorder=2)
        #linef, = plt.plot([p1.x_fwhm,p1.x_fwhm],[0.,1.2*p1.totflux],color='g',zorder=0)
        self.ax2.set_xlabel("radius")
        self.ax2.set_ylabel("flux")
        self.axtitle = self.ax2.set_title('{} {:10.1f} {}'.format(self.imagename,self.totflux,'(growing curve)'))
        self.ax2.set_ylim([self.ylim0,self.ylim1])

        self.click2 = False
        cursor = Cursor( self.ax2, useblit=False, color='g', linewidth=1 )
        cid2 = fig2.canvas.mpl_connect('button_press_event',self.growing_curve)
        #cid3 = fig2.canvas.mpl_connect('close_event',self.growing_curve_close)
        cid4 = fig2.canvas.mpl_connect('key_press_event',self.keypress)
        if self.verbose>=2 :
          print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'totflux', 'sky', 'r0', 'r1' ))
        self.outfile.write('{:<20} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'.format('#','method','xc', 'yc', 'totflux', 'sky', 'r0', 'r1' ))

        fig2.show()

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
        # fitlist =  opt.curve_fit(sky, rfit,growfit, x0, sigma)
        # A = np.vstack([rfit**2*np.pi,np.ones(len(rfit))]).T
        A = np.vstack([areafit,np.ones(len(rfit))]).T
        a2,a0 = np.linalg.lstsq(A, growfit,rcond=None)[0]
        rx = np.linspace(0.,100,100)
        self.linef.set_xdata(rx)
        self.linef.set_ydata(a2*rx**2*np.pi+a0)
        self.linef.figure.canvas.draw_idle()
        self.totflux = a0
        #print(growfit[0]-a2*areafit[0])
        self.sky = self.sky0 + a2
        new_grow = self.grow - a2*self.area
        self.linec.set_ydata(new_grow)
        self.lineh.set_xdata([self.r[0],self.r[-1]])
        self.lineh.set_ydata([a0,a0])
        new_grow_min = np.min(new_grow)
        new_grow_max = np.max(new_grow)
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
          print('----------------------------------------------------------------------------')
        if self.verbose>=1 :         
          print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f}'.format('grow.curve',self.xc,self.yc,self.totflux,self.sky,self.rad0,self.rad1))
        if (event.dblclick) :
          print('----------------------------------------------------------------------------')
          self.outfile.write('{:<20} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f}\n'.format(self.imagename,'grow.curve',self.xc,self.yc,self.totflux,self.sky,self.rad0,self.rad1))

  #def growing_curve_close(self, event) :
      #print('----------------------------------------------------------------------------')
      #if self.rad0 < self.rad1 :
      #  if self.verbose>=2 :
      #    print('{:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format('#','xc', 'yc', 'totflux', 'sky', 'r0', 'r1' ))
      #  if self.verbose>=1 :         
      #    print('{:10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f}'.format('grow.curve',self.xc,self.yc,self.totflux,self.sky,self.rad0,self.rad1))
      #  self.outfile.write('{:<20} {:<10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'.format('#','method','xc', 'yc', 'totflux', 'sky', 'r0', 'r1' ))
      #  self.outfile.write('{:<20} {:<10} {:10.2f} {:10.2f} {:10.1f} {:10.3f} {:10.1f} {:10.1f}\n'.format(self.imagename,'grow.curve',self.xc,self.yc,self.totflux,self.sky,self.rad0,self.rad1))

        
########################### MAIN

parser = argparse.ArgumentParser(description='showfits.py (2020-10-302) by Arno Riffeser (arri@usm.lmu.de)\n')
parser.add_argument(              'imagelist',       nargs='*',                   help='image list')
parser.add_argument('-figsize',   dest='figsize',    type=str,   default='7,7',   help='[%(default)s] figsize' )
parser.add_argument('-v',         dest='verbose',    type=int,   default=2,       help='[%(default)s] verbose' )
parser.add_argument('-box',       dest='boxsize',    type=int,   default='9',     help='[%(default)s] boxsize' )
args = parser.parse_args()

figsize = np.array(args.figsize.split(','),dtype=float)
verbose = args.verbose
boxsize = args.boxsize

print("boxsize=",boxsize)

if args.imagelist==[] :
    print("  usage:")
    print("    1+click - get cutsor")
    print("    2+click - center gauss")
    print("    3+click - center moffat")
    print("    4+click - groving curve")
    print("    r - reset (zoom home)")
    print("    p - pan (move center)")
    print("    s - save image")
    print("    c - zoom back")
    print("    v - zoom forward")
    print("    left  - zoom back")
    print("    right - zoom forward")
    print("    k - x axis log")
    print("    l - y axis log" )
    print("    q - quit")
    exit(-1)

    
# fig, ax = plt.subplots(figsize=(figsize[0],figsize[1]),facecolor='w')
# ia = Image_Analyzer(ax)
# plt.show(block=False)
# plt.ioff()
# print("HALLO")
# fig.canvas.draw()
# import time

# https://matplotlib.org/3.3.0/users/interactive_guide.html

outfile  = open( 'showfits.tab' , 'w' )

show=True
print(args.imagelist)
for imagename in args.imagelist :

  if verbose>=2 :
    print('=============================================================================================================')
  if verbose>=1 :
    print("imagename  = ", imagename)
    
  imamat = np.transpose(pyfits.getdata(imagename))
  ima_med = np.median(imamat)
  ima_std = np.std(imamat)
  ima_min = np.min(imamat)
  ima_max = np.max(imamat)
  ima_del = ima_max - ima_min
  #cut0 = ima_min - 0.1 * ima_del
  #cut1 = ima_max - 0.1 * ima_del
  cut0 = ima_med - 3. * ima_std
  cut1 = ima_med + 5. * ima_std
    
  fig, ax = plt.subplots(figsize=(figsize[0],figsize[1]),facecolor='w')
  #cursor = Cursor(ax, useblit=False, color='g', linewidth=1 )
  ia = Image_Analyzer(ax,boxsize,verbose,outfile)
  ia.load(imamat, cut0, cut1, imagename)
  ia.connect()
  plt.show()
  ia.disconnect()

outfile.close()
  
