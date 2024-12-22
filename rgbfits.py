#!/Users/arri/miniconda3/bin/python
###!/opt/miniconda3/bin/python
###!/Users/arri/miniconda3/bin/python

import matplotlib as mpl
mpl.use('TkAgg')   
#mpl.use('MacOSX')

import numpy as np
import matplotlib.pyplot as plt
from   astropy.io            import fits 
#from   astropy.visualization import make_lupton_rgb
import astropy.visualization as av

import argparse 
# import scipy.misc
# import textwrap

# from   matplotlib.widgets import Slider, Button
from   matplotlib.widgets import RangeSlider, Button

col_1 = np.array([255.,0.,0.]   )  # R
col_2 = np.array([0.,255.,0.]   )  # G
col_3 = np.array([0.,0.,255.]   )  # B

mycol1 = tuple(col_1/255.)
mycol2 = tuple(col_2/255.)
mycol3 = tuple(col_3/255.)



np.seterr(all='ignore')

help = '''\

examples:
  ./rgbfits.py  -h
  ./rgbfits.py  -r M31_i.fits  -g M31_r.fits  -b M31_g.fits
  ./rgbfits.py  -r M31_i.fits  -g M31_r.fits  -b M31_g.fits -kmin -0.2 -kmax 0.5
  ./rgbfits.py  -r M31_i.fits  -g M31_r.fits  -b M31_g.fits -kmin -0.2 -kmax 0.5 -jpg
'''

parser = argparse.ArgumentParser(epilog=help, formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='rgbfits.py  arri@usm.lmu.de  2024-12-22\n')
parser.add_argument('-lupton', dest='lupton',  action="store_true", default=False,     help='[%(default)s] lupton' )
parser.add_argument('-kmin',   dest='kmin',             type=float, default=-0.1,      help='[%(default)s] kmin' )
parser.add_argument('-kmax',   dest='kmax',             type=float, default=1.0,       help='[%(default)s] kmax' )
parser.add_argument('-r',      dest='r_filename',       type=str,   default='',        help='[%(default)s] r file name' )
parser.add_argument('-g',      dest='g_filename',       type=str,   default='',        help='[%(default)s] g file name' )
parser.add_argument('-b',      dest='b_filename',       type=str,   default='',        help='[%(default)s] b file name' )
parser.add_argument('-rmin',   dest='rmin',             type=float, default=None,      help='[%(default)s] rmin' )
parser.add_argument('-rmax',   dest='rmax',             type=float, default=None,      help='[%(default)s] rmax' )
parser.add_argument('-gmin',   dest='gmin',             type=float, default=None,      help='[%(default)s] gmin' )
parser.add_argument('-gmax',   dest='gmax',             type=float, default=None,      help='[%(default)s] gmax' )
parser.add_argument('-bmin',   dest='bmin',             type=float, default=None,      help='[%(default)s] bmin' )
parser.add_argument('-bmax',   dest='bmax',             type=float, default=None,      help='[%(default)s] bmax' )
parser.add_argument('-o',      dest='region',           type=str,   default='0:0,0:0', help='[%(default)s] region' )
parser.add_argument('-log',    dest='log',     action="store_true", default=False,     help='[%(default)s] log' )
parser.add_argument('-jpg',    dest='jpg',     action="store_true", default=False,     help='[%(default)s] jpg' )
#parser.add_argument('-png',    dest='png',     action="store_true", default=False,     help='[%(default)s] png' )
args = parser.parse_args()


if args.r_filename=='' or args.g_filename=='' or args.b_filename=='' :
  parser.print_help()
  print('\n  no images, use :  ./rgbfits.py  -r r.fits  -g g.fits   -b b.fits')
  exit(0)




region = args.region

xreg=slice(0,-1)
yreg=slice(0,-1)
if region!='0:0,0:0' :
  rr = region.split(',')
  if len(rr)==2 :
    xx = rr[0].split(':')
    yy = rr[1].split(':')
    if len(xx)==2 and len(yy)==2 :
      x0 = int(xx[0])
      x1 = int(xx[1])
      y0 = int(yy[0])
      y1 = int(yy[1])
    else :
      exit(-1)
  else :
    exit(-1)
  if x0>x1 or y0>y1 :
    exit(-1)
  xreg=slice(x0-1,x1-1)
  yreg=slice(y0-1,y1-1)
  print(xreg,yreg)


  
#r_ filename = 's.tfdbo-IC434_H_qhy_221212_130.fits'
hdul = fits.open(args.r_filename)
c1_mat = hdul[0].data

#g_filename = 's.tfdbo-IC434_S_qhy_221212_136.fits'
hdul = fits.open(args.g_filename)
c2_mat = hdul[0].data

#b_filename = 's.tfdbo-IC434_O_qhy_221212_133.fits'
hdul = fits.open(args.b_filename)
c3_mat = hdul[0].data


c1 = c1_mat[xreg,yreg]
c2 = c2_mat[xreg,yreg]
c3 = c3_mat[xreg,yreg]

c1[np.isinf(c1)]=np.nan
c2[np.isinf(c2)]=np.nan
c3[np.isinf(c3)]=np.nan

c1mean = np.nanmean(c1)
c2mean = np.nanmean(c2)
c3mean = np.nanmean(c3)

c1std = np.nanstd(c1)
c2std = np.nanstd(c2)
c3std = np.nanstd(c3)

#kmin = 0.5
kmin  = args.kmin
if args.rmin == None :
  c1min = c1mean + kmin * c1std
else :
  c1min = args.rmin
if args.gmin == None :
  c2min = c2mean + kmin * c2std
else :
  c2min = args.gmin
if args.bmin == None :
  c3min = c3mean + kmin * c3std
else :
  c3min = args.bmin

# kmax = 1.
kmax  = args.kmax
if args.rmax == None :
  c1max = c1mean + kmax * c1std
else :
  c1max = args.rmax
if args.gmax == None :
  c2max = c2mean + kmax * c2std
else :
  c2max = args.gmax
if args.bmax == None :
  c3max = c3mean + kmax * c3std
else :
  c3max = args.bmax

print('color                 : {:>5s} {:>5s}'.format('min','max'))
print('red                   : {:5.1f} {:5.1f}'.format(c1min,c1max))
print('green                 : {:5.1f} {:5.1f}'.format(c2min,c2max))
print('blue                  : {:5.1f} {:5.1f}'.format(c3min,c3max))


if args.lupton :
  image_rgb = av.make_lupton_rgb(c1_mat, c2_mat, c3_mat,
                                 minimum=[c1min,c2min,c3min], Q=10, stretch=0.5)
  fig = plt.figure(figsize=(10,6.5),facecolor='w',dpi=120)
  ax  = fig.add_axes([0.,0.10,1.0,0.9])
  ax.axis('off')
  imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
  plt.savefig('lupton.png')  
  plt.show()
   

# linear
# def normalize(x, vmin=0, vmax=5000) :
#   y = (x-vmin) / (vmax-vmin)
#   return np.where(y<0.,0.,np.where(y>1.,1.,y))

# def normalize(x, vmin=0, vmax=5000) :
#   y = np.where(x-vmin<0.,0.,np.log10(1+x-vmin)/np.log10(1.+vmax-vmin) )
#   return np.where(y>1.,1.,y)

# from   matplotlib.colors import LogNorm
# see https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.Normalize.html
#     https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.AsinhNorm.html
#     https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.LogNorm.html

# from astropy.visualization import LogStretch.html
# see https://docs.astropy.org/en/stable/api/astropy.visualization.PowerDistStretch.html
#     https://docs.astropy.org/en/stable/api/astropy.visualization.PowerStretch.html
#     https://docs.astropy.org/en/stable/api/astropy.visualization.SinhStretch.html
#     https://docs.astropy.org/en/stable/api/astropy.visualization.SqrtStretch.html
#     https://docs.astropy.org/en/stable/api/astropy.visualization.SquaredStretch.html
#     https://docs.astropy.org/en/stable/api/astropy.visualization.LogStretch.html


# http://ds9.si.edu/doc/ref/how.html
def normalize(x, vmin=0, vmax=5000) :
  y = (x-vmin) / (vmax-vmin)
  y01 = np.where(y<0.,0.,np.where(y>1.,1.,y))
  if args.log :
    a = 1000.
    #return np.log(a*y01+1)/np.log(a)
    return np.log(1.+a*y01)/np.log(1.+a)
  else :
    return y01



class myParams :
  def __init__(self) :
    self.ch1_0 = 0.
    self.ch1_1 = 0.
    self.ch2_0 = 0.
    self.ch2_1 = 0.
    self.ch3_0 = 0.
    self.ch3_1 = 0.
    #self.fak1 = 1.
    #self.fak2 = 1.
    #self.fak3 = 1.
    self.jpg = args.jpg
    #self.png = args.png
    self.png = True

par    = myParams()
par.ch1_0 = c1min
par.ch1_1 = c1max
par.ch2_0 = c2min
par.ch2_1 = c2max
par.ch3_0 = c3min
par.ch3_1 = c3max



def calc_image_rgb(par) :

  global c1_mat,c2_mat,c3_mat
  
  '''
  fak_c1 = par.fak1  
  fak_c2 = par.fak2  
  fak_c3 = par.fak3  

  fak_norm_R = ( fak_c1*col_1[0] + fak_c2*col_2[0] + fak_c3*col_3[0] ) / 255.
  fak_norm_G = ( fak_c1*col_1[1] + fak_c2*col_2[1] + fak_c3*col_3[1] ) / 255.
  fak_norm_B = ( fak_c1*col_1[2] + fak_c2*col_2[2] + fak_c3*col_3[2] ) / 255.
  
  # print('fak_norm_r =',fak_norm_R) 
  # print('fak_norm_g =',fak_norm_G) 
  # print('fak_norm_b =',fak_norm_B) 
  
  R_norm = ( ( fak_c1*normalize(c1_mat,par.ch1_0,par.ch1_1)*col_1[0] +
               fak_c2*normalize(c2_mat,par.ch2_0,par.ch2_1)*col_2[0] +
               fak_c3*normalize(c3_mat,par.ch3_0,par.ch3_1)*col_3[0] ) / fak_norm_R
            ).astype(np.uint8)
  
  G_norm = ( ( fak_c1*normalize(c1_mat,par.ch1_0,par.ch1_1)*col_1[1] +
               fak_c2*normalize(c2_mat,par.ch2_0,par.ch2_1)*col_2[1] +
               fak_c3*normalize(c3_mat,par.ch3_0,par.ch3_1)*col_3[1] ) / fak_norm_G
            ).astype(np.uint8)
  
  B_norm = ( ( fak_c1*normalize(c1_mat,par.ch1_0,par.ch1_1)*col_1[2] +
               fak_c2*normalize(c2_mat,par.ch2_0,par.ch2_1)*col_2[2] +
               fak_c3*normalize(c3_mat,par.ch3_0,par.ch3_1)*col_3[2] ) / fak_norm_B
             ).astype(np.uint8)
  '''

  R_norm = (255.*normalize(c1_mat,par.ch1_0,par.ch1_1)).astype(np.uint8)  
  G_norm = (255.*normalize(c2_mat,par.ch2_0,par.ch2_1)).astype(np.uint8)
  B_norm = (255.*normalize(c3_mat,par.ch3_0,par.ch3_1)).astype(np.uint8)

  #print('------------------------------------------------------')
  #print('min(R) = ',np.min(R_norm),'                 max(R) = ',np.max(R_norm))
  #print('min(G) = ',np.min(G_norm),'                 max(G) = ',np.max(G_norm))
  #print('min(B) = ',np.min(B_norm),'                 max(B) = ',np.max(B_norm))
  
  return np.dstack((R_norm,G_norm,B_norm))



####### MAIN #######

image_rgb = calc_image_rgb(par)
  
fig = plt.figure(figsize=(10,6.5),facecolor='w',dpi=120)
ax  = fig.add_axes([0.,0.10,1.0,0.9])

#interpolation: None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16',
#               'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
#               'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
# ax.imshow(image_rgb, origin='lower', interpolation='nearest')
#imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
imax = ax.imshow(image_rgb, origin='lower', interpolation='none')

ax.axis('off')


'''

####  Slider
      
def slider_update_fak1(val):
    global par, ax, fig
    #print(val)
    if val != par.fak1 :
      par.fak1 = val
      print('fak1 = ',par.fak1)
      slider_fak1.set_val(val)
      image_rgb = calc_image_rgb(par)
      imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      fig.canvas.draw_idle()

def slider_update_fak2(val):
    global par, ax, fig
    #print(val)
    if val != par.fak2 :
      par.fak2 = val
      print('fak2 = ',par.fak2)
      slider_fak2.set_val(val)
      image_rgb = calc_image_rgb(par)
      imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      fig.canvas.draw_idle()

def slider_update_fak3(val):
    global par, ax, fig
    #print(val)
    if val != par.fak3 :
      par.fak3 = val
      print('fak3 = ',par.fak3)
      slider_fak3.set_val(val)
      image_rgb = calc_image_rgb(par)
      imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      fig.canvas.draw_idle()


def slider_update_ch1_0(val):
    global par, ax, fig
    #print(val)
    if val != par.ch1_0 :
      par.ch1_0 = val
      print('cuts1_0 = ',par.ch1_0)
      slider_ch1_0.set_val(val)
      image_rgb = calc_image_rgb(par)
      imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      fig.canvas.draw_idle()

def slider_update_ch2_0(val):
    global par, ax, fig
    #print(val)
    if val != par.ch2_0 :
      par.ch2_0 = val
      print('cuts2_0 = ',par.ch2_0)
      slider_ch2_0.set_val(val)
      image_rgb = calc_image_rgb(par)
      imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      fig.canvas.draw_idle()

def slider_update_ch3_0(val):
    global par, ax, fig
    #print(val)
    if val != par.ch3_0 :
      par.ch3_0 = val
      print('cuts3_0 = ',par.ch3_0)
      slider_ch3_0.set_val(val)
      image_rgb = calc_image_rgb(par)
      imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      fig.canvas.draw_idle()


def slider_update_ch1_1(val):
    global par, ax, fig
    #print(val)
    if val != par.ch1_1 :
      par.ch1_1 = val
      print('cuts1_1 = ',par.ch1_1)
      slider_ch1_1.set_val(val)
      image_rgb = calc_image_rgb(par)
      imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      fig.canvas.draw_idle()

def slider_update_ch2_1(val):
    global par, ax, fig
    #print(val)
    if val != par.ch2_1 :
      par.ch2_1 = val
      print('cuts2_1 = ',par.ch2_1)
      slider_ch2_1.set_val(val)
      image_rgb = calc_image_rgb(par)
      imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      fig.canvas.draw_idle()

def slider_update_ch3_1(val):
    global par, ax, fig
    #print(val)
    if val != par.ch3_1 :
      par.ch3_1 = val
      print('cuts3_1 = ',par.ch3_1)
      slider_ch3_1.set_val(val)
      image_rgb = calc_image_rgb(par)
      imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      fig.canvas.draw_idle()

'''


####  RangeSlider

def slider_update_ch1(val):
    global par, ax, fig
    #print(val)
    if val != (par.ch1_0,par.ch1_1) :
      par.ch1_0,par.ch1_1 = val
      print('cuts1 = ',par.ch1_0,par.ch1_1)
      slider_ch1.set_val(val)
      image_rgb = calc_image_rgb(par)
      #imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      imax = ax.imshow(image_rgb, origin='lower')
      fig.canvas.draw_idle()

def slider_update_ch2(val):
    global par, ax, fig
    #print(val)
    if val != (par.ch2_0,par.ch2_1) :
      par.ch2_0,par.ch2_1 = val
      print('cuts2 = ',par.ch2_0,par.ch2_1)
      slider_ch2.set_val(val)
      image_rgb = calc_image_rgb(par)
      #imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      imax = ax.imshow(image_rgb, origin='lower')
      fig.canvas.draw_idle()

def slider_update_ch3(val):
    global par, ax, fig
    #print(val)
    if val != (par.ch3_0,par.ch3_1) :
      par.ch3_0,par.ch3_1 = val
      print('cuts3 = ',par.ch3_0,par.ch3_1)
      slider_ch3.set_val(val)
      image_rgb = calc_image_rgb(par)
      #imax = ax.imshow(image_rgb, origin='lower', interpolation='gaussian')
      imax = ax.imshow(image_rgb, origin='lower')
      fig.canvas.draw_idle()



      
axcolor = 'grey'
startx = 0.03
starty = -0.09
Dfin = 0.03


'''

#  Slider

ax_fak1 = plt.axes([0.8,     starty+0.16, 0.15,      0.02], facecolor=axcolor)
ax_fak2 = plt.axes([0.8,     starty+0.13, 0.15,      0.02], facecolor=axcolor)
ax_fak3 = plt.axes([0.8,     starty+0.10, 0.15,      0.02], facecolor=axcolor)
slider_fak1 = Slider(ax_fak1, 'fak1', valmin=0., valmax=2., valfmt=" %1.2f", valinit=1., color=mycol1, dragging=True)
slider_fak2 = Slider(ax_fak2, 'fak2', valmin=0., valmax=2., valfmt=" %1.2f", valinit=1., color=mycol2, dragging=True)
slider_fak3 = Slider(ax_fak3, 'fak3', valmin=0., valmax=2., valfmt=" %1.2f", valinit=1., color=mycol3, dragging=True)
slider_fak1.label.set_size(8)
slider_fak2.label.set_size(8)
slider_fak3.label.set_size(8)
slider_fak1.on_changed(slider_update_fak1)
slider_fak2.on_changed(slider_update_fak2)
slider_fak3.on_changed(slider_update_fak3)

lenx = 0.30
ax_ch1_0 = plt.axes([startx,     starty+0.16, lenx,      0.02], facecolor=axcolor)
ax_ch2_0 = plt.axes([startx,     starty+0.13, lenx,      0.02], facecolor=axcolor)
ax_ch3_0 = plt.axes([startx,     starty+0.10, lenx,      0.02], facecolor=axcolor)
slider_ch1_0 = Slider(ax_ch1_0, r'R', valmin=c1min-1.*c1std, valmax=c1max+1.*c1std, valfmt="%1.1f", valinit=c1min, color=mycol1, dragging=False)
slider_ch2_0 = Slider(ax_ch2_0, r'G', valmin=c2min-1.*c2std, valmax=c2max+1.*c2std, valfmt="%1.1f", valinit=c2min, color=mycol2, dragging=False)
slider_ch3_0 = Slider(ax_ch3_0, r'B', valmin=c3min-1.*c3std, valmax=c3max+1.*c3std, valfmt="%1.1f", valinit=c3min, color=mycol3, dragging=False)
slider_ch1_0.label.set_size(8)
slider_ch2_0.label.set_size(8)
slider_ch3_0.label.set_size(8)
slider_ch1_0.on_changed(slider_update_ch1_0)
slider_ch2_0.on_changed(slider_update_ch2_0)
slider_ch3_0.on_changed(slider_update_ch3_0)

ax_ch1_1 = plt.axes([startx+0.37,     starty+0.16, lenx,      0.02], facecolor=axcolor)
ax_ch2_1 = plt.axes([startx+0.37,     starty+0.13, lenx,      0.02], facecolor=axcolor)
ax_ch3_1 = plt.axes([startx+0.37,     starty+0.10, lenx,      0.02], facecolor=axcolor)
slider_ch1_1 = Slider(ax_ch1_1, r'R', valmin=c1min+0.*c1std, valmax=c1max+10.*c1std, valfmt="%1.1f", valinit=c1max, color=mycol1, dragging=False)
slider_ch2_1 = Slider(ax_ch2_1, r'G', valmin=c2min+0.*c2std, valmax=c2max+10.*c2std, valfmt="%1.1f", valinit=c2max, color=mycol2, dragging=False)
slider_ch3_1 = Slider(ax_ch3_1, r'B', valmin=c3min+0.*c3std, valmax=c3max+10.*c3std, valfmt="%1.1f", valinit=c3max, color=mycol3, dragging=False)
slider_ch1_1.label.set_size(8)
slider_ch2_1.label.set_size(8)
slider_ch3_1.label.set_size(8)
slider_ch1_1.on_changed(slider_update_ch1_1)
slider_ch2_1.on_changed(slider_update_ch2_1)
slider_ch3_1.on_changed(slider_update_ch3_1)
'''

# RangeSlider

lenx = 0.85
ax_ch1 = plt.axes([startx,     starty+0.16, lenx,      0.02], facecolor=axcolor)
ax_ch2 = plt.axes([startx,     starty+0.13, lenx,      0.02], facecolor=axcolor)
ax_ch3 = plt.axes([startx,     starty+0.10, lenx,      0.02], facecolor=axcolor)
slider_ch1 = RangeSlider(ax_ch1, r'R', valmin=c1min-1.*c1std, valmax=c1max+10.*c1std, valfmt="%1.1f", valinit=[c1min,c1max], color=mycol1, dragging=False)
slider_ch2 = RangeSlider(ax_ch2, r'G', valmin=c2min-1.*c2std, valmax=c2max+10.*c2std, valfmt="%1.1f", valinit=[c2min,c2max], color=mycol2, dragging=False)
slider_ch3 = RangeSlider(ax_ch3, r'B', valmin=c3min-1.*c3std, valmax=c3max+10.*c3std, valfmt="%1.1f", valinit=[c3min,c3max], color=mycol3, dragging=False)
slider_ch1.label.set_size(8)
slider_ch2.label.set_size(8)
slider_ch3.label.set_size(8)
slider_ch1.on_changed(slider_update_ch1)
slider_ch2.on_changed(slider_update_ch2)
slider_ch3.on_changed(slider_update_ch3)





def button_update_save(val):
    global par, ax
    image_rgb = calc_image_rgb(par)
    if par.jpg :
      print('saving image as rgbfits.jpg')
      plt.imsave('rgbfits.jpg', image_rgb, origin='lower',dpi=100)
    elif par.png :
      print('saving image as rgbfits.png')
      plt.imsave('rgbfits.png', image_rgb, origin='lower',dpi=100)


ax_save = plt.axes([0.0, 0.97, 0.05, 0.03], facecolor=axcolor)
button_save = Button(ax_save, 'save')
button_save.label.set_size(8)
button_save.on_clicked(button_update_save)



#plt.savefig('rgbfits.png')
plt.show()
