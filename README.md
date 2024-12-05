# python

## redmine.py

- options

```
usage: redmine.py [-h] [-v LEVEL] [-p PROJID] [-i ISSUEID] [-l LAST] [-a ASSIGNED_TO] [-f]

redmine.py  2024-12-05  (Arno Riffeser, USM)

options:
  -h, --help      show this help message and exit
  -v LEVEL        [1] level
  -p PROJID       [0] projID
  -i ISSUEID      [0] issueID
  -l LAST         [2] nr last journals
  -a ASSIGNED_TO  [] assigned_to
  -f              [False] full list (also closed)
```

- examples:
```
  ./redmine.py -h
  ./redmine.py
  ./redmine.py -p 37
  ./redmine.py -p 37 -a "Arno Riffeser"
  ./redmine.py -p 95
  ./redmine.py -p 95 -f
  ./redmine.py -p 95 -v 3    -l 1    # VERY USEFUL 
  ./redmine.py -i 3617
  ./redmine.py -i 3617 -v 2  -l 1
  ./redmine.py -i 3617 -v 2  -l 0          
  ./redmine.py -i 3617 -v 3  -l 2    # VERY USEFUL 

  ./redmine.py -i 3617 -v 4
  ./redmine.py -i 3617 -v 5
  ./redmine.py -i 3617 -v 6

important:
  use a terminal width of min. 150

```



## showfits.py

- sequential fits viewer
  with some interactive analysis tools like gauss fit, moffat fit, growing curve analysis and masking

- options

```
     > showfits.py
     usage: showfits.py [-h] [-figsize FIGSIZE] [-v VERBOSE] [-precaper PRECAPER] [-maxaper MAXAPER] [-fitradius FITRADIUS]
                        [-zoom ZOOM] [-interpol INTERPOL] [-cmap CMAP] [-full] [-automask] [-autopng] [-autojpg] [-wcslim] [-fixcuts]
                        [-cuts CUTS]
                        [imagelist ...]
     
     showfits.py V1.4.11 (2024-12-05) (c) USM, Arno Riffeser (arri@usm.lmu.de) based on showfits of Claus Goessl and fitsedit of
     Johannes Koppenhoefer
     
     positional arguments:
       imagelist             image list
     
     options:
       -h, --help            show this help message and exit
       -figsize FIGSIZE      [10.0,8.5] figsize
       -v VERBOSE            [2] verbose
       -precaper PRECAPER    [80] precaper
       -maxaper MAXAPER      [80] maxaper
       -fitradius FITRADIUS  [15] fitradius
       -zoom ZOOM            [0] zoom
       -interpol INTERPOL    [1] interpolation
       -cmap CMAP            [102] cmap
       -full                 [False] full output
       -automask             [False] autosave all masks
       -autopng              [False] autosave all images as png
       -autojpg              [False] autosave all images as jpg
       -wcslim               [False] keep lim according WCS
       -fixcuts              [False] keep cuts fixedl
       -cuts CUTS            [0,0] cuts low and high
     
     usage in plot window:
               1 - get cursor    at mouse position
               2 - center gauss  at mouse position
               3 - center moffat at mouse position
               4 - growth curve  at mouse position
               5 - psf profile   at mouse position
               6 - star contour  at mouse position
               7 - mask     satellite with two mouse positions
               8 - mask inf satellite with two mouse positions
               0 - mask circle   at mouse position
               b - big mask circle   at mouse position
               m - save mask
               n - save mask
               j - save mask
               @ - UNDO last mask command
               z - decrease lower cut
               Z - increase lower cut
               x - decrease higher cut
               X - increase higher cut
             h r - reset (zoom home)
             < c - zoom back
             > v - zoom forward
               p - pan (move center)
               o - Zoom-to-rect
               s - save image
               f - full frame
               g - Toggle major grids
               G - Toggle minor grids
             k L - x axis log/linear
               l - y axis log/linear
               q - quit
```

- example

```
   > showfits.py  M51*.fits
   nr_images   3
   1   M51_r_stx_220324_103.fits
   2   M51_r_stx_220324_104.fits
   3   M51_r_stx_220324_105.fits
   =============================================================================================================
   nr          1  of  3
   imagename   M51_r_stx_220324_103.fits
   imagesize   3172       4096      
   cuts        1034.0     2192.0    
   cuts        1034.0     1999.0    
   cuts        1034.0     1838.2    
   cuts        1034.0     1704.1    
   cuts        1034.0     1592.4    
   cuts        1034.0     1499.4    
   cuts        1034.0     1421.8    
   cuts        1034.0     1357.2    
   cuts        1034.0     1303.3    
   #                  xc         yc      value
   pixel             564       3262    1027.00
   #                  xc         yc      value
   pixel             956       1979    3281.00
   #                  xc         yc    totflux        sky          A       sigx       sigy        phi       fwhm
   gauss          956.15    1978.64    57363.2   1065.910    2503.22       1.60       2.28      76.42       4.50
   #                  xc         yc    totflux        sky          A       sigx       sigy        phi       fwhm
   moffat         956.15    1978.64    67269.7   1055.837    2693.50       2.35       3.39      76.51       4.07
   #                  xc         yc    totflux        sky         r0         r1 errtotflux
   grow_curve     956.15    1978.64    63810.5   1060.331       24.0       70.0     1409.5
   grow_curve     956.15    1978.64    63590.8   1060.353       19.0       70.0     1126.3
   grow_curve     956.15    1978.64    63509.3   1060.368       19.0       77.0     1126.3
   ---------------------------------------------------------------------------------------
   grow_curve     956.15    1978.64    63550.2   1060.361       19.0       74.0     1126.3
   ---------------------------------------------------------------------------------------
   ---------------------------------------------------------------------------------------
```



## exposure_time_calculator.py

- simply change the mag_in and SN_in vector (vector because all filters are always calculated at the same time), then texp, Npeak, Nsky are calculated

- if you select tel="WST 2m defocussed", then an example for defocussed images is calculated...

- you can also specify texp_in and mag_in then SN will be calculated and so on...

- example:

```
   tel='WST 2m'
   SN_in   = 10. * np.ones(5)   # S/N in aperture
   mag_in  = []
   texp_in = 100.* np.ones(5)   # exposure time per image, sec
```

```
   ################################################################################
   telescope    = WST 2m
   filter       = [       u   ,       g   ,       r   ,       i   ,       z    ]
   --------------------------------------------------------------------------------
   SN_in        = [     10.00 ,     10.00 ,     10.00 ,     10.00 ,     10.00  ]
   texp_in      = [    100.00 ,    100.00 ,    100.00 ,    100.00 ,    100.00  ]
   --------------------------------------------------------------------------------
   number       = [      1    ,      1    ,      1    ,      1    ,      1     ]
   psf  ["]     = [      1.0  ,      1.0  ,      1.0  ,      1.0  ,      1.0   ]
   Aper ["]     = [      1.1  ,      1.1  ,      1.1  ,      1.1  ,      1.1   ]
   AM           = [      1.1  ,      1.1  ,      1.1  ,      1.1  ,      1.1   ]
   skymu_AB     = [     20.0  ,     20.0  ,     20.0  ,     20.0  ,     20.0   ]
   galmu_AB     = [    100.0  ,    100.0  ,    100.0  ,    100.0  ,    100.0   ]
   --------------------------------------------------------------------------------
   f_gauss      = [      0.57 ,      0.57 ,      0.57 ,      0.57 ,      0.57  ]
   RON          = [      2.6  ,      2.6  ,      2.6  ,      2.6  ,      2.6   ]
   q            = [     17.7% ,     36.7% ,     44.1% ,     27.1% ,     17.9%  ]
   dA           = [      2.70 ,      2.70 ,      2.70 ,      2.70 ,      2.70  ]
   ZP           = [     24.1  ,     25.4  ,     25.4  ,     24.7  ,     24.1   ]
   ================================================================================
   LIMITING MAGNITUDE:
      INPUT:
         SN     = [     10.0  ,     10.0  ,     10.0  ,     10.0  ,     10.0   ]
         texp   = [    100.0  ,    100.0  ,    100.0  ,    100.0  ,    100.0   ]
         number = [      1    ,      1    ,      1    ,      1    ,      1     ]
      RESULT:
         mag    = [     20.70 ,     21.82 ,     21.91 ,     21.53 ,     21.24  ]
         Npeak  = [     45.7  ,     79.8  ,     79.8  ,     57.7  ,     45.6   ]
         Nobj   = [   1125.8  ,   2301.9  ,   2300.6  ,   1662.3  ,   1315.7   ]
         Nsky   = [    197.1  ,    657.1  ,    656.3  ,    328.4  ,    196.6   ]
   ################################################################################
```
