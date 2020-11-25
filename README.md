# usmpy


## showfits

- fits viewer with gauss fit, moffat fit and interactive growing curve analysis

- options

```
  > showfits.py -h

  usage: showfits.py [-h] [-figsize FIGSIZE] [-v VERBOSE] [-precaper PRECAPER]
                     [-maxaper MAXAPER] [-fitradius FITRADIUS]
                     [imagelist [imagelist ...]]
  
  showfits.py (2020-11-18) by Arno Riffeser (arri@usm.lmu.de)
  
  positional arguments:
    imagelist             image list
  
  optional arguments:
    -h, --help            show this help message and exit
    -figsize FIGSIZE      [8,6] figsize
    -v VERBOSE            [2] verbose
    -precaper PRECAPER    [80] precaper
    -maxaper MAXAPER      [80] maxaper
    -fitradius FITRADIUS  [7] fitradius
```

- runtime commands

```
  > showfits.py 
  usage:
    1+click - get cursor
    2+click - center gauss
    3+click - center moffat
    4+click - growing curve
    0+click - masking
    r - reset (zoom home)
    p - pan (move center)
    s - save image
    c - zoom back
    v - zoom forward
    left  - zoom back
    right - zoom forward
    k - x axis log
    l - y axis log
    q - quit
```

- example

```
  > showfits.py  -figsize 7,5  -v 2  -fitradius 7  HD12303_r.fit
  ['HD12303_r.fit']
  =============================================================================================================
  imagename      HD12303_r.fit
  image(nx,ny)        2328       1760
  cuts               611.0      816.3
  #
  #                  xc         yc      value
  pixel            1299        982     664.00
  #                  xc         yc    totflux        sky          A       sigx       sigy        phi       fwhm
  gauss         1157.72     879.03    40930.1    752.097    4807.42       1.22       1.11      20.35       2.74
  #                  xc         yc    totflux        sky          A       sigx       sigy        phi       fwhm
  moffat        1157.72     879.03    48267.8    721.033    5292.47       1.80       1.61      16.68       2.46
  #                  xc         yc    totflux        sky         r0         r1 errtotflux
  grow.curve    1157.72     879.03    62003.1    680.147       22.0       55.0     1046.7
  grow.curve    1157.72     879.03    61948.3    680.161       22.0       57.0     1046.7
  ---------------------------------------------------------------------------------------
  grow.curve    1157.72     879.03    61948.3    680.161       22.0       57.0     1046.7
```

## exposure time calculator

- simply change the mag_in and SN_in vector (vector because all filters are always calculated at the same time), then texp, Npeak, Nsky are calculated

- if you select tel="WST 2m defocussed", then an example for defocussed images is calculated...

- you can also specify texp_in and mag_in then SN will be calculated and so on...

- example:

```
   SN_in   = 200. * np.ones(5)            # S/N in aperture
   mag_in  = 10. * np.ones(5)             # object AB [AB mag]
   texp_in = [5.79, 1.09, 0.94, 1.7, 0]   # exposure time in sec 
```

```
   ################################################################################
   telescope    = WST 40cm
   filter       = [       u   ,      g   ,      r   ,      i   ,      z    ]
   --------------------------------------------------------------------------------
   mag_in       = [     10.00 ,    10.00 ,    10.00 ,    10.00 ,    10.00  ]
   SN_in        = [    200.00 ,   200.00 ,   200.00 ,   200.00 ,   200.00  ]
   texp_in      = [      5.79 ,     1.09 ,     0.94 ,     1.70 ,     0.00  ]
   --------------------------------------------------------------------------------
   number       = [      1    ,     1    ,     1    ,     1    ,     1     ]
   psf  ["]     = [      1.5  ,     1.5  ,     1.5  ,     1.5  ,     1.5   ]
   Aper ["]     = [      2.0  ,     2.0  ,     2.0  ,     2.0  ,     2.0   ]
   AM           = [      1.0  ,     1.0  ,     1.0  ,     1.0  ,     1.0   ]
   skymu_AB     = [     20.0  ,    20.0  ,    20.0  ,    20.0  ,    20.0   ]
   galmu_AB     = [    100.0  ,   100.0  ,   100.0  ,   100.0  ,   100.0   ]
   --------------------------------------------------------------------------------
   f_gauss      = [      0.7  ,     0.7  ,     0.7  ,     0.7  ,     0.7   ]
   RON          = [      5.0  ,     5.0  ,     5.0  ,     5.0  ,     5.0   ]
   q            = [     15.0% ,    35.0% ,    45.0% ,    30.0% ,     0.0%  ]
   dA           = [      0.12 ,     0.12 ,     0.12 ,     0.12 ,     0.12  ]
   ZP           = [     20.5  ,    22.0  ,    22.0  ,    21.4  ,    -inf   ]
   ================================================================================
   (S/N):
      INPUT:
         mag    = [     10.0  ,    10.0  ,    10.0  ,    10.0  ,    10.0   ]
         texp   = [      5.8  ,     1.1  ,     0.9  ,     1.7  ,     0.0   ]
         number = [      1    ,     1    ,     1    ,     1    ,     1     ]
      RESULT:
         SN     = [    200.05 ,   200.31 ,   199.52 ,   199.51 ,     0.00  ]
         Npeak  = [   8177.4  ,  8196.3  ,  8131.6  ,  8131.2  ,     0.0   ]
         Nsky   = [      3.8  ,     2.7  ,     2.5  ,     2.4  ,     0.0   ]
   ================================================================================
   LIMITING MAGNITUDE:
      INPUT:
         SN     = [    200.0  ,   200.0  ,   200.0  ,   200.0  ,   200.0   ]
         texp   = [      5.8  ,     1.1  ,     0.9  ,     1.7  ,     0.0   ]
         number = [      1    ,     1    ,     1    ,     1    ,     1     ]
      RESULT:
         mag    = [     10.00 ,    10.00 ,     9.99 ,     9.99 ,     -inf  ]
         Npeak  = [   8173.0  ,  8171.2  ,  8170.9  ,  8170.8  ,     nan   ]
         Nsky   = [      3.8  ,     2.7  ,     2.5  ,     2.4  ,     0.0   ]
   ================================================================================
   EXPOSURE TIME:
      INPUT:
         mag    = [     10.0  ,    10.0  ,    10.0  ,    10.0  ,    10.0   ]
         SN     = [    200.0  ,   200.0  ,   200.0  ,   200.0  ,   200.0   ]
         number = [      1    ,     1    ,     1    ,     1    ,     1     ]
      RESULT:
         texp   = [      5.79 ,     1.09 ,     0.94 ,     1.70 ,      nan  ]
         Npeak  = [   8173.0  ,  8171.2  ,  8170.9  ,  8170.8  ,     nan   ]
         Nsky   = [      3.8  ,     2.7  ,     2.5  ,     2.4  ,     nan   ]

```
