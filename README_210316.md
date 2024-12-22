# usmpy

## showfits

- fits viewer with gauss fit, moffat fit and interactive growing curve analysis

- options

```
  > showfits.py

  usage: showfits.py [-h] [-figsize FIGSIZE] [-v VERBOSE] [-precaper PRECAPER]
                     [-maxaper MAXAPER] [-fitradius FITRADIUS] [-zoom ZOOM]
                     [-divvy]
                     [imagelist [imagelist ...]]
  
  showfits.py V1.0 (2021-03-16) (c) USM, Arno Riffeser (arri@usm.lmu.de)
  
  positional arguments:
    imagelist             image list
  
  optional arguments:
    -h, --help            show this help message and exit
    -figsize FIGSIZE      [8.75,8.75] figsize
    -v VERBOSE            [2] verbose
    -precaper PRECAPER    [80] precaper
    -maxaper MAXAPER      [80] maxaper
    -fitradius FITRADIUS  [15] fitradius
    -zoom ZOOM            [0] zoom
    -divvy                [False] divvy
  
  usage in plot window:
      1+click - get cursor
      2+click - center gauss
      3+click - center moffat
      4+click - growing curve
      0+click - masking
            z - decrease lower cut
            x - increase lower cut
            Z - decrease higher cut
            X - increase higher cut
            p - pan (move center)
            r - reset (zoom home)
            f - full frame
            c - zoom back
            v - zoom forward
         left - zoom back
        right - zoom forward
            k - x axis log
            l - y axis log
            s - save image
            q - quit
```

- example

```
   > showfits.py  *.fits
   nr_images   6
   1   17P_r_071029_095.fits
   2   17P_r_071029_095_copy.fits
   3   data.fits
   4   emtveb_dskq1wc-V151220-m31_f1_r_151220_230.fits
   5   res.fits
   6   test.fits
   =============================================================================================================
   nr          1  of  6
   imagename   17P_r_071029_095.fits
   imagesize   1100       1040      
   cuts        -6747.5    15317.8   
   cuts        -6747.5    19730.8   
   cuts        -6747.5    25026.5   
   cuts        -6747.5    31381.3   
   #                  xc         yc      value
   pixel             126        202   34883.00
   #                  xc         yc    totflux        sky          A       sigx       sigy        phi       fwhm
   gauss          125.89     201.97   300566.1    615.891   33148.02       1.31       1.10      59.48       2.83
   #                  xc         yc    totflux        sky          A       sigx       sigy        phi       fwhm
   moffat         125.88     201.96   348154.6    566.523   36104.29       1.91       1.61      58.88       2.53
   #                  xc         yc    totflux        sky         r0         r1 errtotflux
   grow.curve     125.89     201.97   334546.2    576.293       13.0       39.0      805.7
   grow.curve     125.89     201.97   342713.0    570.957       10.0       32.0      724.2
   ---------------------------------------------------------------------------------------
   grow.curve     125.89     201.97   342713.0    570.957       10.0       32.0      724.2
   ---------------------------------------------------------------------------------------
   ---------------------------------------------------------------------------------------
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
