# python

## rgbfits.py

- options
```
usage: rgbfits.py [-h] [-lupton] [-kmin KMIN] [-kmax KMAX] [-r R_FILENAME] [-g G_FILENAME] [-b B_FILENAME] [-rmin RMIN]
                  [-rmax RMAX] [-gmin GMIN] [-gmax GMAX] [-bmin BMIN] [-bmax BMAX] [-o REGION] [-log] [-jpg]

rgbfits.py  arri@usm.lmu.de  2024-12-22

options:
  -h, --help     show this help message and exit
  -lupton        [False] lupton
  -kmin KMIN     [-0.1] kmin
  -kmax KMAX     [1.0] kmax
  -r R_FILENAME  [] r file name
  -g G_FILENAME  [] g file name
  -b B_FILENAME  [] b file name
  -rmin RMIN     [None] rmin
  -rmax RMAX     [None] rmax
  -gmin GMIN     [None] gmin
  -gmax GMAX     [None] gmax
  -bmin BMIN     [None] bmin
  -bmax BMAX     [None] bmax
  -o REGION      [0:0,0:0] region
  -log           [False] log
  -jpg           [False] jpg
```

- examples

```
  ./rgbfits.py  -h
  ./rgbfits.py  -r M31_i.fits  -g M31_r.fits  -b M31_g.fits
  ./rgbfits.py  -r M31_i.fits  -g M31_r.fits  -b M31_g.fits -kmin -0.2 -kmax 0.5
  ./rgbfits.py  -r M31_i.fits  -g M31_r.fits  -b M31_g.fits -kmin -0.2 -kmax 0.5 -jpg
```




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
