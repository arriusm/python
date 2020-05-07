#!/Users/arri/anaconda3/bin/python

##### version 2020-05-07 Arno Riffeser arri@usm.lmu.de #####

import numpy as np
np.seterr(all='ignore')

#  all calculated using AB-magnitudes:   AB_nu  = -2.5 lg (f_nu) - 48.60            with f_nu in ergs/cm**2/s/Hz
#                                  or:   f_nu/Jy  = 3631* dex(-0.4 AB_nu)       (see e.g. Fukugita et al AJ 111, 1748 (1996)) 
#   
#  dN_nu     =  f_nu/(h nu) dnu dA dt           where dN_nu is the number of photens and dA is the collecting area of the telescope
#  dN_nu,e   =  T(nu)  f_nu/(h nu) dnu dA dt    where dN_nu,e is the number of registered electrons and T(nu) is the throughput
#   
#  N_f  =  5.48E10  *  10^(-0.4 Ab_nu)  *  (dA/m^2)  *  (dt/s)  *  INTEGRAL[ T(nu) d ln(lambda) ]
#  N_f  ~  5.48E10  *  10^(-0.4 AB_nu)  *  (dA/m^2)  *  (dt/s)  *  Q
#  N_f  ~  5.48E10  *  10^(-0.4 AB_nu)  *  (dA/m^2)  *  (dt/s)  *  q  * lambda_width / lambda_central
#  where  N_f  is the number of photons through a rectangular filter with lambda_cen and delta_lambda,
#  q is the approximate mean transmission of the total system in the filter, and the source has constant f_nu
#  (for the correct SDSS Q, see  Gunn et al. 1998, AJ 116, 3040)
#   
#  ZP = 2.5*LOG10(Q*dA/[m^2]*5.48E10)


print('################################################################################')

############################ primary inputs 

#tel='WST 2m'
tel='WST 40cm'
#tel='WST 40cm defocussed'
#tel='WST 2m defocussed'

filter   = ['u', 'g', 'r', 'i', 'z']

print('telescope    =',tel)
print('filter       = [  {:>6s}   , {:>6s}   , {:>6s}   , {:>6s}   , {:>6s}    ]'.format(filter[0],filter[1],filter[2],filter[3],filter[4]))


SN_in   = []
SN_in   = 200. * np.ones(5)   # S/N in aperture

mag_in  = []
# mag_in  = np.array([ 24.878, 25.460, 24.998, 24.426, 23.536]) #  object AB [AB mag]
mag_in  = 10. * np.ones(5)  #  object AB [AB mag]

texp_in = []
# texp_in   = np.array([ 360, 360, 360, 360, 360]) # exposure time per image, sec
# texp_in   = 10.* np.ones(5)  # exposure time per image, sec
# texp_in   = [    4.32 ,   1.05 ,   1.02 ,   1.58 ,   3.60  ]
# texp_in   = [5.78685715, 1.08666208, 0.94453801, 1.70828441, 0]

print('--------------------------------------------------------------------------------')
if len(mag_in)==5 :
  print('mag_in       = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(mag_in[0],mag_in[1],mag_in[2],mag_in[3],mag_in[4]))
if len(SN_in)==5 :
  print('SN_in        = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(SN_in[0],SN_in[1],SN_in[2],SN_in[3],SN_in[4]))
if len(texp_in)==5  :
  print('texp_in      = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(texp_in[0],texp_in[1],texp_in[2],texp_in[3],texp_in[4]))

  
#################################### secondary inputs

n        = 1 * np.ones(5)  # number of images []
aperphot = True
if tel=='WST 2m' :
  psf         = 0.90 * np.ones(5) #  seeing/fwhm [arcsec]
  Aper_arcsec = 1.10 * np.ones(5) # aperture diameter [arcsec]
if tel=='WST 2m defocussed' :
  psf         = 10. * np.ones(5) #  seeing/fwhm [arcsec]
  Aper_arcsec = 20. * np.ones(5) # aperture diameter [arcsec]
if tel=='WST 40cm' :
  psf         = 1.50 * np.ones(5) #  seeing/fwhm [arcsec]
  Aper_arcsec = 2.00 * np.ones(5) # aperture diameter [arcsec]
if tel=='WST 40cm defocussed' :
  psf         = 10. * np.ones(5) #  seeing/fwhm [arcsec]
  Aper_arcsec = 20. * np.ones(5) # aperture diameter [arcsec]
AM       = 1.00 * np.ones(5)   # airmass
#skymu_AB = np.array([ 22.80, 21.90, 20.85, 20.15, 19.26])      # nightsky AB [AB mag / arcsec^2]
#skymu_AB = np.array([ 21.00, 21.00, 21.00, 21.00, 21.00])      # nightsky AB [AB mag / arcsec^2]
skymu_AB = np.array([ 20.00, 20.00, 20.00, 20.00, 20.00])      # nightsky AB [AB mag / arcsec^2]
galmu_AB = np.array([ 100.00, 100.00, 100.00, 100.00, 100.00]) # SB galaxy (e.g. M31) background mu_AB  [AB mag / arcsec^2], 100 = no galaxy light
ext      = np.array([ 0.56, 0.18, 0.10, 0.08, 0.07])           #  extinction
# AB_Vega  = np.array([ 0.000, -0.090, 0.145, 0.363, 0.521])     #  mag_AB-mag_Vega
# m_Vega   = mag_in - AB_Vega                            # object Vega [Vega mag]

print('--------------------------------------------------------------------------------')
print('number       = [  {:5.0f}    , {:5.0f}    , {:5.0f}    , {:5.0f}    , {:5.0f}     ]'.format(n[0],n[1],n[2],n[3],n[4]))
print('psf  ["]     = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(psf[0],psf[1],psf[2],psf[3],psf[4]))
print('Aper ["]     = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(Aper_arcsec[0],Aper_arcsec[1],Aper_arcsec[2],Aper_arcsec[3],Aper_arcsec[4]))
print('AM           = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(AM[0],AM[1],AM[2],AM[3],AM[4]))
print('skymu_AB     = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(skymu_AB[0],skymu_AB[1],skymu_AB[2],skymu_AB[3],skymu_AB[4]))
print('galmu_AB     = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(galmu_AB[0],galmu_AB[1],galmu_AB[2],galmu_AB[3],galmu_AB[4]))


#################################### hardware

if tel=='WST 2m' or tel=='WST 2m defocussed' :
  D_m =    2.00 * np.ones(5) # mirror radius, m
  s_m = 0.14056 * np.ones(5) # secondary mirror blockageso that dA=2.7
  pxsc =   0.20 * np.ones(5) # pixelscale [arcsec/px]
  RN_px =  2.60 * np.ones(5) # read-noise px  [N_e / pixel]
  q = np.array([ 0.177, 0.367, 0.441, 0.271, 0.179 ]) # mean throughput in band [q]  
  
if tel=='WST 40cm' or tel=='WST 40cm defocussed' :
  # https://planewave.com/products-page/telescopes/17-inch-cdk-optical-tube-assembly/#.XrL1dprgp3k
  D_m =           0.432 * np.ones(5)  # mirror radius, m
  s_m = 0.190**2/D_m**2 * np.ones(5)  # secondary mirror blockageso that dA=2.7
  pxsc =           0.63 * np.ones(5)  # pixelscale [arcsec/px]
  RN_px =          5.00 * np.ones(5)  # read-noise px  [N_e / pixel]
  q = np.array([ 0.15, 0.35, 0.45, 0.30, 0.0 ]) # mean throughput in band [q]  

#################################### calculations

dA = (D_m/2)**2 *np.pi*(1-s_m) #  effective telescope area mirror area [m^2]

# filters
lam_c = np.array([ 3540.00, 4770.00, 6180.00, 7590.00, 8970.00]) # filter, lambda_central
lam_w = np.array([ 600.00, 1300.00, 1400.00, 1400.00, 1500.00]) # filter, lambda_width
#
Q = q * lam_w / lam_c # throughput 
ZPcalc = 2.5*np.log10(Q*dA*54800000000) # ZPcalc [AB mag] 
ZP = ZPcalc
# ZP = np.array([ 24.25, 25.41, 25.36, 24.87, 23.96]) #  ZPmanual [AB mag]


# psf_at_AM = np.where(AM>1, psf*AM**0.6, psf)                # seeing/fwhm [arcsec], =psf*AM^0.6
f_gauss = (1-np.exp(-1*(Aper_arcsec/2)**2/(psf/2.35)**2/2)) # fraction of Gaussian PSF in aperture [N_e,ap / N_e,tot] 
Aeff_arcsec = psf*2/np.sqrt(np.log(4))                      # eff. aperture diameter for psf phot [arcsec]
Aarea = np.pi*Aper_arcsec**2/np.sqrt(np.log(4))             # eff. aperture area [arcsec^2]

fak_px = pxsc**2
fak_ap = (Aper_arcsec/2)**2 * np.pi
fak_psf = (Aeff_arcsec/2)**2 * np.pi 

sqpx_I_sqarcsec = (1/pxsc)**2   # [px^2 / arcsec^2]

RN_ap =  np.sqrt(RN_px**2 * fak_ap * sqpx_I_sqarcsec)  #  read-noise aperture per image [N_e / aperture]
RN_psf = np.sqrt(RN_px**2 * fak_psf * sqpx_I_sqarcsec) # read-noise psfphot per image [N_e / eff. aperture]

Nsky_sec = 10**(-0.4*(skymu_AB       -ZP)) * AM   # nightsky [N_e / sec]
Ngal_sec = 10**(-0.4*(galmu_AB+AM*ext-ZP))        # galaxy [N_e / sec]

print('--------------------------------------------------------------------------------')
print('f_gauss      = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(f_gauss[0],f_gauss[1],f_gauss[2],f_gauss[3],f_gauss[4]))
print('RON          = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(RN_px[0],RN_px[1],RN_px[2],RN_px[3],RN_px[4]))
print('q            = [  {:7.1f}% , {:7.1f}% , {:7.1f}% , {:7.1f}% , {:7.1f}%  ]'.format(100.*q[0],100.*q[1],100.*q[2],100.*q[3],100.*q[4]))
#print('Q            = [  {:7.1f}% , {:7.1f}% , {:7.1f}% , {:7.1f}% , {:7.1f}%  ]'.format(100.*Q[0],100.*Q[1],100.*Q[2],100.*Q[3],100.*Q[4]))
print('dA           = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(dA[0],dA[1],dA[2],dA[3],dA[4]))
#print('ZPcalc       = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(ZPcalc[0],ZPcalc[1],ZPcalc[2],ZPcalc[3],ZPcalc[4]))
print('ZP           = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(ZP[0],ZP[1],ZP[2],ZP[3],ZP[4]))


####################################


if len(mag_in)==5 :
  Nobj_sec = 10**(-0.4*(mag_in+AM*ext-ZP))               # object [N_e / sec]

  
if len(texp_in)==5 :
  ttot_in  = n*texp_in  # total exposure time [sec], =3*365*6/12*(1-bad)*(t_d_dark*1/2+t_d_gray*1/4+t_d_bright*1/4)
  Nsky_px  = Nsky_sec * ttot_in * fak_px  # nightsky px [N_e / pixel]
  Nsky_ap  = Nsky_sec * ttot_in * fak_ap  # nightsky aperture [N_e / aperture]
  Nsky_psf = Nsky_sec * ttot_in * fak_psf # nightsky PSF  [N_e / eff. aperture]
  Ngal_px  = Ngal_sec * ttot_in * fak_px        # galaxy px [N_e / pixel] 
  Ngal_ap  = Ngal_sec * ttot_in * fak_ap        # galaxy [ N_e / aperture]
  Ngal_psf = Ngal_sec * ttot_in * fak_psf       # galaxy, N_e / eff. aperture, =10^(-0.4*(g_ab-zp))*d'*('Aeff_arcsec'/2)^2*3.14
  

  
if len(texp_in)==5 and len(mag_in)==5 :

  print('================================================================================')
  print('(S/N):')
  print('   INPUT:')
  print('      mag    = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(mag_in[0],mag_in[1],mag_in[2],mag_in[3],mag_in[4]))
  print('      texp   = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(texp_in[0],texp_in[1],texp_in[2],texp_in[3],texp_in[4]))
  print('      number = [  {:5.0f}    , {:5.0f}    , {:5.0f}    , {:5.0f}    , {:5.0f}     ]'.format(n[0],n[1],n[2],n[3],n[4]))
  #print('      ttot   = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(ttot_in[0],ttot_in[1],ttot_in[2],ttot_in[3],ttot_in[4]))

  Nsky_px  = Nsky_sec * texp_in * fak_px   # nightsky px [N_e / pixel]
  Ngal_px  = Ngal_sec * texp_in * fak_px   # galaxy px [N_e / pixel]
  Nobj_peak = Nobj_sec * texp_in * (1-np.exp(-1*(np.sqrt(1/np.pi)*2*pxsc/2)**2/(psf/2.35)**2/2))  # object  [N_e / pixel]
  #Ntot_px = Nobj_peak + Nsky_px + Ngal_px
  
  Nobj_ap  = Nobj_sec * ttot_in * f_gauss  # object [N_e / aperture]   
  Nobj_psf = Nobj_sec * ttot_in            # object [N_e full psf]     stimmt nicht ganz!!
  SN_aper = Nobj_ap   / np.sqrt(Nobj_ap +Ngal_ap +Nsky_ap +n*RN_ap**2)      # S/N in aperture
  SN_psf  = Nobj_psf / np.sqrt(Nobj_psf+Ngal_psf+Nsky_psf+n*RN_psf**2)  # S/N in psfphot

  print('   RESULT:')
  if aperphot :
    print('      SN     = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(SN_aper[0],SN_aper[1],SN_aper[2],SN_aper[3],SN_aper[4]))
  else :
    print('      SNpsf  = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(SN_psf[0],SN_psf[1],SN_psf[2],SN_psf[3],SN_psf[4]))
  print('      Npeak  = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(Nobj_peak[0],Nobj_peak[1],Nobj_peak[2],Nobj_peak[3],Nobj_peak[4]))
  print('      Nsky   = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(Nsky_px[0],Nsky_px[1],Nsky_px[2],Nsky_px[3],Nsky_px[4]))

    
if len(texp_in)==5 and  len(SN_in)==5 :

  print('================================================================================')
  print('LIMITING MAGNITUDE:')
  print('   INPUT:')
  print('      SN     = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(SN_in[0],SN_in[1],SN_in[2],SN_in[3],SN_in[4]))
  print('      texp   = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(texp_in[0],texp_in[1],texp_in[2],texp_in[3],texp_in[4]))
  print('      number = [  {:5.0f}    , {:5.0f}    , {:5.0f}    , {:5.0f}    , {:5.0f}     ]'.format(n[0],n[1],n[2],n[3],n[4]))
  #print('      ttot   = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(ttot_in[0],ttot_in[1],ttot_in[2],ttot_in[3],ttot_in[4]))

  # SN_ap = Nobj_sec*texp* f_gauss * np.sqrt(n) / np.sqrt(Nobj_sec*texp* f_gauss  + Ngal_sec*texp*fak_ap + Nsky_sec*texp*fak_ap + RN_ap**2) # S/N in psfphot
  # SN_ap**2*(Nobj_sec*texp* f_gauss  + Ngal_sec*texp*fak_ap + Nsky_sec*texp*fak_ap + RN_ap**2) = Nobj_sec**2 * f_gauss**2 * texp**2 * n 
  # Nobj_sec**2 * n * f_gauss**2 * texp**2 - SN_ap**2*(Nobj_sec*f_gauss + Ngal_sec*fak_ap + Nsky_sec*fak_ap)*texp - SN_ap**2*RN_ap**2 = 0

  Nobj_ap_SN =(SN_in**2+SN_in*np.sqrt(SN_in**2+4*(Nsky_ap+Ngal_ap+n*RN_ap**2)))/2  # object [N_e / aperture] 
  Nobj_psf_SN = (SN_in**2+SN_in*np.sqrt(SN_in**2+4*(Nsky_psf+Ngal_psf+n*RN_psf**2)))/2  # object [ N_e  full psf ]
  mag_ap = -2.5*np.log10(Nobj_ap_SN/f_gauss/ttot_in)-AM*ext+ZP # limiting mag, AB mag
  mag_psf = -2.5*np.log10(Nobj_psf_SN/ttot_in)-AM*ext+ZP # limiting mag, AB mag

  if aperphot :
    mag = mag_ap
  else :
    mag = mag_psf
  Nobj_sec = 10**(-0.4*(mag+AM*ext-ZP))               # object [N_e / sec]
  Nobj_peak = Nobj_sec * texp_in * (1-np.exp(-1*(np.sqrt(1/np.pi)*2*pxsc/2)**2/(psf/2.35)**2/2))  # object  [N_e / pixel]
  Nsky_px  = Nsky_sec * texp_in * fak_px   # nightsky px [N_e / pixel]
  Ngal_px  = Ngal_sec * texp_in * fak_px   # galaxy px [N_e / pixel] 
  #Ntot_px = Nobj_peak + Nsky_px + Ngal_px

  print('   RESULT:')
  if aperphot :
    print('      mag    = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(mag_ap[0],mag_ap[1],mag_ap[2],mag_ap[3],mag_ap[4]))
  else :
    print('      magpsf = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(mag_psf[0],mag_psf[1],mag_psf[2],mag_psf[3],mag_psf[4]))
  print('      Npeak  = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(Nobj_peak[0],Nobj_peak[1],Nobj_peak[2],Nobj_peak[3],Nobj_peak[4]))
  print('      Nsky   = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(Nsky_px[0],Nsky_px[1],Nsky_px[2],Nsky_px[3],Nsky_px[4]))
  
if len(SN_in)==5 and len(mag_in)==5 :

  print('================================================================================')
  print('EXPOSURE TIME:')
  print('   INPUT:')
  print('      mag    = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(mag_in[0],mag_in[1],mag_in[2],mag_in[3],mag_in[4]))
  print('      SN     = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(SN_in[0],SN_in[1],SN_in[2],SN_in[3],SN_in[4]))
  print('      number = [  {:5.0f}    , {:5.0f}    , {:5.0f}    , {:5.0f}    , {:5.0f}     ]'.format(n[0],n[1],n[2],n[3],n[4]))

  # SN_ap = Nobj_sec*texp* f_gauss * np.sqrt(n) / np.sqrt(Nobj_sec*texp* f_gauss  + Ngal_sec*texp*fak_ap + Nsky_sec*texp*fak_ap + RN_ap**2) # S/N in psfphot
  # SN_ap**2*(Nobj_sec*texp* f_gauss  + Ngal_sec*texp*fak_ap + Nsky_sec*texp*fak_ap + RN_ap**2) = Nobj_sec**2 * f_gauss**2 * texp**2 * n 
  # Nobj_sec**2 * n * f_gauss**2 * texp**2 - SN_ap**2*(Nobj_sec*f_gauss + Ngal_sec*fak_ap + Nsky_sec*fak_ap)*texp - SN_ap**2*RN_ap**2 = 0
  a =  Nobj_sec**2 * n * f_gauss**2 
  b = -SN_in**2*(Nobj_sec*f_gauss + Ngal_sec*fak_ap + Nsky_sec*fak_ap)
  c = -SN_in**2*RN_ap**2
  det = b**2-4*a*c
  texp = (-b+np.sqrt(det))/(2.*a)
  ttot = n*texp
  
  Nsky_px  = Nsky_sec * texp * fak_px   # nightsky px [N_e / pixel]
  Ngal_px  = Ngal_sec * texp * fak_px   # galaxy px [N_e / pixel] 
  Nobj_peak = Nobj_sec * texp * (1-np.exp(-1*(np.sqrt(1/np.pi)*2*pxsc/2)**2/(psf/2.35)**2/2))  # object  [N_e / pixel]
  #Ntot_px = Nobj_peak + Nsky_px + Ngal_px
  
  print('   RESULT:')
  #print('      texp   =',texp)
  print('      texp   = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(texp[0],texp[1],texp[2],texp[3],texp[4]))
  print('      Npeak  = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(Nobj_peak[0],Nobj_peak[1],Nobj_peak[2],Nobj_peak[3],Nobj_peak[4]))
  print('      Nsky   = [  {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}  , {:7.1f}   ]'.format(Nsky_px[0],Nsky_px[1],Nsky_px[2],Nsky_px[3],Nsky_px[4]))
  #print('      ttot   = [  {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f} , {:8.2f}  ]'.format(ttot[0],ttot[1],ttot[2],ttot[3],ttot[4]))

      
# print('   M31:')
# dmod = 2.5*np.log10(770000**2/10**2) # distance modulus galaxy M31
# print('dmod=',dmod)
# lim_M_AB_ap = mag_ap - dmod           # limiting mag absolute in M31, AB mag, 
# lim_M_AB_psf = mag_psf - dmod          # limiting mag absolute in M31, AB mag, 
# lim_M_Vega_ap = lim_M_AB_ap - AB_Vega   # limiting mag_Vega  absolute M31, Vega mag,
# lim_M_Vega_psf = lim_M_AB_psf - AB_Vega  # limiting mag_Vega  absolute M31, Vega mag, 
# print('lim_M_Vega_psf=',lim_M_Vega_psf)
  
print('################################################################################')

 

