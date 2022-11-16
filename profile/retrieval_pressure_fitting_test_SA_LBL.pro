;---------------------------
function forward_LBL, x, p
;start = [SP, TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input, Result_Fit0(0), Result_Fit0(1), Result_Fit0(2)]
;---------------------------
openw,lun,'x_tmp.dat',/get_lun
for j = n_elements(x)-1, 0, -1 do printf, lun, 1d4/x(j)
free_lun,lun

restore, '/work1/LUT/Common/HTP_tmp.sav'

T1 = P(1)
T2 = P(2)
SZA = P(3)
EA = P(4)
PA = P(5)
Dust = P(6)
Waterice = P(7)
Albedo = P(8)

;run ARS
;factor_CO2 = 1.0
openw,lun,'arsm_Mars.ini',/get_lun
printf,lun," 0 = .Info_Level."                                                                                                             
printf,lun,"' ' = .File_Info."
printf,lun,"' ' = .File_Status. ! debug and multiple processes run control purpose"
printf,lun,"0 = .Intel_VML. ! 0 or 1"
printf,lun,"                       "
printf,lun," ===   Output Files   ===  "
printf,lun,"                       "
printf,lun," ''  = .File_Transmittance. "
printf,lun," ''  = .File_Monochromatic_Transmittance. "
printf,lun," 'arsm_rad.dat'  = .File_Monochromatic_Radiance. "
printf,lun," 'arsm_rad_conv.dat'  = .File_Radiance.           "
printf,lun,"                                               "
printf,lun,"  ===   Wavenumber range and radiative transfer option   ================    "
printf,lun,"                             "
printf,lun,"   4500  = .Vmin.             "
printf,lun,"   5600  = .Vmax.            "
printf,lun,"   1     = .Geometry. ! 1 only"
printf,lun,"   1     = .Planck_Source.   "
printf,lun,"   1     = .Solar_Source.    "
printf,lun,"   4000 =  .V_Planck_Max.     "
printf,lun,"   1000 =  .V_Solar_Min.      "
printf,lun,"   2 4 2 = .Radiative_Transfer_Options. ! see the manual, e.g.:     "
printf,lun,"                                       ! 1 1 0, 1 2 0, 2 1 1, 2 2 1, 2 3 1, 2 4 1, 3 0 0  "
printf,lun,"                                       !2 4 1 ottima per marte.   "
printf,lun,"   16384 = .Computational_Subinterval_Size.      "
printf,lun,"   1.1   = .Gas_Subinterval_Factor.   "
printf,lun,"                                         "
printf,lun,"  ===   Corr-I Correction   ============================================== "
printf,lun,"                                                         "
printf,lun,"   0     = .Cor_I_Partition.                  ! 10 is ok "
printf,lun,"   1000   = .Cor_I_Max_Subinterval_Size.       ! cm-1    "
printf,lun,"   0      = .Interpolate_Between_Subintervals.           "
printf,lun,"                                      "
printf,lun,"  =====   Model condition parameters   ================================== "
printf,lun,"   300.00   = .Observer_Altitude.             "
printf,lun,"   " + strcompress(SZA,/REMOVE_ALL)  + " = .Theta_Sun. "
printf,lun,"   " + strcompress(EA,/REMOVE_ALL)  + " = .Theta. "
printf,lun,"   0.00  = .Phi_Sun.             "
printf,lun,"   " + strcompress(PA,/REMOVE_ALL)  + " = .Phi. "
printf,lun,"   1.5 = .Heliocentric_Distance.   " ;
printf,lun,"   0.0 = .Surface_Temperature."                                                                                                 
printf,lun,"                             "
printf,lun,"  '/work1/LUT/Common/pfsolspec_hr.dat' = .File_Solar_Radiance.    "
printf,lun,"  'htp_Mars.atmos' = .File_HTP."
printf,lun,"  'HTP' = .HTP_Level_Selection."
printf,lun,"                                 "
printf,lun,"   " +STRCOMPRESS(Albedo,/REMOVE_ALL)+" = .Surface_Albedo.   "
printf,lun,"  ' ' = .File_Surface_Albedo.   "
printf,lun,"                     "
printf,lun,"  =====   Monochromatic wavenumber grid   ================================  "
printf,lun,"                                                                   "
printf,lun,"  'wn.v' = .File_Wavenumber_Grid.  "
printf,lun,"                                                          "
printf,lun,"  =====   Gases   ========================================================  "
printf,lun,"                                              "
printf,lun,"                                          "
printf,lun,"   1 = .Number_of_Species.                 "
printf,lun,"                                          "
printf,lun,"   2 0 1                                  ! molecule and isotope codes, add.broad."
printf,lun,"  ' '          ! spectral lines file name"
printf,lun,"  'wn.v'                ! specific wavenumber grid file name"
printf,lun,"  'Mars.k'          ! absorption coefficients file name"
printf,lun,"     -1 'Mars.atmos' "+STRCOMPRESS(P(0),/REMOVE_ALL)
printf,lun,"   1 0 0 0 0                                        ! line shape         "
printf,lun,"                                          "
printf,lun,"                                          "
printf,lun,"                           "
printf,lun,"                           "
printf,lun,"  =====   Aerosols   =====================================================  "
printf,lun,"                                 "
printf,lun,"   0 = .Number_of_Aerosols.     "
printf,lun,"                                "
printf,lun,"   3 '/work1/LUT/Common/dust_2500_8500.aero'"
printf,lun,"   1 'Profile_dust.dat' 1"
printf,lun,"                                "
printf,lun,"   3 '/work1/LUT/Common/waterice_2500_8500.aero'"
printf,lun,"   1 'Profile_waterice.dat' 1"
printf,lun,"  =====   Other settings   ===============================================  "
printf,lun,"                                                     "
printf,lun,"   1     = .Observation_Geometry.                    "
printf,lun,"   1     = .Intel_VML.                               "
printf,lun,"  ' '   = .File_Info.                                "
printf,lun,"  ======   V-computation   ===============================================   "
printf,lun,"                     "
printf,lun,"   0.5d+0 = .Alpha.       "
printf,lun,"   1.5d+0 = .Beta.      "
printf,lun,"   100.d+0 = .Delta.       "
printf,lun,"   0.01d+0 = .Step_in_Wing.   "
printf,lun,"   1.e-5 = .Special_Grid_Taumin. "
printf,lun,"                              "
printf,lun,"  ======   K-computation   ===============================================     "
printf,lun,"                             "
printf,lun,"   1 = .Line_Contour.        "
printf,lun,"   0 = .Normalization.       "
printf,lun,"   50 = .Line_Cutoff.       "
printf,lun,"   0 = .Taumin.             "
printf,lun,"   1 = .Inform_Interval.     "
printf,lun,"   0 = .Grid_Centered.       "
printf,lun,"   0 = .HTP_Grid_Reduction.   "
printf,lun,"                      "
printf,lun,"  ===   Limb sounding parameters (Unused)  ===============================    "
printf,lun,"                                                                  "
printf,lun,"   0     = .Tangent_Sounding_Altitude.  [km]    ! for limb sounding only "
printf,lun,"   43.4   = .Mean_Molecular_Weight.      [g//mole]! ---//---    "
printf,lun,"   372.1  = .Gravitational_Acceleration. [cm//s2] ! ---//---     "
printf,lun,"   3387.1 = .Radius_of_Planet.           [km]    ! ---//---     "
printf,lun,"   0      = .Exponential_Tail.                   ! ---//---     "
printf,lun,"                                    "
printf,lun,"1 = .FBJ1J2_Warnings."
printf,lun,"0 = .Continuum_Absorption."
printf,lun,"' ' = .File_Continuum_Absorption."
printf,lun,"  ======================================================================== "
printf,lun,"                                            "
printf,lun,"   3000 = .Sounding_Altitude.      "
printf,lun,"   0    = .Rayleigh_Scattering.     "
printf,lun,"   1    = .Delta_M.                "
printf,lun,"   64    = .Number_of_Streams.      "
printf,lun,"   32   =  .Number_of_Moments.      "
printf,lun,"   0    = .Set_Zero_Scattering.    "
printf,lun,"                                  "
printf,lun,"  ===   Instrument specification   =======================================   "
printf,lun,"                                                      "
printf,lun,"  ! 1-Rect, 2-Trian, 3-Gauss, 4-Sinc, 5-Sinc2, 6-Hamming, 7-Cos2, 8-User,File  "
printf,lun,"                             "
printf,lun,"   3 = .Instrumental_Function.     "
printf,lun,"  '' = .File_Instrumental_Function. "
printf,lun,"   32.5 = .Resolution.    [Ideal PFS: Sinc:1.211376319487D+0, Hamming:1.822245687D+0]   "
;printf,lun,"   30. = .Resolution.    [Ideal PFS: Sinc:1.211376319487D+0, Hamming:1.822245687D+0]   "
printf,lun,"                      [IRIS: Hamming => 2.1254448]    "
printf,lun,"                                         "
printf,lun,"  -100. = .Spectral_Step.                   "
printf,lun,"'x_tmp.dat' = .File_Spectral_Channels.  ! if Spectral_Step < 0"                            
printf,lun,"                                        "
printf,lun,"  ===   Output format   ==================================================   "
printf,lun,"                                                        "
printf,lun,"  'e' = .Units_Output.                          ! for thermal region only  "
printf,lun,"   1 = .Two_Pass_Transmittance_Flag.                "
printf,lun,"  'f12.7' 'f9.4' 'f10.6' '1p,e15.7,0p' 'f8.3' = .Output_Formats. [V,C,T,I,K]  "
printf,lun,"   2   = .Transmittance_Path.                 "
printf,lun,"                                                                       "
printf,lun,"  === Unused parameters (from other programs) ============================ "
printf,lun,"                             "
printf,lun,"   1048576 = .Buffer_Size.     "
printf,lun,"   1 = .Curtis_Godson_Modification.  "
printf,lun,"                                                                           "
printf,lun,"  ===  Other output  =====================================================  "
printf,lun,"                                            "
printf,lun,"                                            "
printf,lun,"                                            "
printf,lun,"  ' '    = .File_Monochromatic_Reflectance. "
printf,lun,"  ' '    = .File_Reflectance."
printf,lun,"                            "
printf,lun,"   0 = .Monochromatic_Dumps."
printf,lun,"' '    = .File_Monochromatic_Optical_Depth."
printf,lun,"' '    = .File_Monochromatic_Flux_Up."
printf,lun,"' '    = .File_Monochromatic_Flux_Down."
printf,lun,"' '    = .File_LibRadtran."
printf,lun,"' '    = .File_Medium_Properties."
printf,lun,""
free_lun,lun

spawn, '/work1/bin/arsm arsm_Mars.ini'

nf = n_elements(x)
f = dblarr(nf)
a1 = 0.d
a2 = 0.d

openr,lun,'arsm_rad_conv.dat',/get_lun
for k = 0l, nf-1 do begin
  readf,lun,a1,a2
  f(k) = a2
endfor
free_lun,lun

f = reverse(f)

restore, '/work1/LUT/Common/specmars_CO2.sav'
F = (F/(x*x*1d-4*1d-4))*1d-7
F = F / specmars 

;F = F * pyroxenes(x, P(9:11))
F = F * poly(findgen(nf), P(9:10))

return, F
end

;---------------------------
function pyroxenes, x, p
;---------------------------
y1 = 1d - gauss1(x, [1.9d, 0.5d / (2 * sqrt(2*alog(2))), P(0)], /peak)
y2 = 1d - gauss1(x, [2.3d, 0.56d / (2 * sqrt(2*alog(2))), P(1)], /peak)
y = (y1 + y2)/2d * P(2)
return, y
end

Pro retrieval_pressure_fitting_test_SA_LBL

Set_Plot, 'x'
device, retain=1, decomposed=0
loadct, 39

path = '/data2/omega/sav/'
path2 = '/work1/LUT/SP/table/absorption/'
restore, path+'ORB0363_3.sav'
restore, path + 'specmars.sav'

;reference spectrum
ref = dblarr(2,3539)
openr, lun, '/work1/LUT/Common/psg_trn.txt', /get_lun
for i = 0, 3539-1 do begin
  readf, lun, a, b
  ref(0,i) = a
  ref(1,i) = b
endfor
free_lun,lun

ind=where_xyz(longi ge 360-48.238 and longi le 360-48.236 and lati ge 22.705 and lati le 22.706, xind=xind, yind=yind)

;test getting surface feature 
x0 = reform(wvl(0:127))
y0 = reform(jdat(xind,0:127,yind)/specmars(0:127))
nanserch=where(y0 ge 0 and y0 le 0.0001)
y0(nanserch)=!VALUES.F_NAN

;interpolate reference spectrum
ref_Mars = interpol(ref(1,*), ref(0,*), x0, /nan)
good = where(x0 ge 1.2 and x0 le 2.6 and ref_Mars ge 0.995)

;window,1
;plot, x0(good), y0(good), thick=3, back=255, color=0, xs=1, ys=1, psym=1, symsize=2, xr=[1.2, 2.6]

pi = replicate({step:0d, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 3)
start = [0.06, 0.06, median(y0(good))]
Result_Fit0 = MPFITFUN('pyroxenes', x0(good), y0(good), y0(good)*1d-2, start, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit, /nan)
F0 = yfit

;oplot, x0(good), F0, color=60, thick=3
;
;window,0
;plot, x0, y0, xr=[1.2, 2.6], thick=3, back=255, color=0, xs=1, ys=1, psym=-1
;;oplot, x0, ref_Mars*0.32, color=254, thick=3
;oplot, x0(good), y0(good), color=100, psym=1, thick=3
;oplot, x0, pyroxenes(x0, Result_Fit0), color=200, thick=3


; CO2 absorption emission line
;CO2=where(wvl gt 1.81 and wvl lt 2.19)
CO2=where(wvl gt 1.81 and wvl lt 2.20)
wvl=wvl[CO2]
jdat=jdat(*,CO2,*)
specmars = specmars(CO2)
save, specmars, filename='/work1/LUT/Common/specmars_CO2.sav'

nanserch=where(jdat ge 0 and jdat le 0.0001)
jdat(nanserch)=!VALUES.F_NAN

; ----- MCD ------
lat = median(lati(xind,yind),/double)
lon = median(longi(xind,yind),/double)
Ls = SOLAR_LONGITUDE
Loct = 12. + (lon - SUB_SOLAR_LONGITUDE)/15. ;! TBC!!! !
dust = 2  ; our best guess MY24 scenario, with solar average conditions
hrkey = 1 ; set high resolution mode on (hrkey=0 to set high resolution off)
zkey = 3    ; specify that xz is the altitude above surface (m)
xz = 0. ; (array of altitude: step=2000 m start: 5m)
datekey = 1       ; <integer> type of input date (1=Mars date)
xdate = ls        ; <double precision> date (IF datekey = 1 : Value of Ls [deg.])
dset ='/work1/LUT/MCD/MCD5.3/data/' ;‘MCD_DATA/’  ; <character*50> data set
scena = 1         ; <integer> scenario (1 = Climatology ave solar)
perturkey = 1     ; <integer>  perturbation type (1= none)
seedin   = 7.0    ; <real>
gwlength = 0.0    ; <real>  for small scale (ie: gravity wave) perturbations;
extvarkeys = LONARR(101) ; add 1 element because indexes in IDL start at 0
for i0 = 0, 100 do extvarkeys[i0] = 1 ; <integer> array output type (extvar(i) = 0 : don’t compute)
meanvar = FLTARR(6) ; add 1 element because indexes in IDL start at 0
for i0 = 0, 5 do meanvar[i0] = 0.0 ; <real> mean unperturbed values (array of 5)
extvar = FLTARR(101) ; add 1 element because indexes in IDL start at 0
for i0 = 0, 100 do extvar[i0] = 0.0 ; <real>  extra variables (array of 100)
pres = 0.0        ; <real> atmospheric pressure (Pa)
dens = 0.0        ; <real> atmospheric density (kg/m^3)
temp = 0.0        ; <real> atmospheric temperature (K)
zonwind = 0.0     ; <real> zonal wind component (East-West)
merwind = 0.0     ; <real> meridional wind component (North-South)
seedout = 0.0     ; <real> current value of the seed of the random number generator
ierr = 0          ; <integer> error flag (0 = OK)

xz = findgen(101)
T_MCD = fltarr(101)
H_MCD = fltarr(101)
P_MCD = fltarr(101)
VMR_MCD = fltarr(101)
Dust_MCD = fltarr(101)
WaterI_MCD = fltarr(101)

cd, '/work1/LUT/MCD/MCD5.3/mcd/idl/'
for i = 0, 100 do begin
  xz = float(i)*1e3
  a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
    dset,scena,perturkey,seedin,gwlength,extvarkeys, $
    pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)

  H_MCD(i) = xz
  T_MCD(i) = temp
  P_MCD(i) = pres
  VMR_MCD(i) = extvar(57)

  V_dust = 4d/3d*3.14d*((1.5d-6)^3.0d)
  Dust_MCD(i) = (dens * extvar(38) / 2500.0d+00 / V_dust * 1.0d-6)

  V_ice = 4d/3d*3.14d*((1.36d-6)^3d)
  WaterI_MCD(i) = pres / (temp * 8.314d) * (extvar(44)) * 18.0d-3 / 917.0d+00 / V_ice  * 1.0d-6
endfor

SZA = reform(geocube(xind,8,yind))*1.e-4
EA = reform(geocube(xind,9,yind))*1.e-4
PA = reform(geocube(xind,10,yind))*1.e-4
Albedo_input = 0.29; dat(xind,0,yind)/specmars(0) / cos(geocube(xind,8,yind)*1e-4*!DTOR)
save, H_MCD, T_MCD, P_MCD, VMR_MCD, Dust_MCD, WaterI_MCD, filename='/work1/LUT/Common/HTP_tmp.sav'

openw,lun,'/work1/LUT/test/Profile_dust.dat',/get_lun
for i0 = 0, 100 do begin
  printf,lun,float(i0),float(Dust_MCD(i0))
endfor
free_lun,lun

openw,lun,'/work1/LUT/test/Profile_waterice.dat',/get_lun
for i0 = 0, 100 do begin
  printf,lun,float(i0),float(WaterI_MCD(i0))
endfor
free_lun,lun

;cal ABS
cd, '/work1/LUT/test/'
openw,lun,'atm_tmp.dat',/get_lun
for j = 0, 100 do printf, lun, T_MCD(j), P_MCD(j), VMR_MCD(j)
free_lun,lun
;spawn, 'rm -r -f wn.v'
;spawn, 'rm -r -f Mars.k'
;spawn, 'rm -r -f Mars.atmos'
;spawn, 'rm -r -f htp_Mars.atmos'
;spawn, './ABS_Mars.out'

x = wvl
y = reform(jdat(xind(0), *, yind(0)))
y = y/specmars

y(0:2) = -0d/0d
y(7) = -0d/0d
y(23) = -0d/0d

err = y*1d-3
err(*) = median(y)*1d-3

;retrieval
;pi = replicate({step:0d, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 12)
;start = [1d, 0., 0., SZA, EA, PA, 0, 0, Albedo_input, Result_Fit0(0), Result_Fit0(1), Result_Fit0(2)]
;;pi(1:10).fixed = 1
;;pi(0:8).fixed = 1
;pi(1:8).fixed = 1

pi = replicate({step:0d, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 11)
start = [1d, 0., 0., SZA, EA, PA, 0, 0, Albedo_input, 1, 0]
pi(1:8).fixed = 1

pi(0).step = 0.01d


Result_Fit = MPFITFUN('forward_LBL', x, y, err, start, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit, /nan)
F = yfit

good = where(FINITE(y) eq 1)
window,6, xs=800,ys=800
plot, x(good), y(good), yr=[-0.1,0.3], back=255, color=0, thick=3, psym=-1, xr=[1.85, 2.2], xs=1
oplot, x(good), f(good), color=0, thick=3, psym=-1, linestyle=2
oplot, x(good), (y(good)-f(good))*5, color=200, thick=2, psym=-1
;oplot, x(good), pyroxenes(x(good), Result_Fit(9:11)), color=100
xyouts, 2.05, 0.1, 'SP='+strcompress(P_MCD(0)*Result_Fit(0)), charsize=2, color=0
stop

window,6, xs=800,ys=800
plot, x(good), pyroxenes(x(good), Result_Fit(9:11)), yr=[0,1], back=255, color=0, thick=3, psym=-1, xr=[1.85, 2.2], xs=1
;print, P_MCD(0)*Result_Fit(0)
stop

end
