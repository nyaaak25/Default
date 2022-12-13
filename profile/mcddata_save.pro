;------------------------------------------------------------------------------------------------------------------------
;
; mcd data をsaveする
; dust opacity map
;
; create by Akira Kazama
;
; mcddata_save　 ::2022.12.12 Mon 18:54:00
; 
; +++
;------------------------------------------------------------------------------------------------------------------------


Pro mcddata_save

start_time = systime(1)
; define color platte
Set_Plot, 'x'
device, retain=1, decomposed=0
loadct, 39

; path 
path = '/data2/omega/sav/'
path2 = '/work1/LUT/SP/table/absorption/'

; 相対誤差の比較のためのORB [Forget+ 2007] Fig 13
;restore, path+'ORB0920_3.sav'
restore, path+'ORB0931_3.sav'

restore, path + 'specmars.sav'

;ORB0920_3
;ind = where_xyz(longi ge 273 and longi le 277 and lati ge 53 and lati le 56, xind=xind, yind=yind)
ind = where_xyz(longi ge 274 and longi le 277 and lati ge 54 and lati le 56, xind=xind, yind=yind)

; ORB0931_3
;ind = where_xyz(longi ge 272 and longi le 277 and lati ge 50 and lati le 61, xind=xind, yind=yind)

; 全部のやつ
ip = n_elements(LATI(*,0))
io = n_elements(LATI(0,*))

; loop変数
ip_b = min(xind)
ip_a = max(xind)
io_b = min(yind)
io_a = max(yind)

; create array
MCDpressure = dblarr(ip,io)
dustmap = dblarr(ip,io)
TAmap = dblarr(ip,io)
altitude = dblarr(ip,io)
albedo = dblarr(ip,io)

for l = ip_b, ip_a -1 do begin ;loop for slit scan
  for k = io_b, io_a -1 do begin ;test getting surface feature

    start_time = systime(1)
    ;   ----- MCD ------
    lat = lati(l,k)
    lon = longi(l,k)
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
    ;    cd, ‘/Users/Shohei/Documents/Programs/MCD_ver4.3/mcd/idl/’
    ;   mcd_idl, Ls, Loct, lat, lon, dust, zkey, hrkey, xz, meanvarz, extvarz
    cd, '/work1/LUT/MCD/MCD5.3/mcd/idl/'
    a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
      dset,scena,perturkey,seedin,gwlength,extvarkeys, $
      pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
    scaleH = extvar(13)
    Z1 = scaleH * 0.1d
    Z2 = scaleH * 4d
    dust_opacity = extvar(36)*1.2090217E+00/4.6232791E+00
    ice_opacity = 0.d;*4.2410289E-01/***  ;!TBD!
    SP = PRES
    xz=[Z1] ; (array of altitude with step=2000 m)
    a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
      dset,scena,perturkey,seedin,gwlength,extvarkeys, $
      pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
    TA = temp
    xz=[Z2] ; (array of altitude with step=2000 m)
    a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
      dset,scena,perturkey,seedin,gwlength,extvarkeys, $
      pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
    TB = temp
    ;;-----------------

    CO2=where(wvl ge 1.8 and wvl le 2.2)
    wvl_CO2=wvl[CO2]
    jdat_CO2=jdat(*,CO2,*)
    specmars_CO2 = specmars(CO2)

    altitude(l,k) = reform(geocube(l,12,k))
    Albedo_input = jdat_CO2(l,0,k)/specmars_CO2(0) / cos(geocube(l,8,k)*1e-4*!DTOR)
    dust = dust_opacity

    albedo(l,k) = Albedo_input
    MCDpressure(l,k) = SP
    dustmap(l,k) = dust 
    TAmap(l,k) = TA

    end_time = systime(1)
    save, lati, longi, altitude, dustmap, TAmap, MCDpressure, albedo, filename='/work1/LUT/SP/table/absorption/mcddata_ORB0931_3.sav'

  endfor
endfor


stop
end
