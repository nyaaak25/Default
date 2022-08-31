;------------------------------------------------------------------------------------------------------------------------
function ret_pressure, TA, TB, SZA, EA, PA, Dust, Waterice, Albedo
;------------------------------------------------------------------------------------------------------------------------
; fitting function

  T1 = TA
  T2 = TB

  ;result
  pressure_CD = -999d
  Error_U = -999d
  Error_L = -999d

  ;skip
  if T1 gt 285. then begin
    print,'Warning: T1 > limit'
    goto, skip
  endif
  if T1 lt 135. then begin
    print,'Warning: T1 < limit'
    goto, skip
  endif
  if T2 gt 200. then begin
    print,'Warning: T2 > limit'
    goto, skip
  endif
  if T2 lt 80. then begin
    print,'Warning: T2 < limit'
    goto, skip
  endif
  if SZA gt 75. then begin
    print,'Warning: SZA > limit'
    goto, skip
  endif
  if EA gt 50. then begin
    print,'Warning: EA > limit'
    goto, skip
  endif
  if PA gt 180. then begin
    print,'Warning: PA > limit'
    goto, skip
  endif
  if Dust gt 1.5 then begin
    print,'Warning: DUST > limit'
    goto, skip
  endif
  if Waterice gt 1.0 then begin
    print,'Warning: ICE > limit'
    goto, skip
  endif
  if Albedo gt 0.6 then begin
    print,'Warning: Albedo > limit'
    goto, skip
  endif
  if Albedo lt 0.05 then begin
    print,'Warning: Albedo < limit'
    goto, skip
  endif

  x1a = dblarr(5) ; T1
  x1a(0) = (135.d)^1.5d
  x1a(1) = (160.d)^1.5d
  x1a(2) = (213.d)^1.5d
  x1a(3) = (260.d)^1.5d
  x1a(4) = (285.d)^1.5d

  x2a = dblarr(3) ; T2
  x2a(0) = (80.d)^1.5d
  x2a(1) = (146.d)^1.5d
  x2a(2) = (200.d)^1.5d

  x3a = dblarr(6) ; SZA
  x3a(0) = exp(-cos(0d/180d*!dpi))
  x3a(1) = exp(-cos(15d/180d*!dpi))
  x3a(2) = exp(-cos(30d/180d*!dpi))
  x3a(3) = exp(-cos(45d/180d*!dpi))
  x3a(4) = exp(-cos(60d/180d*!dpi))
  x3a(5) = exp(-cos(75d/180d*!dpi))

  x4a = dblarr(5) ; EA
  x4a(0) = exp(-cos(0d/180d*!dpi))
  x4a(1) = exp(-cos(5d/180d*!dpi))
  x4a(2) = exp(-cos(10d/180d*!dpi))
  x4a(3) = exp(-cos(30d/180d*!dpi))
  x4a(4) = exp(-cos(50d/180d*!dpi))

  x5a = dblarr(5) ; PA
  x5a(0) = -cos(0d/180d*!dpi) + 1.d
  x5a(1) = -cos(45d/180d*!dpi) + 1.d
  x5a(2) = -cos(90d/180d*!dpi) + 1.d
  x5a(3) = -cos(135d/180d*!dpi) + 1.d
  x5a(4) = -cos(180d/180d*!dpi) + 1.d

  x6a = dblarr(6) ; Dust
  x6a(0) = 0.0d
  x6a(1) = 0.3d
  x6a(2) = 0.6d
  x6a(3) = 0.9d
  x6a(4) = 1.2d
  x6a(5) = 1.5d

  x7a = dblarr(3) ; Water-ice
  x7a(0) = 0.0d
  x7a(1) = 0.5d
  x7a(2) = 1.0d

  x8a = dblarr(7) ; Surface Albedo
  x8a(0) = 0.05d
  x8a(1) = 0.1d
  x8a(2) = 0.2d
  x8a(3) = 0.3d
  x8a(4) = 0.4d
  x8a(5) = 0.5d
  x8a(6) = 0.6d

  ;Surface pressure grid
  Pressure_grid = dblarr(15)
  Pressure_grid(0) = alog(50d)
  Pressure_grid(1) = alog(150d)
  Pressure_grid(2) = alog(180d)
  Pressure_grid(3) = alog(215d)
  Pressure_grid(4) = alog(257d)
  Pressure_grid(5) = alog(308d)
  Pressure_grid(6) = alog(369d)
  Pressure_grid(7) = alog(442d)
  Pressure_grid(8) = alog(529d)
  Pressure_grid(9) = alog(633d)
  Pressure_grid(10) = alog(758d)
  Pressure_grid(11) = alog(907d)
  Pressure_grid(12) = alog(1096d)
  Pressure_grid(13) = alog(1300d)
  Pressure_grid(14) = alog(1500d)

  y = dblarr(15)
  X1 = (T1)^1.5d ;T1 as test
  X2 = (T2)^1.5d ;T2 as test
  x3 = exp(-cos(SZA/180d*!dpi))
  x4 = exp(-cos(EA/180d*!dpi))
  x5 = -cos(PA/180d*!dpi) + 1.d
  x6 = Dust ;reduce_dust(loop)*1.86
  x7 = Waterice; reduce_waterice(loop)*3.87
  x8 = Albedo

  for I = 0, 5-1 do begin
    F = x1a(I) - X1
    if (F gt 0.0d0) then j1 = I
    if (F gt 0.0d0) then goto, skip1
  endfor
  skip1:

  for I = 0, 3-1 do begin
    F = x2a(I) - X2
    if (F gt 0.0d0) then j2 = I
    if (F gt 0.0d0) then goto, skip2
  endfor
  skip2:

  for I = 0, 6-1 do begin
    F = x3a(I) - X3
    if (F gt 0.0d0) then j3 = I
    if (F gt 0.0d0) then goto, skip3
  endfor
  skip3:

  for I = 0, 5-1 do begin  ; !!TBD
    F = x4a(I) - X4
    if (F gt 0.0d0) then j4 = I
    if (F gt 0.0d0) then goto, skip4
  endfor
  skip4:

  for I = 0, 5-1 do begin
    F = x5a(I) - X5
    if (F gt 0.0d0) then j5 = I
    if (F gt 0.0d0) then goto, skip5
  endfor
  skip5:

  for I = 0, 6-1 do begin
    F = x6a(I) - X6
    if (F gt 0.0d0) then j6 = I
    if (F gt 0.0d0) then goto, skip6
  endfor
  skip6:

  for I = 0, 3-1 do begin
    F = x7a(I) - X7
    if (F gt 0.0d0) then j7 = I
    if (F gt 0.0d0) then goto, skip7
  endfor
  skip7:

  for I = 0, 7-1 do begin
    F = x8a(I) - X8
    if (F gt 0.0d0) then j8 = I
    if (F gt 0.0d0) then goto, skip8
  endfor
  skip8:

;  j1 = j1 -1
;  j2 = j2 -1
;  j3 = j3 -1
;  j4 = j4 -1
;  j5 = j5 -1
;  j6 = j6 -1
;  j7 = j7 -1
;  j8 = j8 -1
;
;  if j1 lt 0 then stop
;  if j2 lt 0 then stop
;  if j3 lt 0 then stop
;  if j4 lt 0 then stop
;  if j5 lt 0 then stop
;  if j6 lt 0 then stop
;  if j7 lt 0 then stop
;  if j8 lt 0 then stop
;  if j1+1 ge 5 then stop
;  if j2+1 ge 3 then stop
;  if j3+1 ge 6 then stop
;  if j4+1 ge 5 then stop
;  if j5+1 ge 5 then stop
;  if j6+1 ge 6 then stop
;  if j7+1 ge 3 then stop
;  if j8+1 ge 6 then stop

  print, j1,j2,j3,j4,j5,j6,j7,j8

  for I = 0, 15-1 do begin
    if i eq 0 then path = '/work1/LUT/SP/table/output_SP1/'
    if i eq 1 then path = '/work1/LUT/SP/table/output_SP2/'
    if i eq 2 then path = '/work1/LUT/SP/table/output_SP3/'
    if i eq 3 then path = '/work1/LUT/SP/table/output_SP4/'
    if i eq 4 then path = '/work1/LUT/SP/table/output_SP5/'
    if i eq 5 then path = '/work1/LUT/SP/table/output_SP6/'
    if i eq 6 then path = '/work1/LUT/SP/table/output_SP7/'
    if i eq 7 then path = '/work1/LUT/SP/table/output_SP8/'
    if i eq 8 then path = '/work1/LUT/SP/table/output_SP9/'
    if i eq 9 then path = '/work1/LUT/SP/table/output_SP10/'
    if i eq 10 then path = '/work1/LUT/SP/table/output_SP11/'
    if i eq 11 then path = '/work1/LUT/SP/table/output_SP12/'
    if i eq 12 then path = '/work1/LUT/SP/table/output_SP13/'
    if i eq 13 then path = '/work1/LUT/SP/table/output_SP14/'
    if i eq 14 then path = '/work1/LUT/SP/table/output_SP15/'


    wn=dblarr(27)
    rad=dblarr(27)
    y = dblarr(27,15)

    file = path + 'SP' +  STRCOMPRESS(fix(I+1),/REMOVE_AL) $
                          + '_TA' + STRCOMPRESS(fix(j1),/REMOVE_AL) $
                          + '_TB' + STRCOMPRESS(fix(j2),/REMOVE_AL) $
                          + '_SZA' + STRCOMPRESS(fix(j3),/REMOVE_AL) $
                          + '_EA' + STRCOMPRESS(fix(j4),/REMOVE_AL) $
                          + '_PA' + STRCOMPRESS(fix(j5),/REMOVE_AL) $
                          + '_Dust' + STRCOMPRESS(fix(j6),/REMOVE_AL) $
                          + '_WaterI' + STRCOMPRESS(fix(j7),/REMOVE_AL) $
                          + '_SurfaceA' + STRCOMPRESS(fix(j8),/REMOVE_AL) + '_rad.dat'
                          
    print, file
    ;read files
    openr,lun,file,/get_lun
    for i = 0, 27-1 do begin
      readf,lun,a,b
      wn(i)=a
      rad(i)=b
    endfor
    free_lun,lun

    wn = (1/wn)*10000
    wn = reverse(wn)
    rad = reverse(rad)

    y(*,I-1) = rad

  endfor


stop

skip:
return, rad

end



Pro retrieval_pressure_fitting
; 必要なものだけを取り出す
;observation OMEGA data path
path = '/data2/omega/sav/'
path_output = '/work1/LUT/SP/table/absorption/'

;Loop start
for Loop = 0, 3 do begin
  if loop eq 0 then file = path+'ORB0931_3.sav'
  if loop eq 1 then file = path+'ORB0030_1.sav'
  if loop eq 2 then file = path+'ORB0920_3.sav'
  if loop eq 3 then file = path+'ORB0313_4.sav'
  restore, file

  ip = n_elements(LATI(*,0))
  io = n_elements(LATI(0,*))

  CO2=where(wvl gt 1.81 and wvl lt 2.19)
  wvl=wvl[CO2]
  jdat=jdat(*,CO2,*)

  restore, path + 'specmars.sav'
  specmars = specmars(CO2)
  
  nanserch=where(jdat ge 0 and jdat le 0.0001)
  jdat(nanserch)=!VALUES.F_NAN
  
  ;latitude grid
  span_lati = max(lati(0,*)) - min(lati(0,*))
  nlay_span = ceil(span_lati/1.875)

  for n = 0, nlay_span-1 do begin
    points = where(lati(0,*) ge min(lati(0,*)) + float(n)*1.875 and lati(0,*) lt min(lati(0,*)) + float(n+1)*1.875, count)
    count_nscan = count
 
;    ; ----- MCD ------
      lat = median(lati(*,points),/double)
      lon = median(longi(*,points),/double)
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
      ;    ;-----------------

    for i = 0, count_nscan-1 do begin ;loop for slit scan
      for j = 0, ip-1 do begin  

        SZA = reform(geocube(j,8,points(i)))*1.e-4
        EA = reform(geocube(j,9,points(i)))*1.e-4
        PA = reform(geocube(j,10,points(i)))*1.e-4
        
        ; reflectance factor ( I/F / cos(SZA) )
        Albedo_input = jdat(j,0,points(i))/specmars(0) / cos(geocube(j,8,points(i))*1e-4*!DTOR)
        pressure(j,points(i)) = ret_pressure(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input)

      endfor
    endfor
  endfor
endfor

; result = MPFITEXPR(expr,t,r)
end