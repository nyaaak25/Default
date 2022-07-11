Pro retrieval_water_test_vevr2

;++++++++++++++++++++++++
;Written by SHOHEI AOKI
;++++++++++++++++++++++++

;================================================
;path
;================================================
;path = '/Users/Shohei/tmp/test_case/'
path = '/Users/shohei/tmp/'
restore, path + 'specmars.sav'

;restore
restore,'/Users/shohei/tmp/Table_water_EW_v1.sav'

;para
kb = 1.38d-23
R = 192.0d+00
g = 3.72d+00
P_layer = dblarr(51)
H_layer = dindgen(51)*1d3
T_layer = dblarr(51)
nn = 0l
Water_CD_tot = dblarr(248584)
Lat_tot = dblarr(248584)

;================================================
;Loop start
;================================================
for Loop = 0, 4 do begin
  
  if loop eq 0 then file = path+'ORB1023_3.sav'
  if loop eq 1 then file = path+'ORB1023_4.sav'
  if loop eq 2 then file = path+'ORB1023_5.sav'
  if loop eq 3 then file = path+'ORB1023_6.sav'
  if loop eq 4 then file = path+'ORB1023_7.sav'
  if loop eq 5 then file = path+'ORB1983_1.sav'
;  if loop eq 5 then file = path+'ORB0018_4.sav'
;  if loop eq 6 then file = path+'ORB0024_4.sav'
;  if loop eq 7 then file = path+'ORB0049_4.sav'
;  if loop eq 8 then file = path+'ORB0171_3.sav'

  restore, file
  ip = n_elements(LATI(*,0))
  io = n_elements(LATI(0,*))
  
  ;wvl(113) 2.51 um
  ;wvl(115) 2.54 um
  ;wvl(120) 2.60 um
  ;wvl(123) 2.64 um
  ;*wvl(121) 2.48768 um - bad pixle

  trans = reform(jdat(*,0,*))
  water_VMR = trans
  water_CD = trans
  water_VMR(*,*) = -0d/0d
  water_CD(*,*) = -0d/0d
  X = [wvl(115), wvl(123)]
  
  ;latitude grid
  span_lati = max(lati(0,*)) - min(lati(0,*))
  nlay_span = ceil(span_lati/1.875)

  for n = 0, nlay_span-1 do begin
  ;for n = nlay_span-1, 0, -1 do begin
  ;for n = 0, 0 do begin
    print, 'nlay_span('+strcompress(nlay_span)+'):  '+strcompress(n)
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
    dset ='/Users/shohei/tmp/MCD5.3/data/' ;'MCD_DATA/'  ; <character*50> data set
    scena = 1         ; <integer> scenario (1 = Climatology ave solar)
    perturkey = 1     ; <integer>  perturbation type (1= none)
    seedin   = 7.0    ; <real>
    gwlength = 0.0    ; <real>  for small scale (ie: gravity wave) perturbations;
    extvarkeys = LONARR(101) ; add 1 element because indexes in IDL start at 0
    for i0 = 0, 100 do extvarkeys[i0] = 1 ; <integer> array output type (extvar(i) = 0 : don't compute)
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

;    cd, '/Users/Shohei/Documents/Programs/MCD_ver4.3/mcd/idl/'
;   mcd_idl, Ls, Loct, lat, lon, dust, zkey, hrkey, xz, meanvarz, extvarz
    cd, '/Users/shohei/tmp/MCD5.3/mcd/idl/'
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
;    print, TB
;    ;-----------------

    for i = 0, count_nscan-1 do begin ;loop for slit scan
    ;for i = 0, 0 do begin ;loop for slit scan
      ;print, 'count_nscan('+strcompress(count_nscan)+'):  '+strcompress(i)
      for j = 0, ip-1 do begin  
      
        X_local = [wvl(113), wvl(114), wvl(115), wvl(116), wvl(117), wvl(118), wvl(119), wvl(120), wvl(122), wvl(123)]
        Y_local = [jdat(j,113,points(i)), jdat(j,114,points(i)), jdat(j,115,points(i)), jdat(j,116,points(i)), jdat(j,117,points(i)), $
                  jdat(j,118,points(i)), jdat(j,119,points(i)), jdat(j,120,points(i)), jdat(j,122,points(i)), jdat(j,123,points(i))]

        X_cont = [wvl(113), wvl(114), wvl(115), wvl(123)]
        Y_cont = [jdat(j,113,points(i)), jdat(j,114,points(i)), jdat(j,115,points(i)), jdat(j,123,points(i))]

        coef = linfit(X_cont, Y_cont)
        cont = coef(0) + coef(1)*x_local

        SZA = reform(geocube(j,8,points(i)))*1.e-4
        EA = reform(geocube(j,9,points(i)))*1.e-4
        PA = reform(geocube(j,10,points(i)))*1.e-4

        Albedo_input = jdat(j,113,points(i))/specmars(113) / cos(geocube(j,8,points(i))*1e-4*!DTOR)
        trans(j,points(i)) = total(1.0 - Y_local/cont)
        ;water(j,points(i)) = ret_water(trans(j,points(i)), SP, TA, TB, SZA, EA, PA, dust_opacity, Albedo_input)

        ;-------
          trans_omega = trans(j,points(i))
          Dust = dust_opacity
          Albedo =  Albedo_input
          T1 = TA
          T2 = TB
        
        
          ;result
;          Water_CD = -999d
          Error_U = -999d
          Error_L = -999d
        
          ;skip
          if SP gt 1500. then begin
            print,'Warning: SP > limit'
            goto, skip
          endif
          if SP lt 50. then begin
            print,'Warning: SP < limit'
            goto, skip
          endif
          if T1 gt 285. then begin
            print,'Warning: T1 > limit'
            goto, skip
          endif
          if T1 lt 160. then begin
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
          if SZA gt 60. then begin
            print,SZA
            print,'Warning: SZA > limit'
            goto, skip
          endif
          if SZA lt 0. then begin
            print,'Warning: SZA < limit'
            goto, skip
          endif
          if EA ge 10. then begin
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
;          if Waterice ne 1.0 then begin
;            print,'Warning: ICE > limit'
;            stop
;            goto, skip
;          endif
          if Albedo gt 0.5 then begin
            print,'Warning: Albedo > limit'
            goto, skip
          endif
          if Albedo lt 0.05 then begin
            print,'Warning: Albedo < limit'
            goto, skip
          endif
        
          x1a = dblarr(15) ; SP
          x1a(0) = alog(50d)
          x1a(1) = alog(150d)
          x1a(2) = alog(180d)
          x1a(3) = alog(215d)
          x1a(4) = alog(257d)
          x1a(5) = alog(308d)
          x1a(6) = alog(369d)
          x1a(7) = alog(442d)
          x1a(8) = alog(529d)
          x1a(9) = alog(633d)
          x1a(10) = alog(758d)
          x1a(11) = alog(907d)
          x1a(12) = alog(1096d)
          x1a(13) = alog(1300d)
          x1a(14) = alog(1500d)
        
          x2a = dblarr(4) ; T1
          x2a(0) = (135.d)^1.5d
          x2a(1) = (160.d)^1.5d
          x2a(2) = (213.d)^1.5d
          x2a(3) = (285.d)^1.5d
        
          x3a = dblarr(3) ; T2
          x3a(0) = (80.d)^1.5d
          x3a(1) = (146.d)^1.5d
          x3a(2) = (200.d)^1.5d
        
          x4a = dblarr(5) ; SZA
          x4a(0) = exp(-cos(0d/180d*!dpi))
          x4a(1) = exp(-cos(20d/180d*!dpi))
          x4a(2) = exp(-cos(40d/180d*!dpi))
          x4a(3) = exp(-cos(60d/180d*!dpi))
          x4a(4) = exp(-cos(80d/180d*!dpi))
        
          x5a = dblarr(3) ; EA
          x5a(0) = exp(-cos(0d/180d*!dpi))
          x5a(1) = exp(-cos(5d/180d*!dpi))
          x5a(2) = exp(-cos(10d/180d*!dpi))
        
          x6a = dblarr(4) ; PA
          x6a(0) = -cos(0d/180d*!dpi) + 1.d
          x6a(1) = -cos(60d/180d*!dpi) + 1.d
          x6a(2) = -cos(120d/180d*!dpi) + 1.d
          x6a(3) = -cos(180d/180d*!dpi) + 1.d
        
          x7a = dblarr(5) ; Dust
          x7a(0) = 0.0d
          x7a(1) = 0.3d
          x7a(2) = 0.6d
          x7a(3) = 0.9d
          x7a(4) = 1.2d
        
          ;x7a = dblarr(3) ; Water-ice
          ;x7a(0) = 0.0d
          ;x7a(1) = 0.5d
          ;x7a(2) = 1.0d
        
          x8a = dblarr(6) ; Surface Albedo
          x8a(0) = 0.05d
          x8a(1) = 0.1d
          x8a(2) = 0.2d
          x8a(3) = 0.3d
          x8a(4) = 0.4d
          x8a(5) = 0.5d
        
          ;water abundance grid
          Water_grid = dblarr(15)
          Water_grid(0) = 0.0*15.*1e-6
          Water_grid(1) = 5.0*15.*1e-6
          Water_grid(2) = 10.0*15.*1e-6
          Water_grid(3) = 10.0*15.*1e-6
          Water_grid(4) = 20.0*15.*1e-6
          Water_grid(5) = 25.0*15.*1e-6
          Water_grid(6) = 30.0*15.*1e-6
          Water_grid(7) = 35.0*15.*1e-6
          Water_grid(8) = 40.0*15.*1e-6
          Water_grid(9) = 45.0*15.*1e-6
          Water_grid(10) = 50.0*15.*1e-6
          Water_grid(11) = 60.0*15.*1e-6
          Water_grid(12) = 70.0*15.*1e-6
          Water_grid(13) = 80.0*15.*1e-6
          Water_grid(14) = 100.0*15.*1e-6
        
          y = dblarr(15)
          X1 = alog(SP) ;SP as test
          X2 = (T1)^1.5d ;T1 as test
          X3 = (T2)^1.5d ;T2 as test
          x4 = exp(-cos(SZA/180d*!dpi))
          x5 = exp(-cos(EA/180d*!dpi))
          x6 = -cos(PA/180d*!dpi) + 1.d
          x7 = Dust ;reduce_dust(loop)*1.86
          ;x7 = Waterice; reduce_waterice(loop)*3.87
          x8 = Albedo
        
          for ii = 0, 15-1 do begin
            F = x1a(ii) - X1
            if (F gt 0.0d0) then j1 = ii
            if (F gt 0.0d0) then goto, skip1
          endfor
          skip1:
        
          for ii = 0, 4-1 do begin
            F = x2a(ii) - X2
            if (F gt 0.0d0) then j2 = ii
            if (F gt 0.0d0) then goto, skip2
          endfor
          skip2:
        
          for ii = 0, 3-1 do begin
            F = x3a(ii) - X3
            if (F gt 0.0d0) then j3 = ii
            if (F gt 0.0d0) then goto, skip3
          endfor
          skip3:
        
          for ii = 0, 5-1 do begin
            F = x4a(ii) - X4
            if (F gt 0.0d0) then j4 = ii
            if (F gt 0.0d0) then goto, skip4
          endfor
          skip4:
        
          for ii = 0, 3-1 do begin
            F = x5a(ii) - X5
            if (F gt 0.0d0) then j5 = ii
            if (F gt 0.0d0) then goto, skip5
          endfor
          skip5:
        
          for ii = 0, 4-1 do begin
            F = x6a(ii) - X6
            if (F gt 0.0d0) then j6 = ii
            if (F gt 0.0d0) then goto, skip6
          endfor
          skip6:
        
          for ii = 0, 5-1 do begin
            F = x7a(ii) - X7
            if (F gt 0.0d0) then j7 = ii
            if (F gt 0.0d0) then goto, skip7
          endfor
          skip7:
        
          ;for I = 0, 3-1 do begin
          ;    F = x8a(I) - X8
          ;    if (F gt 0.0d0) then j8 = I
          ;    if (F gt 0.0d0) then goto, skip8
          ;endfor
          ;skip8:
          ;
          for ii = 0, 6-1 do begin
            F = x8a(ii) - X8
            if (F gt 0.0d0) then j8 = ii
            if (F gt 0.0d0) then goto, skip8
          endfor
          skip8:
        
          j1 = j1 -1
          j2 = j2 -1
          j3 = j3 -1
          j4 = j4 -1
          j5 = j5 -1
          j6 = j6 -1
          j7 = j7 -1
          j8 = j8 -1
        
          if j1 lt 0 then stop
          if j2 lt 0 then stop
          if j3 lt 0 then stop
          if j4 lt 0 then stop
          if j5 lt 0 then stop
          if j6 lt 0 then stop
          if j7 lt 0 then stop
          if j8 lt 0 then stop
          if j2+1 ge 4 then stop
          if j3+1 ge 3 then stop
          if j4+1 ge 5 then stop
          if j5+1 ge 3 then stop
          if j6+1 ge 4 then stop
          if j7+1 ge 5 then stop
          ;if j7+1 ge 4 then stop
          if j8+1 ge 6 then stop
        
          for ii = 0, 15-1 do begin
        
            ;8-dimensional_interpolation
            if ii eq 0 then Table_Equivalent_width = Table_Equivalent_width1
            if ii eq 1 then Table_Equivalent_width = Table_Equivalent_width2
            if ii eq 2 then Table_Equivalent_width = Table_Equivalent_width3
            if ii eq 3 then Table_Equivalent_width = Table_Equivalent_width4
            if ii eq 4 then Table_Equivalent_width = Table_Equivalent_width5
            if ii eq 5 then Table_Equivalent_width = Table_Equivalent_width6
            if ii eq 6 then Table_Equivalent_width = Table_Equivalent_width7
            if ii eq 7 then Table_Equivalent_width = Table_Equivalent_width8
            if ii eq 8 then Table_Equivalent_width = Table_Equivalent_width9
            if ii eq 9 then Table_Equivalent_width = Table_Equivalent_width10
            if ii eq 10 then Table_Equivalent_width = Table_Equivalent_width11
            if ii eq 11 then Table_Equivalent_width = Table_Equivalent_width12
            if ii eq 12 then Table_Equivalent_width = Table_Equivalent_width13
            if ii eq 13 then Table_Equivalent_width = Table_Equivalent_width14
            if ii eq 14 then Table_Equivalent_width = Table_Equivalent_width15
            y11111111=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7,j8)
            y21111111=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7,j8)
            y12111111=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7,j8)
            y11211111=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7,j8)
            y11121111=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7,j8)
            y11112111=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7,j8)
            y11111211=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7,j8)
            y22111111=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7,j8)
            y21211111=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7,j8)
            y21121111=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7,j8)
            y21112111=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7,j8)
            y21111211=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7,j8)
            y12211111=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7,j8)
            y12121111=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7,j8)
            y12112111=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7,j8)
            y12111211=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7,j8)
            y11221111=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7,j8)
            y11212111=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7,j8)
            y11211211=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7,j8)
            y11122111=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7,j8)
            y11121211=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7,j8)
            y11112211=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7,j8)
            y22211111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7,j8)
            y22121111=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7,j8)
            y22112111=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7,j8)
            y22111211=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7,j8)
            y21221111=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7,j8)
            y21212111=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7,j8)
            y21211211=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7,j8)
            y21122111=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7,j8)
            y21121211=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7,j8)
            y21112211=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7,j8)
            y12221111=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7,j8)
            y12212111=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7,j8)
            y12211211=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7,j8)
            y12122111=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7,j8)
            y12121211=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7,j8)
            y12112211=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7,j8)
            y11222111=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7,j8)
            y11221211=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7,j8)
            y11212211=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7,j8)
            y11122211=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7,j8)
            y11222211=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7,j8)
            y12122211=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7,j8)
            y12212211=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7,j8)
            y12221211=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7,j8)
            y12222111=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7,j8)
            y21122211=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7,j8)
            y21212211=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7,j8)
            y21221211=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7,j8)
            y21222111=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7,j8)
            y22112211=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7,j8)
            y22121211=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7,j8)
            y22122111=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7,j8)
            y22211211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7,j8)
            y22212111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7,j8)
            y22221111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7,j8)
            y12222211=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,j8)
            y21222211=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7,j8)
            y22122211=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7,j8)
            y22212211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7,j8)
            y22221211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7,j8)
            y22222111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7,j8)
            y22222211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,j8)
            y11111121=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7+1,j8)
            y21111121=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7+1,j8)
            y12111121=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7+1,j8)
            y11211121=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7+1,j8)
            y11121121=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7+1,j8)
            y11112121=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7+1,j8)
            y11111221=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7+1,j8)
            y22111121=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7+1,j8)
            y21211121=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7+1,j8)
            y21121121=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7+1,j8)
            y21112121=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7+1,j8)
            y21111221=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7+1,j8)
            y12211121=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7+1,j8)
            y12121121=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7+1,j8)
            y12112121=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7+1,j8)
            y12111221=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7+1,j8)
            y11221121=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7+1,j8)
            y11212121=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7+1,j8)
            y11211221=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7+1,j8)
            y11122121=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7+1,j8)
            y11121221=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7+1,j8)
            y11112221=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7+1,j8)
            y22211121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7+1,j8)
            y22121121=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7+1,j8)
            y22112121=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7+1,j8)
            y22111221=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7+1,j8)
            y21221121=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7+1,j8)
            y21212121=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7+1,j8)
            y21211221=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7+1,j8)
            y21122121=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7+1,j8)
            y21121221=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7+1,j8)
            y21112221=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7+1,j8)
            y12221121=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7+1,j8)
            y12212121=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7+1,j8)
            y12211221=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7+1,j8)
            y12122121=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7+1,j8)
            y12121221=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7+1,j8)
            y12112221=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7+1,j8)
            y11222121=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7+1,j8)
            y11221221=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7+1,j8)
            y11212221=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7+1,j8)
            y11122221=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7+1,j8)
            y11222221=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,j8)
            y12122221=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,j8)
            y12212221=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,j8)
            y12221221=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,j8)
            y12222121=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,j8)
            y21122221=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7+1,j8)
            y21212221=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7+1,j8)
            y21221221=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7+1,j8)
            y21222121=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7+1,j8)
            y22112221=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7+1,j8)
            y22121221=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7+1,j8)
            y22122121=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7+1,j8)
            y22211221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7+1,j8)
            y22212121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7+1,j8)
            y22221121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7+1,j8)
            y12222221=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,j8)
            y21222221=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,j8)
            y22122221=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,j8)
            y22212221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,j8)
            y22221221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,j8)
            y22222121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,j8)
            y22222221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,j8)
        
            y11111112=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7,j8+1)
            y21111112=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7,j8+1)
            y12111112=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7,j8+1)
            y11211112=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7,j8+1)
            y11121112=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7,j8+1)
            y11112112=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7,j8+1)
            y11111212=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7,j8+1)
            y22111112=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7,j8+1)
            y21211112=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7,j8+1)
            y21121112=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7,j8+1)
            y21112112=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7,j8+1)
            y21111212=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7,j8+1)
            y12211112=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7,j8+1)
            y12121112=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7,j8+1)
            y12112112=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7,j8+1)
            y12111212=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7,j8+1)
            y11221112=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7,j8+1)
            y11212112=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7,j8+1)
            y11211212=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7,j8+1)
            y11122112=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7,j8+1)
            y11121212=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7,j8+1)
            y11112212=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7,j8+1)
            y22211112=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7,j8+1)
            y22121112=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7,j8+1)
            y22112112=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7,j8+1)
            y22111212=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7,j8+1)
            y21221112=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7,j8+1)
            y21212112=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7,j8+1)
            y21211212=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7,j8+1)
            y21122112=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7,j8+1)
            y21121212=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7,j8+1)
            y21112212=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7,j8+1)
            y12221112=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7,j8+1)
            y12212112=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7,j8+1)
            y12211212=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7,j8+1)
            y12122112=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7,j8+1)
            y12121212=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7,j8+1)
            y12112212=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7,j8+1)
            y11222112=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7,j8+1)
            y11221212=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7,j8+1)
            y11212212=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7,j8+1)
            y11122212=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7,j8+1)
            y11222212=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7,j8+1)
            y12122212=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7,j8+1)
            y12212212=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7,j8+1)
            y12221212=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7,j8+1)
            y12222112=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7,j8+1)
            y21122212=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7,j8+1)
            y21212212=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7,j8+1)
            y21221212=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7,j8+1)
            y21222112=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7,j8+1)
            y22112212=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7,j8+1)
            y22121212=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7,j8+1)
            y22122112=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7,j8+1)
            y22211212=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7,j8+1)
            y22212112=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7,j8+1)
            y22221112=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7,j8+1)
            y12222212=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,j8+1)
            y21222212=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7,j8+1)
            y22122212=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7,j8+1)
            y22212212=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7,j8+1)
            y22221212=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7,j8+1)
            y22222112=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7,j8+1)
            y22222212=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,j8+1)
            y11111122=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7+1,j8+1)
            y21111122=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7+1,j8+1)
            y12111122=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7+1,j8+1)
            y11211122=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7+1,j8+1)
            y11121122=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7+1,j8+1)
            y11112122=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7+1,j8+1)
            y11111222=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7+1,j8+1)
            y22111122=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7+1,j8+1)
            y21211122=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7+1,j8+1)
            y21121122=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7+1,j8+1)
            y21112122=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7+1,j8+1)
            y21111222=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7+1,j8+1)
            y12211122=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7+1,j8+1)
            y12121122=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7+1,j8+1)
            y12112122=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7+1,j8+1)
            y12111222=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7+1,j8+1)
            y11221122=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7+1,j8+1)
            y11212122=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7+1,j8+1)
            y11211222=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7+1,j8+1)
            y11122122=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7+1,j8+1)
            y11121222=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7+1,j8+1)
            y11112222=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7+1,j8+1)
            y22211122=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7+1,j8+1)
            y22121122=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7+1,j8+1)
            y22112122=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7+1,j8+1)
            y22111222=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7+1,j8+1)
            y21221122=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7+1,j8+1)
            y21212122=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7+1,j8+1)
            y21211222=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7+1,j8+1)
            y21122122=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7+1,j8+1)
            y21121222=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7+1,j8+1)
            y21112222=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7+1,j8+1)
            y12221122=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7+1,j8+1)
            y12212122=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7+1,j8+1)
            y12211222=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7+1,j8+1)
            y12122122=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7+1,j8+1)
            y12121222=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7+1,j8+1)
            y12112222=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7+1,j8+1)
            y11222122=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7+1,j8+1)
            y11221222=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7+1,j8+1)
            y11212222=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7+1,j8+1)
            y11122222=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7+1,j8+1)
            y11222222=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,j8+1)
            y12122222=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,j8+1)
            y12212222=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,j8+1)
            y12221222=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,j8+1)
            y12222122=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,j8+1)
            y21122222=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7+1,j8+1)
            y21212222=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7+1,j8+1)
            y21221222=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7+1,j8+1)
            y21222122=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7+1,j8+1)
            y22112222=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7+1,j8+1)
            y22121222=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7+1,j8+1)
            y22122122=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7+1,j8+1)
            y22211222=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7+1,j8+1)
            y22212122=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7+1,j8+1)
            y22221122=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7+1,j8+1)
            y12222222=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,j8+1)
            y21222222=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,j8+1)
            y22122222=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,j8+1)
            y22212222=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,j8+1)
            y22221222=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,j8+1)
            y22222122=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,j8+1)
            y22222222=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,j8+1)
        
            t1 = (x1-x1a(j1)) / (x1a(j1+1)-x1a(j1))
            t2 = (x2-x2a(j2)) / (x2a(j2+1)-x2a(j2))
            t3 = (x3-x3a(j3)) / (x3a(j3+1)-x3a(j3))
            t4 = (x4-x4a(j4)) / (x4a(j4+1)-x4a(j4))
            t5 = (x5-x5a(j5)) / (x5a(j5+1)-x5a(j5))
            t6 = (x6-x6a(j6)) / (x6a(j6+1)-x6a(j6))
            t7 = (x7-x7a(j7)) / (x7a(j7+1)-x7a(j7))
            t8 = (x8-x8a(j7)) / (x8a(j8+1)-x8a(j8))
        
            y(ii) = y11111111*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y21111111*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y12111111*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y11211111*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y11121111*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y11112111*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y11111211*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y22111111*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y21211111*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y21121111*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y21112111*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y21111211*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y12211111*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y12121111*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y12112111*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y12111211*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y11221111*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y11212111*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y11211211*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y11122111*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y11121211*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y11112211*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y22211111*t1*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y22121111*t1*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y22112111*t1*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y22111211*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y21221111*t1*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y21212111*t1*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y21211211*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y21122111*t1*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y21121211*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y21112211*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y12221111*(1.d - t1)*t2*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y12212111*(1.d - t1)*t2*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y12211211*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y12122111*(1.d - t1)*t2*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y12121211*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y12112211*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y11222111*(1.d - t1)*(1.d - t2)*t3*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y11221211*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y11212211*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y11122211*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y11222211*(1.d - t1)*(1.d - t2)*t3*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y12122211*(1.d - t1)*t2*(1.d - t3)*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y12212211*(1.d - t1)*t2*t3*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y12221211*(1.d - t1)*t2*t3*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y12222111*(1.d - t1)*t2*t3*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y21122211*t1*(1.d - t2)*(1.d - t3)*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y21212211*t1*(1.d - t2)*t3*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y21221211*t1*(1.d - t2)*t3*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y21222111*t1*(1.d - t2)*t3*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y22112211*t1*t2*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y22121211*t1*t2*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y22122111*t1*t2*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y22211211*t1*t2*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y22212111*t1*t2*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y22221111*t1*t2*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y12222211*(1.d - t1)*t2*t3*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y21222211*t1*(1.d - t2)*t3*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y22122211*t1*t2*(1.d - t3)*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y22212211*t1*t2*t3*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y22221211*t1*t2*t3*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
              y22222111*t1*t2*t3*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
              y22222211*t1*t2*t3*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
              y11111121*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y21111121*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y12111121*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y11211121*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y11121121*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y11112121*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y11111221*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y22111121*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y21211121*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y21121121*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y21112121*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y21111221*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y12211121*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y12121121*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y12112121*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y12111221*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y11221121*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y11212121*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y11211221*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y11122121*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y11121221*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y11112221*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
              y22211121*t1*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y22121121*t1*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y22112121*t1*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y22111221*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y21221121*t1*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y21212121*t1*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y21211221*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y21122121*t1*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y21121221*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y21112221*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
              y12221121*(1.d - t1)*t2*t3*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y12212121*(1.d - t1)*t2*t3*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y12211221*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y12122121*(1.d - t1)*t2*(1.d - t3)*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y12121221*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y12112221*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
              y11222121*(1.d - t1)*(1.d - t2)*t3*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y11221221*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y11212221*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
              y11122221*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*t6*t7*(1.d - t8) +   $
              y11222221*(1.d - t1)*(1.d - t2)*t3*t4*t5*t6*t7*(1.d - t8) +   $
              y12122221*(1.d - t1)*t2*(1.d - t3)*t4*t5*t6*t7*(1.d - t8) +   $
              y12212221*(1.d - t1)*t2*t3*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
              y12221221*(1.d - t1)*t2*t3*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y12222121*(1.d - t1)*t2*t3*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y21122221*t1*(1.d - t2)*(1.d - t3)*t4*t5*t6*t7*(1.d - t8) +   $
              y21212221*t1*(1.d - t2)*t3*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
              y21221221*t1*(1.d - t2)*t3*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y21222121*t1*(1.d - t2)*t3*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y22112221*t1*t2*(1.d - t3)*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
              y22121221*t1*t2*(1.d - t3)*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y22122121*t1*t2*(1.d - t3)*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y22211221*t1*t2*t3*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y22212121*t1*t2*t3*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y22221121*t1*t2*t3*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
              y12222221*(1.d - t1)*t2*t3*t4*t5*t6*t7*(1.d - t8) +   $
              y21222221*t1*(1.d - t2)*t3*t4*t5*t6*t7*(1.d - t8) +   $
              y22122221*t1*t2*(1.d - t3)*t4*t5*t6*t7*(1.d - t8) +   $
              y22212221*t1*t2*t3*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
              y22221221*t1*t2*t3*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
              y22222121*t1*t2*t3*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
              y22222221*t1*t2*t3*t4*t5*t6*t7*(1.d - t8) +   $
              y11111112*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y21111112*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y12111112*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y11211112*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y11121112*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y11112112*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y11111212*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y22111112*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y21211112*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y21121112*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y21112112*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y21111212*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y12211112*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y12121112*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y12112112*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y12111212*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y11221112*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y11212112*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y11211212*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y11122112*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y11121212*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y11112212*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
              y22211112*t1*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y22121112*t1*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y22112112*t1*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y22111212*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y21221112*t1*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y21212112*t1*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y21211212*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y21122112*t1*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y21121212*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y21112212*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
              y12221112*(1.d - t1)*t2*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y12212112*(1.d - t1)*t2*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y12211212*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y12122112*(1.d - t1)*t2*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y12121212*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y12112212*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
              y11222112*(1.d - t1)*(1.d - t2)*t3*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y11221212*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y11212212*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
              y11122212*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*t6*(1.d - t7)*t8 +   $
              y11222212*(1.d - t1)*(1.d - t2)*t3*t4*t5*t6*(1.d - t7)*t8 +   $
              y12122212*(1.d - t1)*t2*(1.d - t3)*t4*t5*t6*(1.d - t7)*t8 +   $
              y12212212*(1.d - t1)*t2*t3*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
              y12221212*(1.d - t1)*t2*t3*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y12222112*(1.d - t1)*t2*t3*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y21122212*t1*(1.d - t2)*(1.d - t3)*t4*t5*t6*(1.d - t7)*t8 +   $
              y21212212*t1*(1.d - t2)*t3*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
              y21221212*t1*(1.d - t2)*t3*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y21222112*t1*(1.d - t2)*t3*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y22112212*t1*t2*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
              y22121212*t1*t2*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y22122112*t1*t2*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y22211212*t1*t2*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y22212112*t1*t2*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y22221112*t1*t2*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
              y12222212*(1.d - t1)*t2*t3*t4*t5*t6*(1.d - t7)*t8 +   $
              y21222212*t1*(1.d - t2)*t3*t4*t5*t6*(1.d - t7)*t8 +   $
              y22122212*t1*t2*(1.d - t3)*t4*t5*t6*(1.d - t7)*t8 +   $
              y22212212*t1*t2*t3*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
              y22221212*t1*t2*t3*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
              y22222112*t1*t2*t3*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
              y22222212*t1*t2*t3*t4*t5*t6*(1.d - t7)*t8 +   $
              y11111122*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y21111122*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y12111122*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y11211122*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y11121122*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y11112122*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
              y11111222*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
              y22111122*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y21211122*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y21121122*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y21112122*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
              y21111222*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
              y12211122*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y12121122*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y12112122*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
              y12111222*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
              y11221122*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y11212122*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
              y11211222*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
              y11122122*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*t7*t8 +   $
              y11121222*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*t7*t8 +   $
              y11112222*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*t7*t8 +   $
              y22211122*t1*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y22121122*t1*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y22112122*t1*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
              y22111222*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
              y21221122*t1*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y21212122*t1*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
              y21211222*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
              y21122122*t1*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*t7*t8 +   $
              y21121222*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*t7*t8 +   $
              y21112222*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*t7*t8 +   $
              y12221122*(1.d - t1)*t2*t3*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y12212122*(1.d - t1)*t2*t3*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
              y12211222*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
              y12122122*(1.d - t1)*t2*(1.d - t3)*t4*t5*(1.d - t6)*t7*t8 +   $
              y12121222*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*t6*t7*t8 +   $
              y12112222*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*t6*t7*t8 +   $
              y11222122*(1.d - t1)*(1.d - t2)*t3*t4*t5*(1.d - t6)*t7*t8 +   $
              y11221222*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*t6*t7*t8 +   $
              y11212222*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*t6*t7*t8 +   $
              y11122222*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*t6*t7*t8 +   $
              y11222222*(1.d - t1)*(1.d - t2)*t3*t4*t5*t6*t7*t8 +   $
              y12122222*(1.d - t1)*t2*(1.d - t3)*t4*t5*t6*t7*t8 +   $
              y12212222*(1.d - t1)*t2*t3*(1.d - t4)*t5*t6*t7*t8 +   $
              y12221222*(1.d - t1)*t2*t3*t4*(1.d - t5)*t6*t7*t8 +   $
              y12222122*(1.d - t1)*t2*t3*t4*t5*(1.d - t6)*t7*t8 +   $
              y21122222*t1*(1.d - t2)*(1.d - t3)*t4*t5*t6*t7*t8 +   $
              y21212222*t1*(1.d - t2)*t3*(1.d - t4)*t5*t6*t7*t8 +   $
              y21221222*t1*(1.d - t2)*t3*t4*(1.d - t5)*t6*t7*t8 +   $
              y21222122*t1*(1.d - t2)*t3*t4*t5*(1.d - t6)*t7*t8 +   $
              y22112222*t1*t2*(1.d - t3)*(1.d - t4)*t5*t6*t7*t8 +   $
              y22121222*t1*t2*(1.d - t3)*t4*(1.d - t5)*t6*t7*t8 +   $
              y22122122*t1*t2*(1.d - t3)*t4*t5*(1.d - t6)*t7*t8 +   $
              y22211222*t1*t2*t3*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
              y22212122*t1*t2*t3*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
              y22221122*t1*t2*t3*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
              y12222222*(1.d - t1)*t2*t3*t4*t5*t6*t7*t8 +   $
              y21222222*t1*(1.d - t2)*t3*t4*t5*t6*t7*t8 +   $
              y22122222*t1*t2*(1.d - t3)*t4*t5*t6*t7*t8 +   $
              y22212222*t1*t2*t3*(1.d - t4)*t5*t6*t7*t8 +   $
              y22221222*t1*t2*t3*t4*(1.d - t5)*t6*t7*t8 +   $
              y22222122*t1*t2*t3*t4*t5*(1.d - t6)*t7*t8 +   $
              y22222222*t1*t2*t3*t4*t5*t6*t7*t8
              ;----------------------------------------------------------------------------------
              bad = where(Table_Equivalent_width(j1:j1+1,j2:j2+1,j3:j3+1,j4:j4+1,j5:j5+1,j6:j6+1,j7:j7+1,j8:j8+1) eq 0, count)
              if count ge 1 then Y(ii) = -0d/0d
;              print, float(256d - count)/256d*100d
              ;----------------------------------------------------------------------------------
        
        
          endfor
        
          water_VMR(j,points(i)) = interpol(Water_grid,Y,trans_omega,/nan)

          H = R*TA/g
          P_layer = SP * exp(-H_layer/(H))
          H1 = H*0.1d+00
          H2 = H*4.0d+00
          a = (TA - TB) / (H1 - H2)
          b = TA - a*H1
          For ii = 0, 50 do begin
            T_layer(ii) = a*H_layer(ii) + b
            IF(H_layer(ii) ge H*4.0) then T_layer(ii) = TB
          endfor
          water_CD(j,points(i)) = 0d
          for ii = 0, 50 do water_CD(j,points(i)) = water_CD(j,points(i)) + P_layer(ii)/(T_layer(ii) * Kb) * 1d-6 * 1d5
          water_CD(j,points(i)) = water_CD(j,points(i)) * water_VMR(j,points(i)) / 3.34d18
;          print,water_VMR(j,points(i))*1d6
;          print,water_CD(j,points(i))
          skip:
          
      endfor
    endfor
  endfor

;      Set_Plot, 'Z'
;      Device, Set_Resolution=[1000,800], Set_Pixel_Depth=24, Decomposed=0

;      if loop eq 0 then restore, path+'work_ORB1023_3.sav'
;      if loop eq 1 then restore, path+'work_ORB1023_4.sav'
;      if loop eq 2 then restore, path+'work_ORB1023_5.sav'
;      if loop eq 3 then restore, path+'work_ORB1023_6.sav'
;      if loop eq 4 then restore, path+'work_ORB1023_7.sav'
;      if loop eq 5 then restore, path+'work_ORB1983_1.sav'

      loadct, 39
      Set_Plot, 'x'
      device, retain=1, decomposed=0
      !p.multi=[0,0,0]

      A = FINDGEN(17) * (!PI*2/16.)
      USERSYM, COS(A), SIN(A), /FILL
      
      !Y.OMargin = [3,0]
      window,2,xs=800,ys=800
      ;!p.multi=[0,2,1]
      range = max([max(longi(*,*))-min(longi(*,*)),max(lati(*,*))-min(lati(*,*))])
      min_water = 0;min(trans(*,*),/nan);0.
      max_water = max(trans(*,*),/nan)
      rn_water = max_water - min_water
      plot,longi(*,*),lati(*,*),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
        title='Water Depth: '+'Ls:'+STRCOMPRESS(Solar_longitude),charsize=2;.,$
;        xr=[min(min(longi(*,*))),min(longi(*,*))+range],yr=[min(min(lati(*,*))),min(lati(*,*))+range]
      for i = 0, io-1 do begin
        for j = 0, ip-1 do begin
          color = (trans(j,i)-min_water)/rn_water*254.
          if color gt 254 then color = 254
          if color lt 0 then color = 0
          plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
        endfor
      endfor
      colorbar_u2,format='(f4.2)',charsize=1.5,charthick=1,range=[0,max_water],DIVISIONS=4,position=[0.8/21., 4/29.7, 1.25/21., 28.2/29.7],title='Line depth',color=0
      snapshot = TVRD(True=1)
      Write_JPEG, path+file_basename(file,'.sav')+'_'+STRCOMPRESS(n, /remove_all)+'_Depth.jpg', snapshot, True=1, Quality=75
      
      window,2,xs=800,ys=800
      min_water = 0
;      max_water = 1000
      max_water = 1200;max(water_VMR(*,*)*1d6,/nan)
      rn_water = max_water - min_water
      plot,longi,lati,xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
        title='Water VMR: '+'Ls:'+STRCOMPRESS(Solar_longitude),charsize=2;.,$
;        xr=[min(min(longi(*,*))),min(longi(*,*))+range],yr=[min(min(lati(*,*))),min(lati(*,*))+range]
      for i = 0, io-1 do begin
        for j = 0, ip-1 do begin
          color = (water_VMR(j,i)*1d6-min_water)/rn_water*254.
          if color gt 254 then color = 254
          if FINITE(water_VMR(j,i)) eq 0 then color = 255
          if color lt 0 then color = 0
          plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
        endfor
      endfor
      colorbar_u2,format='(f5.0)',charsize=1.5,charthick=1,range=[0,max_water],DIVISIONS=4,position=[0.8/21., 4/29.7, 1.25/21., 28.2/29.7],title='ppm',color=0
      snapshot = TVRD(True=1)
      Write_JPEG, path+file_basename(file,'.sav')+'_'+STRCOMPRESS(n, /remove_all)+'_VMR.jpg', snapshot, True=1, Quality=75

      window,2,xs=800,ys=800
      min_water = 0
;      max_water = 80
      max_water = 80.;max(water_CD(*,*),/nan)
      rn_water = max_water - min_water
      plot,longi,lati,xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
        title='Water CD: '+'Ls:'+STRCOMPRESS(Solar_longitude),charsize=2;.,$
;        xr=[min(min(longi(*,*))),min(longi(*,*))+range],yr=[min(min(lati(*,*))),min(lati(*,*))+range]
      for i = 0, io-1 do begin
        for j = 0, ip-1 do begin
          color = (water_CD(j,i)-min_water)/rn_water*254.
          if color gt 254 then color = 254
          if FINITE(water_CD(j,i)) eq 0 then color = 255
          if color lt 0 then color = 0
          plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
          
          if FINITE(water_CD(j,i)) ne 0 then begin
            Water_CD_tot(nn) =  water_CD(j,i)
            lat_tot(nn) =  lati(j,i)
            nn = nn + 1l
          endif
        endfor
      endfor
      colorbar_u2,format='(f5.1)',charsize=1.5,charthick=1,range=[0,max_water],DIVISIONS=4,position=[0.8/21., 4/29.7, 1.25/21., 28.2/29.7],title='Pr-um',color=0
      snapshot = TVRD(True=1)
      Write_JPEG, path+file_basename(file,'.sav')+'_'+STRCOMPRESS(n, /remove_all)+'_CD.jpg', snapshot, True=1, Quality=75

      good = where(FINITE(water_CD) ne 0, count)
;      print, count
;      sza = reform(geocube(*,8,*)*1.e-4)
;      ea = reform(geocube(*,9,*)*1.e-4)
;      Albedo = jdat(*,113,*)/specmars(113) / cos(geocube(*,8,*)*1e-4*!DTOR)
;      bad = FINITE(water_CD) eq 0 
;      a = where(sza(bad) ge 60, count)
;      print,count
;      b = where(ea(bad) ge 10, count)
;      print, count
;      c = where(Albedo(bad) le 0.05, count)
;      print, count
;      d = where(Albedo(bad) ge 0.5, count)
;      print, count
;      stop

;      if loop eq 0 then file = path+'ORB1023_3.sav'
;      if loop eq 1 then file = path+'ORB1023_4.sav'
;      if loop eq 2 then file = path+'ORB1023_5.sav'
;      if loop eq 3 then file = path+'ORB1023_6.sav'
;      if loop eq 4 then file = path+'ORB1023_7.sav'
;      if loop eq 5 then file = path+'ORB0018_4.sav'
;      if loop eq 6 then file = path+'ORB0024_4.sav'
;      if loop eq 7 then file = path+'ORB0049_4.sav'
;      if loop eq 8 then file = path+'ORB0171_3.sav'

      if loop eq 0 then save,filename=path+'work_ORB1023_3.sav',/all
      if loop eq 1 then save,filename=path+'work_ORB1023_4.sav',/all
      if loop eq 2 then save,filename=path+'work_ORB1023_5.sav',/all
      if loop eq 3 then save,filename=path+'work_ORB1023_6.sav',/all
      if loop eq 4 then save,filename=path+'work_ORB1023_7.sav',/all
      if loop eq 5 then save,filename=path+'work_ORB1983_1.sav',/all
;      if loop eq 5 then save,filename=path+'work_ORB0018_4.sav',/all
;      if loop eq 6 then save,filename=path+'work_ORB0024_4.sav',/all
;      if loop eq 7 then save,filename=path+'work_ORB0049_4.sav',/all
;      if loop eq 8 then save,filename=path+'work_ORB0171_3.sav',/all
;      stop
endfor

final_lat = dindgen(181) - 90d
final_water = dblarr(181)
for i = 0, 180 do begin
  good = where(lat_tot ge final_lat(i)-1. and lat_tot le final_lat(i)+1., count)
  if count gt 1 then begin
    final_water(i) = mean(Water_CD_tot(good))
  endif
endfor

good = where(final_water ne 0)
window,0
plot, findgen(10), back=255, color=0, thick=3, xs=1, ys=1, xr=[-90,90], yr=[0,80], /nodata
oplot, lat_tot, Water_CD_tot, psym=8, symsize=0.2, color=0
oplot, final_lat(good), final_water(good), color=254, thick=3
stop
end
