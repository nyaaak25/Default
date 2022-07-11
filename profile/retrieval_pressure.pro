Pro retrieval_pressure

;++++++++++++++++++++++++
;Written by SHOHEI AOKI
;++++++++++++++++++++++++

;================================================
;path
;================================================
path = '/data2/omega/sav/'
restore, path + 'specmars.sav'

;restore
restore,'/work1/LUT/SP/table/absorption/Table_SP_Trans_calc.sav'

;para
kb = 1.38d-23
R = 192.0d+00
g = 3.72d+00
P_layer = dblarr(51)
H_layer = dindgen(51)*1d3
T_layer = dblarr(51)
nn = 0l

;================================================
;Loop start
;================================================
for Loop = 0, 2 do begin
  
  if loop eq 0 then file = path+'ORB0920_3.sav'
  if loop eq 1 then file = path+'ORB0931_3.sav'
  if loop eq 2 then file = path+'ORB0313_4.sav'

  restore, file
  ip = n_elements(LATI(*,0))
  io = n_elements(LATI(0,*))
  trans = reform(jdat(*,0,*))
  pressure = trans
  
  CO2=where(wvl gt 1.81 and wvl lt 2.19)
  wvl=wvl[CO2]
  jdat=jdat(*,CO2,*)
  specmars = specmars(CO2)
  band=where(wvl gt 1.9 and wvl lt 2.1)

  x = [wvl(0),wvl(3),wvl(5),wvl(23),wvl(24),wvl(25)]

  
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
      ;    print, TB
      ;    ;-----------------

    for i = 0, count_nscan-1 do begin ;loop for slit scan
      for j = 0, ip-1 do begin  
      

        Y = [jdat(j,0,points(i)), jdat(j,3,points(i)), jdat(j,5,points(i)), jdat(j,23,points(i)), jdat(j,24,points(i)), jdat(j,25,points(i))]
        coef = linfit(X,Y)
        cont = coef(0) + coef(1)*wvl

        SZA = reform(geocube(j,8,points(i)))*1.e-4
        EA = reform(geocube(j,9,points(i)))*1.e-4
        PA = reform(geocube(j,10,points(i)))*1.e-4
        
        Albedo_input = jdat(j,0,points(i))/specmars(0) / cos(geocube(j,8,points(i))*1e-4*!DTOR)
        width = 1.0 - jdat(j,*,points(i))/cont
        trans(j,points(i)) = total(width[band], /nan)
        pressure(j,points(i)) = ret_pressure(trans(j,points(i)), TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input)



        ;-------
          trans_omega = trans(j,points(i))
          Dust = dust_opacity
          Albedo =  Albedo_input
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
          if EA ge 10. then begin
            print,'Warning: EA > limit'
            goto, skip
          endif
          if PA gt 180. then begin
            print,'Warning: PA > limit'
            goto, skip
          endif
          if Dust gt 1.2 then begin
            print,'Warning: DUST > limit'
            goto, skip
          endif
          if Waterice gt 1.0 then begin
            print,'Warning: ICE > limit'
            goto, skip
          endif
          if Albedo gt 0.5 then begin
            print,'Warning: Albedo > limit'
            goto, skip
          endif
          if Albedo lt 0.1 then begin
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
          
          x4a = dblarr(3) ; EA
          x4a(0) = exp(-cos(0d/180d*!dpi))
          x4a(1) = exp(-cos(5d/180d*!dpi))
          x4a(2) = exp(-cos(10d/180d*!dpi))
          
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
          
          x8a = dblarr(6) ; Surface Albedo
          x8a(0) = 0.05d
          x8a(1) = 0.1d
          x8a(2) = 0.2d
          x8a(3) = 0.3d
          x8a(4) = 0.4d
          x8a(5) = 0.5d
          
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
          
          for I = 0, 3-1 do begin
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
          
          for I = 0, 6-1 do begin
            F = x8a(I) - X8
            if (F gt 0.0d0) then j8 = I
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
          if j1+1 ge 4 then stop
          if j2+1 ge 3 then stop
          if j3+1 ge 5 then stop
          if j4+1 ge 3 then stop
          if j5+1 ge 4 then stop
          if j6+1 ge 5 then stop
          if j7+1 ge 4 then stop
          if j8+1 ge 5 then stop
          
          for I = 0, 15-1 do begin
          
              ;8-dimensional_interpolation
              if i eq 0 then Table_Equivalent_width = Table_Equivalent_Pressure1
              if i eq 1 then Table_Equivalent_width = Table_Equivalent_Pressure2
              if i eq 2 then Table_Equivalent_width = Table_Equivalent_Pressure3
              if i eq 3 then Table_Equivalent_width = Table_Equivalent_Pressure4
              if i eq 4 then Table_Equivalent_width = Table_Equivalent_Pressure5
              if i eq 5 then Table_Equivalent_width = Table_Equivalent_Pressure6
              if i eq 6 then Table_Equivalent_width = Table_Equivalent_Pressure7
              if i eq 7 then Table_Equivalent_width = Table_Equivalent_Pressure8
              if i eq 8 then Table_Equivalent_width = Table_Equivalent_Pressure9
              if i eq 9 then Table_Equivalent_width = Table_Equivalent_Pressure10
              if i eq 10 then Table_Equivalent_width = Table_Equivalent_Pressure11
              if i eq 11 then Table_Equivalent_width = Table_Equivalent_Pressure12
              if i eq 12 then Table_Equivalent_width = Table_Equivalent_Pressure13
              if i eq 13 then Table_Equivalent_width = Table_Equivalent_Pressure14
              if i eq 14 then Table_Equivalent_width = Table_Equivalent_Pressure15
                          
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
              
              y(i) = y11111111*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
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
        
          pressure(j,points(i)) = interpol(Pressure_grid,Y,trans_omega,/nan)

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
          Pressure_CD(j,points(i)) = 0d
          for ii = 0, 50 do Pressure_CD(j,points(i)) = Pressure_CD(j,points(i)) + P_layer(ii)/(T_layer(ii) * Kb) * 1d-6 * 1d5
          Pressure_CD(j,points(i)) = Pressure_CD(j,points(i)) * Pressure_VMR(j,points(i)) / 3.34d18
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
