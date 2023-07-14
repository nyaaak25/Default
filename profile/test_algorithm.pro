;------------------------------------------------------------------------------------------------------------------------
function WHERE_XYZ, Array_expression, Count, XIND=xind, YIND=yind, ZIND=zind
  ;------------------------------------------------------------------------------------------------------------------------
  ; works for 1, 2 or 3 dimensional arrays
  ;
  ; Returns the 1D indices (same as WHERE)
  ;
  ; ARGUMENTS
  ;  - same as WHERE (see WHERE)
  ;
  ; KEYWORDS
  ; - Optionally returns X, Y, Z locations through:
  ;
  ; XIND: Output keyword, array of locations
  ; along the first dimension
  ; YIND: Output keyword, array of locations
  ; along the second dimension (if present)
  ; ZIND: Output keyword, array of locations
  ; along the third dimension (if present)
  ;
  ; If no matches where found, then XIND returns -1
  ;
  index_array=where(Array_expression, Count, /L64)
  dims=size(Array_expression,/dim)
  xind=index_array mod dims[0]
  case n_elements(dims) of
    2: yind=index_array / dims[0]
    3: begin
      yind=index_array / dims[0] mod dims[1]
      zind=index_array / dims[0] / dims[1]
    end
    else:
  endcase
  return, index_array
end

function BT,B,lambda
  C1 = 1.1911D-08
  C2 = 1.4387D+00
  T = 1.d + C1*(lambda^3) / B * 1e3
  T = 1.d / alog(T)
  T = T * C2 * lambda
  return,T
end

function check_continous, index

n = size(index)
n = fix(n(1))
flag = 0

for i = 0, n-3 do begin
  if (index(i) ne 0) and abs(index(i+1) - index(i)) eq 1 and abs(index(i+2) - index(i+1)) eq 1 then flag = 1
endfor
return, flag
end

Pro test_algorithm

;++++++++++++++++++++++++
;Written by SHOHEI AOKI
;++++++++++++++++++++++++

;================================================
;path
;================================================
path = '/Users/Shohei/tmp/test_case/'
restore, path + 'specmars.sav'

;================================================
;search file
;================================================
files = file_search(path+'*.sav',count=count)

;================================================
;Loop start
;================================================
for Loop = 0, count-1 do begin
  
  file = files(Loop)
  print, file
  restore, file
  ip = n_elements(LATI(*,0))
  io = n_elements(LATI(0,*))
  
  st1 = reform(bt(1d7*jdat(*,251,*)/(1d4/wvl(251))^2d,1d4/wvl(251)))
  st2 = reform(bt(1d7*jdat(*,252,*)/(1d4/wvl(252))^2d,1d4/wvl(252)))
  st3 = reform(bt(1d7*jdat(*,253,*)/(1d4/wvl(253))^2d,1d4/wvl(253)))
  st4 = reform(bt(1d7*jdat(*,254,*)/(1d4/wvl(254))^2d,1d4/wvl(254)))
  st5 = reform(bt(1d7*jdat(*,255,*)/(1d4/wvl(255))^2d,1d4/wvl(255)))

  albedo = reform(jdat(*,26,*)/specmars(26))
  rad_2777 = reform(jdat(*,140,*))
  ice_index = reform( (jdat(*,170,*)/specmars(170)) / (jdat(*,176,*)/specmars(176)) )

  ;wvl(115) 2.54 um
  ;wvl(120) 2.60 um
  ;wvl(123) 2.64 um
  trans = reform(jdat(*,0,*))
  trans_norm = reform(jdat(*,0,*))
  CO = reform(jdat(*,0,*))
  CO2 = reform(jdat(*,0,*))
  st = reform(jdat(*,0,*))
  LST = reform(jdat(*,0,*))
  R_surface_CO2ice = reform(jdat(*,0,*))
  h = 10. ;scale height
  flag = intarr(io)
  EM = reform(geocube(*, 9, *))*1.e-4
  IA = reform(geocube(*, 8, *))*1.e-4
  alt = reform(geocube(*,12,*)*1.e-3)
  flag(*) = 0
  X = [wvl(115), wvl(123)]
  ;X = [wvl(114), wvl(115), wvl(123), wvl(124)]
  X_CO = [wvl(94), wvl(108)]
  X_CO2 = [wvl(61), wvl(62), wvl(63), wvl(85), wvl(86), wvl(87)]
  CO2_wl = where(findgen(100) ge 61 and findgen(100) le 87 and findgen(100) ne 78)

  ;CO2 ice on the ground
   ;1.385 um = wvl(32)
   ;1.429 um = wvl(35)
   ;1.443 um = wvl(36)
  

  for i = 0, io-1 do begin
      for j = 0, ip-1 do begin 
        Y = [jdat(j,115,i), jdat(j,123,i)]
        ;Y = [jdat(j,114,i), jdat(j,115,i), jdat(j,123,i), jdat(j,124,i)]
        Y_CO = [jdat(j,94,i), jdat(j,108,i)]
        Y_CO2 = [jdat(j,61,i)/specmars(61), jdat(j,62,i)/specmars(62), jdat(j,63,i)/specmars(63), $
                 jdat(j,85,i)/specmars(85), jdat(j,86,i)/specmars(86), jdat(j,87,i)/specmars(87)]
        coef = linfit(X,Y)
        coef_CO = linfit(X_CO,Y_CO)
        coef_CO2 = linfit(X_CO2,Y_CO2)
        cont = coef(0) + coef(1)*wvl(120)
        cont_CO = coef_CO(0) + coef_CO(1)*wvl(100)
        cont_CO2 = coef_CO2(0) + coef_CO2(1)*wvl(CO2_wl)
        CO2(j,i) = total(1.0 - (jdat(j,CO2_wl,i)/specmars(CO2_wl))/cont_CO2,/double,/nan)
        RF1429 = jdat(j,35,i)/specmars(35)*cos(geocube(j,8,i)*1e-4*!DTOR)
        RF1385 = jdat(j,32,i)/specmars(32)*cos(geocube(j,8,i)*1e-4*!DTOR)
        RF1443 = jdat(j,36,i)/specmars(36)*cos(geocube(j,8,i)*1e-4*!DTOR)
        R_surface_CO2ice(j,i) = RF1429/((RF1385^0.5d)*(RF1443^0.5d))
 
;        plot, wvl(CO2_wl), jdat(j,CO2_wl,i)/specmars(CO2_wl),psym=-1,xs=1
;        oplot, wvl(CO2_wl), cont_CO2, color=254
;        stop
;        
        trans(j,i) = 1.0 - jdat(j,120,i)/cont
        CO(j,i) = 1.0 - jdat(j,100,i)/cont_CO
;        airmass1 = 1.0 / cos(geocube(j,4,i)*1e-4*!DTOR)
;        airmass2 = 1.0 / cos(geocube(j,5,i)*1e-4*!DTOR)
        airmass1 = 1.0 / cos(geocube(j,8,i)*1e-4*!DTOR)
        ;airmass1 = 1.0 / (1.d - cos(geocube(j,8,i)*1e-4*!DTOR))
        airmass2 = 1.0 / cos(geocube(j,9,i)*1e-4*!DTOR)
        albedo(j,i) = albedo(j,i); /cos(geocube(j,8,i)*1e-4*!DTOR)
        trans(j,i) = trans(j,i); / (airmass1 + airmass2)
        CO(j,i) = CO(j,i); / (airmass1 + airmass2)
        CO2(j,i) = CO2(j,i) ;/ (airmass1 + airmass2)
        trans_norm(j,i) = trans(j,i) * exp(alt(j,i)/H)
        st_sets = [st1(j,i),st2(j,i),st3(j,i),st4(j,i),st5(j,i)]
        st(j,i) = median(st_sets)
        LST(j,i) = (-SUB_SOLAR_LONGITUDE + longi(j,i)) / 15. + 12.
        if EM(j,i) ge 60. then flag(i) = 1
        if IA(j,i) ge 90. then flag(i) = 1
        if JDAT(j,115,i) ge 2. then flag(i) = 1
        if JDAT(j,115,i) le 0.02 then flag(i) = 1
        if alt(j,i) ge 65.535 then flag(i) = 1
      endfor
  endfor
;  stop

  ;latitude grid
  span_lati = max(lati(0,*)) - min(lati(0,*))
  nlay_span = ceil(span_lati/5.)
    
  ;search H2O enhancement 
  flag_detection = 0
  for n = 0, nlay_span-1 do begin
    hw_line = intarr(1000)
    flag_hw_1 = 0 
    flag_hw_2 = 0
    i_keep = 0
    points = where(lati(0,*) ge min(lati(0,*)) + float(n)*5. and lati(0,*) lt min(lati(0,*)) + float(n+1)*5., count)
    count_nscan = count
    base_value = median(trans(*,points))  
    stddev_value = stddev(trans(*,points),/nan)

    for i = 0, count_nscan-1 do begin ;loop for slit scan 

      ;skip the bad data
      if flag(points(i)) eq 1 then goto, skip0

      ;search H2O enhancement
      for j = 0, ip-1 do begin ;loop for pixel 

        if trans(j,points(i)) ge base_value + stddev_value*4. then begin
          flag_hw_1 = flag_hw_1 + 1
          if i ne i_keep then begin
            hw_line(flag_hw_2) = i
            flag_hw_2 = flag_hw_2 + 1
          endif
          i_keep = i
        endif

      endfor

      skip0:
    endfor

    flag_detection = check_continous(hw_line)
    iscan_detected = where(flag_detection eq 1, count)
    if count ge 1 then goto, detected
    if count le 1 then goto, nodetection
  
    detected:
    max_point = where_xyz(trans(iscan_detected,*) eq max(trans(iscan_detected,*),/nan), XIND=xind, YIND=yind)
    ;print, xind, yind
    max_lat = lati(iscan_detected(xind),yind)
    max_lon = longi(iscan_detected(xind),yind)
    max_Ls = Solar_longitude(iscan_detected(xind),yind)
    
      
      ;set up to create JPG file without X window
;      Set_Plot, 'Z'
;;      Device, Set_Resolution=[2500,2000], Set_Pixel_Depth=24, Decomposed=0
;      Device, Set_Resolution=[2500,1000], Set_Pixel_Depth=24, Decomposed=0
      loadct, 39
      Set_Plot, 'x'
      device, retain=1, decomposed=0
      window,0,xs=1250,ys=800
      !p.multi=[0,6,2]
      ;!p.multi=[0,5,2]
      !Y.OMargin = [5,0]
      
      min_water = min(trans(*,points),/nan);0.
      max_water = max(trans(*,points),/nan);0.1; max(trans(*,where(flag(points) eq 0)),/nan)*1.1
      rn_water = max_water - min_water
      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
        title='Water: '+'Ls:'+STRCOMPRESS(Solar_longitude)+' Max: '+string(max_water,format='(e0.1)' ),charsize=2.
      for i = 0, count_nscan-1 do begin
       if flag(points(i)) eq 1 then goto, skip1
        for j = 0, ip-1 do begin
          color = (trans(j,points(i))-min_water)/rn_water*254.
          if color gt 254 then color = 254
          if color lt 0 then color = 0
          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
        endfor
        skip1:
      endfor
      
      min_CO2 = min(CO2(*,points),/nan)
      max_CO2 = max(CO2(*,points),/nan)
      rn_CO2 = max_CO2 - min_CO2
      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
        title='CO2'+' Max: '+string(max_CO2,format='(e0.1)' )+' Min: '+string(min_CO2,format='(e0.1)' ),charsize=2.
      for i = 0, count_nscan-1 do begin
        if flag(points(i)) eq 1 then goto, skip4
        for j = 0, ip-1 do begin
          color = (CO2(j,points(i))-min_CO2)/rn_CO2*254.
          if color gt 254 then color = 254
          if color lt 0 then color = 0
          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
        endfor
        skip4:
      endfor
      
      min_CO = min(CO(*,points),/nan);0.
      max_CO = max(CO(*,points),/nan);0.05; max(CO(*,where(flag(points) eq 0)),/nan)*1.1
      rn_CO = max_CO - min_CO
      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
        title='CO'+' Max: '+string(max_CO,format='(e0.1)' ),charsize=2.
      for i = 0, count_nscan-1 do begin
        if flag(points(i)) eq 1 then goto, skip3
        for j = 0, ip-1 do begin
          color = (CO(j,points(i))-min_CO)/rn_CO*254.
          if color gt 254 then color = 254
          if color lt 0 then color = 0
          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
        endfor
        skip3:
      endfor
  
      max_albedo = max(albedo(*,points),/nan)
      min_albedo = min(albedo(*,points),/nan)
      rn_albedo = max_albedo - min_albedo
      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
        title='Albedo (I/F@1.3um): '+string(min_albedo,format='(f5.3)')+' - '+string(max_albedo,format='(f5.3)'),charsize=2.
      for i = 0, count_nscan-1 do begin
        for j = 0, ip-1 do begin
          color = (albedo(j,points(i))-min_albedo)/rn_albedo*254.
          if color gt 254 then color = 254
          if color lt 0 then color = 0
          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
        endfor
      endfor
  
      max_ice = max(ice_index(*,points));1.0
      min_ice = min(ice_index(*,points));0.5
      rn_ice = max_ice - min_ice
      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
        title='Ice index (3.40/3.52): '+string(min_ice,format='(F3.1)')+' - '+string(max_ice,format='(F3.1)'),charsize=2.
      for i = 0, count_nscan-1 do begin
        for j = 0, ip-1 do begin
          color = (ice_index(j,points(i))-min_ice)/rn_ice*254.
          if color gt 254 then color = 254
          if color lt 0 then color = 0
          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
        endfor
      endfor

      ;
;      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
;        title='MOLA altitude',charsize=2.
;      max_alt = max(alt(*,points),/nan)
;      min_alt = min(alt(*,points),/nan)
;      rn_alt = max_alt - min_alt
;      for i = 0, count_nscan-1 do begin
;        for j = 0, ip-1 do begin
;          color = (alt(j,points(i))-min_alt)/rn_alt*254.
;          if color lt 0 then color = 0
;          if color gt 254 then color = 254
;          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
;        endfor
;      endfor
;  
      max_CO2ice = max(R_surface_CO2ice(*,points));1.0
      min_CO2ice = min(R_surface_CO2ice(*,points));0.5
      rn_CO2ice = max_CO2ice - min_CO2ice
      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
        title='R_surface_CO2ice: '+string(min_CO2ice,format='(F3.1)')+' - '+string(max_CO2ice,format='(F3.1)'),charsize=2.
      for i = 0, count_nscan-1 do begin
        for j = 0, ip-1 do begin
          color = (R_surface_CO2ice(j,points(i))-min_CO2ice)/rn_CO2ice*254.
          if color gt 254 then color = 254
          if color lt 0 then color = 0
          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
        endfor
      endfor

;      max_rad_2777 = max(rad_2777(*,points),/nan)
;      min_rad_2777 = min(rad_2777(*,points),/nan)
;      rn_rad_2777 = max_rad_2777 - min_rad_2777
;      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
;        title='rad@2777 (dust): '+string(min_rad_2777,format='(e0.1)')+' - '+string(max_rad_2777,format='(e0.1)'),charsize=2.
;      for i = 0, count_nscan-1 do begin
;        for j = 0, ip-1 do begin
;          color = (rad_2777(j,points(i))-min_rad_2777)/rn_rad_2777*254.
;          if color gt 254 then color = 254
;          if color lt 0 then color = 0
;          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
;        endfor
;      endfor
;  
;      !Y.OMargin = [5,5]
;      max_st = max(st(*,points),/nan)
;      min_st = min(st(*,points),/nan)
;      rn_st = max_st - min_st
;      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
;        title='Surface Temperature: '+string(min_st,format='(f5.1)')+' - '+string(max_st,format='(f5.1)')+' K',charsize=2.
;      for i = 0, count_nscan-1 do begin
;        for j = 0, ip-1 do begin
;          color = (st(j,points(i))-min_st)/rn_st*254.
;          if color gt 254 then color = 254
;          if color lt 0 then color = 0
;          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
;        endfor
;      endfor

;
;  
;
;      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
;        title='Emission angle',charsize=2.
;      max_EM = max(EM(*,points))
;      min_EM = min(EM(*,points))
;      rn_EM = max_EM - min_EM
;      for i = 0, count_nscan-1 do begin
;        for j = 0, ip-1 do begin
;          color = (EM(j,points(i))-min_EM)/rn_EM*254.
;          if color lt 0 then color = 0
;          if color gt 254 then color = 254
;          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
;        endfor
;      endfor
;
;      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
;        title='Incidence angle',charsize=2.
;      max_IA = max(IA(*,points))
;      min_IA = min(IA(*,points))
;      rn_IA = max_IA - min_IA
;      for i = 0, count_nscan-1 do begin
;        for j = 0, ip-1 do begin
;          color = (IA(j,points(i))-min_IA)/rn_IA*254.
;          if color lt 0 then color = 0
;          if color gt 254 then color = 254
;          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
;        endfor
;      endfor

;      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
;        title='Local Solar Time',charsize=2.
;      max_LST = 12. + 6.
;      min_LST = 12. - 6. 
;      rn_LST = max_LST - min_LST
;      for i = 0, count_nscan-1 do begin
;        for j = 0, ip-1 do begin
;          color = (LST(j,points(i))-min_LST)/rn_LST*254.
;          if color lt 0 then color = 0
;          if color gt 254 then color = 254
;          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
;        endfor
;      endfor

;      max_lll = max(cos(geocube(*,8,points)*!dtor*1.e-4))
;      min_lll = min(cos(geocube(*,8,points)*!dtor*1.e-4))
;      rn_lll = max_lll - min_lll
;      plot,longi(*,points),lati(*,points),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
;        title='local_light_level',charsize=2.
;      for i = 0, count_nscan-1 do begin
;        for j = 0, ip-1 do begin
;          color = (cos(geocube(j,8,points(i))*!dtor*1.e-4)-min_lll)/rn_lll*254.
;          plots, longi(j,points(i)),lati(j,points(i)),color=color,psym=6,symsize=0.5,thick=3
;        endfor
;      endfor

;      COLORBAR, NCOLORS=254, POSITION=[0.05, 0.506, 0.17, 0.516], COLOR=0, MAXRANGE=max_water, MINRANGE=min_water, FORMAT='(F4.2)'
;      COLORBAR, NCOLORS=254, POSITION=[0.25, 0.506, 0.37, 0.516], COLOR=0, MAXRANGE=max_CO, MINRANGE=min_CO, FORMAT='(F4.2)'
;      COLORBAR, NCOLORS=254, POSITION=[0.45, 0.506, 0.57, 0.516], COLOR=0, MAXRANGE=max_albedo, MINRANGE=min_albedo, FORMAT='(F4.2)'
;      COLORBAR, NCOLORS=254, POSITION=[0.65, 0.506, 0.77, 0.516], COLOR=0, MAXRANGE=max_ice, MINRANGE=min_ice, FORMAT='(F4.2)'
;      COLORBAR, NCOLORS=254, POSITION=[0.85, 0.506, 0.97, 0.516], COLOR=0, MAXRANGE=max_rad_2777, MINRANGE=m_rad_2777, FORMAT='(F4.2)'
;
;      COLORBAR, NCOLORS=254, POSITION=[0.05, 0.020, 0.17, 0.030], COLOR=0, MAXRANGE=max_st, MINRANGE=min_st, FORMAT='(F5.1)'
;      COLORBAR, NCOLORS=254, POSITION=[0.25, 0.020, 0.37, 0.030], COLOR=0, MAXRANGE=max_alt, MINRANGE=min_alt, FORMAT='(F5.1)'
;      COLORBAR, NCOLORS=254, POSITION=[0.45, 0.020, 0.57, 0.030], COLOR=0, MAXRANGE=max_EM, MINRANGE=min_EM, FORMAT='(F5.1)'
;      COLORBAR, NCOLORS=254, POSITION=[0.65, 0.020, 0.77, 0.030], COLOR=0, MAXRANGE=max_IA, MINRANGE=min_IA, FORMAT='(F5.1)'
;      COLORBAR, NCOLORS=254, POSITION=[0.85, 0.020, 0.97, 0.030], COLOR=0, MAXRANGE=max_LST, MINRANGE=min_LST, FORMAT='(F5.1)'
;      COLORBAR, NCOLORS=254, POSITION=[0.02, 0.020, 0.16, 0.030], COLOR=0, MAXRANGE=max_water, MINRANGE=min_water, FORMAT='(F5.1)'
;      COLORBAR, NCOLORS=254, POSITION=[0.185, 0.020, 0.325, 0.030], COLOR=0, MAXRANGE=max_CO2, MINRANGE=min_CO2, FORMAT='(F5.1)'
;      COLORBAR, NCOLORS=254, POSITION=[0.35, 0.020, 0.49, 0.030], COLOR=0, MAXRANGE=max_CO, MINRANGE=min_CO, FORMAT='(F5.1)'
;      COLORBAR, NCOLORS=254, POSITION=[0.515, 0.020, 0.655, 0.030], COLOR=0, MAXRANGE=max_albedo, MINRANGE=min_albedo, FORMAT='(F5.1)'
;      COLORBAR, NCOLORS=254, POSITION=[0.68, 0.020, 0.82, 0.030], COLOR=0, MAXRANGE=max_ice, MINRANGE=min_ice, FORMAT='(F4.2)'
;      COLORBAR, NCOLORS=254, POSITION=[0.845, 0.020, 0.985, 0.030], COLOR=0, MAXRANGE=max_alt, MINRANGE=min_alt, FORMAT='(F5.1)'
;
      !Y.OMargin = [0,0]
 ;     !p.multi=0  
      
;      snapshot = TVRD(True=1)
;      Write_JPEG, path+file_basename(file,'.sav')+'_'+STRCOMPRESS(n, /remove_all)+'.jpg', snapshot, True=1, Quality=75
;      Write_JPEG, path+file_basename(file,'.sav')+'_'+STRCOMPRESS(n, /remove_all)+'_short.jpg', snapshot, True=1, Quality=75
      stop
      
;      Set_Plot, 'Z'
;      ;      Device, Set_Resolution=[2500,2000], Set_Pixel_Depth=24, Decomposed=0
;      Device, Set_Resolution=[2000,1000], Set_Pixel_Depth=24, Decomposed=0
;      loadct, 39
;      !p.multi=[0,1,2]
;
;;      Set_Plot, 'X'
;;      device, retain=1, decomposed=0
;;      window,0,xs=1200
;      plot, findgen(100), color=0, back=255, thick=3, ys=1, xs=1, xr=[1,2.65], yr=[0,.2]
;      for i = 0, count_nscan -1 do for j = 0, ip-1 do oplot, wvl(0:124), jdat(j,0:124,points(i))/specmars(0:124), color=float(i)/float(count_nscan)*254.
;
;;      window,1,xs=1200
;      plot, findgen(100), color=0, back=255, thick=3, ys=1, xs=1, xr=[2.5,2.7], yr=[0.5,1.2]
;      for i = 0, count_nscan -1 do begin
;        for j = 0, ip-1 do begin
;                X = [wvl(115), wvl(123)]
;                Y = [jdat(j,115,points(i)), jdat(j,123,points(i))]
;                coef = linfit(X,Y)
;          cont = coef(0) + coef(1)*wvl
;          oplot, wvl(0:124), jdat(j,0:124,points(i))/cont, color=float(i)/float(count_nscan)*254.
;        endfor
;      endfor
;
;      snapshot = TVRD(True=1)
;;      Write_JPEG, path+file_basename(file,'.sav')+'_'+STRCOMPRESS(n, /remove_all)+'_sp.jpg', snapshot, True=1, Quality=75
;;
;;      window,2,xs=500,ys=500
;;      plot,albedo(*,points),trans(*,points),psym=1,color=0,back=255;,yr=[0,0.2],xr=[0,0.5]
;      stop

      nodetection:
    endfor

endfor


 stop
end
