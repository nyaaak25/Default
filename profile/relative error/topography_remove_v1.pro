; ================================================
; EW法を検証するためのpro file
; 高度補正、2軌道の差を取る、ヒストグラムなど、、（TBD）
; ================================================

Pro topography_remove_v1

  ;おまじない
  device, decomposed = 0, retain = 2
  loadct,39

  ;restore file
  ; 2軌道の差を見る
  restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0920_3_LUT1.sav'
  ind = where_xyz(longi ge 272 and longi le 277 and lati ge 50 and lati le 61, xind=xind, yind=yind)
  TAMAP1 = TAMAP
  albedomap1 = inputalbedo
  Alt1 = altitude
  lat1 = lati
  lon1 = longi
  p1 = pressure
  p1 = exp(p1)

;  restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0931_3_LUT5.sav'
;  lat2 = lati
;  lon2 = longi
;  p2 = pressure
;  p2 = exp(p2)
;  p2_mod = p2
;  p2_mod(*,*) = -0d/0d

  restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0931_3_LUT1.sav'
  lat2 = lati
  lon2 = longi
  p2 = pressure
  p2 = exp(p2)
  p2_mod = p2
  p2_mod(*,*) = -0d/0d
  
  ; 1番近いピクセルを探してくる
  for i = 0, 127 do begin
    for j = 0, 595 do begin
      a = WHERE_XYZ(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*)) eq min(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*))), XIND=xind, YIND=yind, count)
      if count ge 1 then begin
        if abs(lat1(i,j)-lat2(xind(0),yind(0))) le 0.03 and abs(lon1(i,j)-lon2(xind(0),yind(0))) le 0.03 then begin
          p2_mod(i,j) = p2(xind(0),yind(0))
        endif
      endif
    endfor
  endfor
  
  stop


  ; プロットの枠組み
  A = FINDGEN(17) * (!PI*2/16.)
  USERSYM, COS(A), SIN(A), /FILL

  ; plot
  window, 1, xs=1200, ys=800
  !P.Multi = [0, 2, 2]
  
  ; SP1のプロット
  
  loadct,39
  plot, findgen(10),xs=1, ys=1, yr=[49.5, 61.5], xr=[271.8, 277], back=255, color=0, /nodata, charsize=2, title='Orb920 SPmap', xtitle='Lon', ytitle='Lat'
  ;plot, findgen(10),xs=1, ys=1, yr=[57.8, 60.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB920_3 SPmap', xtitle='Lon', ytitle='Lat',position=[0.15,0.55,0.5,0.9]
  ;plot, findgen(10),xs=1, ys=1, yr=[51.8, 54.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB920_3 SPmap', xtitle='Lon', ytitle='Lat',position=[0.15,0.55,0.5,0.9]
  ;plot, findgen(10),xs=1, ys=1, yr=[55.8, 58.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB920_3 SPmap', xtitle='Lon', ytitle='Lat',position=[0.15,0.55,0.5,0.9]
  ;plot, findgen(10),xs=1, ys=1, yr=[53.8, 56.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB920_3 Albedo map', xtitle='Lon', ytitle='Lat',position=[0.14,0.55,0.49,0.9]
  ;cgLOADCT, 16
  ;cgColorbar, Divisions=4, Minor=5, Range=[622,834],position=[0.058,0.57,0.065,0.87],title='SPmap', TCHARSIZE=15,/vertical
  ;cgColorbar, Divisions=4, Minor=5, Range=[1000,1600],position=[0.058,0.57,0.065,0.87],title='Albedo', TCHARSIZE=15,/vertical
  loadct, 16
  for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 and lat1(i,j) gt 50 and lat1(i,j) lt 61 and lon1(i,j) gt 272 and lon1(i,j) lt 276.8 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-622.)/212.*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 and lat1(i,j) gt 58 and lat1(i,j) lt 60 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-622.)/212.*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 and lat1(i,j) gt 52 and lat1(i,j) lt 54 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-622.)/212.*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 and lat1(i,j) gt 54 and lat1(i,j) lt 56 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-622.)/212.*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 and lat1(i,j) gt 54 and lat1(i,j) lt 56 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-0.1)/0.2*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 and lat1(i,j) gt 56 and lat1(i,j) lt 58 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-622.)/212.*254., psym=8, symsize=1

  ; SP2のプロット
  loadct, 39
  plot, findgen(10), xs=1, ys=1, yr=[49.5, 61.5], xr=[271.8, 277], back=255, color=0, /nodata, charsize=2, title='Orb931, SPmap', xtitle='Lon', ytitle='Lat
  ;plot, findgen(10),xs=1, ys=1, yr=[57.8, 60.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB931_3 SPmap', xtitle='Lon', ytitle='Lat',position=[0.6,0.55,0.95,0.9]
  ;plot, findgen(10),xs=1, ys=1, yr=[51.8, 54.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB931_3 SPmap', xtitle='Lon', ytitle='Lat',position=[0.6,0.55,0.95,0.9]
  ;plot, findgen(10),xs=1, ys=1, yr=[55.8, 58.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB931_3 SPmap', xtitle='Lon', ytitle='Lat',position=[0.6,0.55,0.95,0.9]
  ;plot, findgen(10),xs=1, ys=1, yr=[53.8, 56.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB931_3 Albedo map', xtitle='Lon', ytitle='Lat',position=[0.6,0.55,0.95,0.9]
  loadct, 16
  for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 50 and lat1(i,j) lt 61 and lon1(i,j) gt 272 and lon1(i,j) lt 276.8 then plots, lon1(i,j), lat1(i,j), color=(p2_mod(i,j)-622.)/212.*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 58 and lat1(i,j) lt 60 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p2_mod(i,j)-622.)/212.*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 52 and lat1(i,j) lt 54 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p2_mod(i,j)-622.)/212.*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 56 and lat1(i,j) lt 58 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p2_mod(i,j)-622.)/212.*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 54 and lat1(i,j) lt 56 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p2_mod(i,j)-0.1)/0.2*254., psym=8, symsize=1

  ; ORB920-0RB931を差し引いた図
  loadct, 39
  plot, findgen(10),  xs=1, ys=1,yr=[49.5, 61.5], xr=[271.8, 277], back=255, color=0, /nodata, charsize=2, title='ORB920 - ORB931', xtitle='Lon', ytitle='Lat'
  ;plot, findgen(10),  xs=1, ys=1,yr=[57.8, 60.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB920 - ORB931', xtitle='Lon', ytitle='Lat',position=[0.15,0.1,0.5,0.4]
  ;plot, findgen(10),  xs=1, ys=1,yr=[51.8, 54.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB920 - ORB931', xtitle='Lon', ytitle='Lat',position=[0.15,0.1,0.5,0.4]
  ;plot, findgen(10),  xs=1, ys=1,yr=[55.8, 58.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB920 - ORB931', xtitle='Lon', ytitle='Lat',position=[0.15,0.1,0.5,0.4]
  ;plot, findgen(10),  xs=1, ys=1,yr=[53.8, 56.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB920 - ORB931', xtitle='Lon', ytitle='Lat',position=[0.15,0.1,0.5,0.4]
  ;cgLOADCT, 33
  ;cgColorbar, Divisions=4, Minor=5, Range=[-15,15],position=[0.058,0.07,0.065,0.37],title='Pressure (Pa)', TCHARSIZE=15,/vertical
  ;cgColorbar, Divisions=4, Minor=5, Range=[-0.08,0.08],position=[0.058,0.07,0.065,0.37],title='Albedo', TCHARSIZE=15,/vertical
  loadct, 33
  ; map
  for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 50 and lat1(i,j) lt 61 and lon1(i,j) gt 272 and lon1(i,j) lt 276.8 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+25)/50.*254., psym=8, symsize=1

  ; regionⅠ
  ;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 58 and lat1(i,j) lt 60 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+15)/30.*254., psym=8, symsize=1

  ; regionⅣ
  ;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 52 and lat1(i,j) lt 54 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+15)/30.*254., psym=8, symsize=1

  ; regionⅡ
  ;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 56 and lat1(i,j) lt 58 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+15)/30.*254., psym=8, symsize=1

  ; regionⅢ
  ;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 54 and lat1(i,j) lt 56 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+15)/30.*254., psym=8, symsize=1
  ;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 54 and lat1(i,j) lt 56 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+0.08)/0.16*254., psym=8, symsize=1

  ; ヒストグラムのプロット
  loadct, 39
  dif = p1-p2_mod
  good = where(abs(dif) ge 0.001 and albedomap1 ge 0.15)
  hist = histogram(dif(good), min=-100, max=100, binsize=1)
  bin = dindgen(201) - 100d
  ;window, 1, xs=1200, ys=800
  ;!P.Multi = [0, 1, 1]
  plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,2000], xr=[-50,50], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa'
  
  
;  good = where(abs(dif) ge 0.00001 and p1 ge 0.15 and SZA1 le 60)
;  hist = histogram(dif(good), min=-0.05, max=0.05, binsize=0.001)
;  bin = dindgen(101)*0.001 - 0.05d
;  plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,600], xr=[-0.03,0.03], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Albedo',position=[0.6,0.1,0.95,0.4]

  stop
  snapshot = TVRD(True=1)
  ;Write_JPEG, 'Albedo_920-931.jpg', snapshot, True=1, Quality=100

  stop
end