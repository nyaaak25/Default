; 2軌道の差をとるプログラム

Pro topography_test

;おまじない
device, decomposed = 0, retain = 2
loadct,39

;restore file
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0920_3.sav'
ind = where_xyz(longi ge 272 and longi le 277 and lati ge 50 and lati le 61, xind=xind, yind=yind)
SZA1 = reform(geocube(*,8,*))*1.e-4
EA1 = reform(geocube(*,9,*))*1.e-4
PA1 = reform(geocube(*,10,*))*1.e-4
Alt1 = reform(alt)

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0931_3.sav'
SZA2 = reform(geocube(*,8,*))*1.e-4
EA2 = reform(geocube(*,9,*))*1.e-4
PA2 = reform(geocube(*,10,*))*1.e-4
Alt2 = reform(alt)
SZA2_mod = SZA2
EA2_mod = EA2

; 2軌道のファイル
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_B2_update.sav'
lat1 = lati
lon1 = longi
p1 = pressure
albedomap1 = albedomap
TAMAP1 = TAMAP

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0931_3_B2_update.sav'
lat2 = lati
lon2 = longi
p2 = pressure
p2_mod = p2
p2_mod(*,*) = -0d/0d
TAMAP2 = TAMAP

  ; 1番近いピクセルを探してくる
for i = 0, 127 do begin
  for j = 0, 595 do begin
    a = WHERE_XYZ(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*)) eq min(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*))), XIND=xind, YIND=yind, count)
    if count ge 1 then begin
      if abs(lat1(i,j)-lat2(xind(0),yind(0))) le 0.03 and abs(lon1(i,j)-lon2(xind(0),yind(0))) le 0.03 then begin
;      if abs(lat1(i,j)-lat2(xind(0),yind(0))) le 0.02 and abs(lon1(i,j)-lon2(xind(0),yind(0))) le 0.02 then begin
        p2_mod(i,j) = p2(xind(0),yind(0))
        SZA2_mod(i,j) = SZA2(xind(0),yind(0))
        EA2_mod(i,j) = EA2(xind(0),yind(0))
      endif
    endif
  endfor
endfor

g = 3.72 ;m/s2
Gconst = 6.67430d-11
MMars = 6.42d23
RMars = 3.4d6
g = -Gconst*MMars/(-RMars*RMars)
goodT = where(TAMAP1 gt 0)
R = 192d ;J/K/kg

;alt; km
p1_slev = p1
;for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then p1_slev(i,j) = p1(i,j) * exp( (alt1(i,j)*1d3) / (R*TAMAP1(i,j)/g) )
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then p1_slev(i,j) = p1(i,j) * exp( (alt1(i,j)*1d3) / (R*TAMAP1(i,j)/(-Gconst*MMars/(-1d*(RMars+alt(i,J)*1d3)*(RMars+alt(i,J)*1d3) ))))

dif = p1-p2_mod
good = where(abs(dif) ge 0.001 and albedomap1 ge 0.15 and SZA1 le 60 and SZA2_mod le 60 and EA1 le 15 and EA2_mod le 15)

; プロットの枠組み
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

; plot
window, 0, xs=1200, ys=800
!P.Multi = [0, 2, 2]

; SP1のプロット
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[53.8, 56.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='ORB920_3 SPmap', xtitle='Lon', ytitle='Lat',position=[0.14,0.55,0.49,0.9]
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, Range=[834,622],position=[0.058,0.57,0.065,0.87],title='Pressure (Pa)', TCHARSIZE=15,/vertical 
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 and lat1(i,j) gt 54 and lat1(i,j) lt 56 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-622.)/212.*254., psym=8, symsize=1

; topography remove
loadct, 39
plot, findgen(10),xs=1, ys=1, yr=[53.8, 56.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=2, title='topography remove', xtitle='Lon', ytitle='Lat',position=[0.64,0.55,0.99,0.9]
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, Range=[470,682],position=[0.556,0.57,0.563,0.87],title='Pressure (Pa)', TCHARSIZE=15,/vertical
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and lat1(i,j) gt 54 and lat1(i,j) lt 56 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1_slev(i,j)-470.)/212.*254., psym=8, symsize=1

; topography remove +-15
loadct, 39
plot, findgen(10),  xs=1, ys=1,yr=[53.8, 56.2], xr=[273.8, 277.2], back=255, color=0, /nodata, charsize=1.5, title='topography remove 576 +- 15 Pa', xtitle='Lon', ytitle='Lat',position=[0.14,0.1,0.49,0.4]
cgLOADCT, 33
cgColorbar, Divisions=4, Minor=5, Range=[-15,15],position=[0.058,0.07,0.065,0.37],title='Pressure (Pa)', TCHARSIZE=15,/vertical 
loadct, 33
for i = 0, 127 do for j = 0, 595 do if p1_slev(i,j) gt 0 and albedomap1(i,j) ge 0.15 and SZA1(i,j) le 60 and SZA2_mod(i,j) le 60 and EA1(i,j) le 15 and EA2_mod(i,j) le 15 and lat1(i,j) gt 54 and lat1(i,j) lt 56 and lon1(i,j) gt 274 and lon1(i,j) lt 277 then plots, lon1(i,j), lat1(i,j), color=(p1_slev(i,j)-median(p1_slev(good))-15.)/30.*254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p1_slev(i,j) gt median(p1_slev(good))+12. and albedomap1(i,j) ge 0.15 and SZA1(i,j) le 60 and SZA2_mod(i,j) le 60 and EA1(i,j) le 15 and EA2_mod(i,j) le 15then plots, lon1(i,j), lat1(i,j), color=254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p1_slev(i,j) lt median(p1_slev(good))-12. and p1_slev(i,j) gt 0. and albedomap1(i,j) ge 0.15 and SZA1(i,j) le 60 and SZA2_mod(i,j) le 60 and EA1(i,j) le 15 and EA2_mod(i,j) le 15then plots, lon1(i,j), lat1(i,j), color=0., psym=8, symsize=1

loadct, 39
;good = where(abs(dif) ge 0.001)
dif2 = p1_slev - median(p1_slev(good))
dif3 = p1 - median(p1(good))
hist = histogram(dif2(good), min=-1000, max=1000, binsize=1)
hist3 = histogram(dif3(good), min=-1000, max=1000, binsize=1)
bin = dindgen(2001) - 1000d
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,350], xr=[-30,30], thick=3, charsize=1.5,  title='Histogram: Orb920, topography removed, ref: 557 Pa', xtitle='Pa',position=[0.64,0.1,0.99,0.4]


snapshot = TVRD(True=1)
Write_JPEG, 'ORB0920-0930_remove.jpg', snapshot, True=1, Quality=100

  stop
end