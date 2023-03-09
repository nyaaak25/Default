; 高度補正をしてみるプログラム

Pro altitude_correction

;おまじない
device, decomposed = 0, retain = 2
loadct,39

;restore file
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0920_3.sav'
ind = where_xyz(longi ge 274 and longi le 277 and lati ge 54 and lati le 56, xind=xind, yind=yind)
SZA1 = reform(geocube(*,8,*))*1.e-4
EA1 = reform(geocube(*,9,*))*1.e-4
PA1 = reform(geocube(*,10,*))*1.e-4
Alt1 = reform(alt)

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0920_3.sav'
ind = where_xyz(longi ge 272 and longi le 277 and lati ge 50 and lati le 61, xind=xind, yind=yind)
SZA2 = reform(geocube(*,8,*))*1.e-4
EA2 = reform(geocube(*,9,*))*1.e-4
PA2 = reform(geocube(*,10,*))*1.e-4
Alt2 = reform(alt)
SZA2_mod = SZA2
EA2_mod = EA2

;高度補正するfile
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/work_CO2_ORB0920_3.sav'
lat1 = lati
lon1 = longi
p1 = pressure
;albedomap1 = albedomap

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/work_CO2_ORB0931_3.sav'
lat2 = lati
lon2 = longi
p2 = pressure
p2_mod = p2
p2_mod(*,*) = -0d/0d

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_B2_update.sav'
TAMAP2 = TAMAP

; 1番近いピクセルを探してくる
for i = 0, 127 do begin
  for j = 0, 595 do begin
    a = WHERE_XYZ(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*)) eq min(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*))), XIND=xind, YIND=yind, count)
    if count ge 1 then begin
      if abs(lat1(i,j)-lat2(xind(0),yind(0))) le 0.03 and abs(lon1(i,j)-lon2(xind(0),yind(0))) le 0.03 then begin
        p2_mod(i,j) = p2(xind(0),yind(0))
        SZA2_mod(i,j) = SZA2(xind(0),yind(0))
        EA2_mod(i,j) = EA2(xind(0),yind(0))
      endif
    endif
  endfor
endfor

; プロットの枠組み
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

; plot
window, 0, xs=2400, ys=1200
!P.Multi = [0, 4, 2]

; SZAのplot
loadct, 39
plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='SZA', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(SZA1(i,j))/75.*254., psym=8, symsize=1

; EAのプロット
plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='EA', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(EA1(i,j))/30.*254., psym=8, symsize=1

; SP1のプロット
plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-622.)/212.*254., psym=8, symsize=1

; SP2のプロット
plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb931, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p2_mod(i,j)-622.)/212.*254., psym=8, symsize=1

; ORB920-0RB931を差し引いた図
loadct, 39
plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920-Orb931, from -50 Pa to 50 Pa', xtitle='Lon', ytitle='Lat'
loadct, 33
;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+12)/24.*254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and albedomap1(i,j) ge 0.15 and SZA1(i,j) le 60 and SZA2_mod(i,j) le 60 and EA1(i,j) le 15 and EA2_mod(i,j) le 15 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+12)/24.*254., psym=8, symsize=1

; ヒストグラムのプロット
loadct, 39
dif = p1-p2_mod
;good = where(abs(dif) ge 0.001)
good = where(abs(dif) ge 0.001 and albedomap1 ge 0.15 and SZA1 le 60 and SZA2_mod le 60 and EA1 le 15 and EA2_mod le 15)
hist = histogram(dif(good), min=-100, max=100, binsize=1)
bin = dindgen(201) - 100d
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,300], xr=[-30,30], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa'

; MOLA高度のマップ
loadct, 39
plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='SZA', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(alt1(i,j)-6d)/5d*254., psym=8, symsize=1

;correction of Topograohical effects --->
g = 3.72d ;m/s
Gconst = 6.67430d-11
MMars = 6.42d23
RMars = 3.4d6
g = -Gconst*MMars/(-RMars*RMars)
goodT = where(TAMAP1 gt 0)
R = 192d ;mks unit
zref = mean(alt1(ind))
p1_slev = p1

; 重力補正なし
;for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then p1_slev(i,j) = p1(i,j) * exp( - ((zref-alt1(i,j)*1d3) / (R*TAMAP1(i,j)/g)) )
; 重力補正あり
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then p1_slev(i,j) = p1(i,j) * exp( (alt1(i,j)*1d3) / (R*TAMAP1(i,j)/(-Gconst*MMars/(-1d*(RMars+alt(i,J)*1d3)*(RMars+alt(i,J)*1d3) ))))

; color rangeが557 Paから±106 Pa
loadct, 39
plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=3, title='Orb920, topography removed, 557 Pa +/- 106 Pa', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1_slev(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1_slev(i,j)-median(p1_slev(good))-106.)/212.*254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p1_slev(i,j) gt median(p1_slev(good))+106 then plots, lon1(i,j), lat1(i,j), color=254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p1_slev(i,j) lt median(p1_slev(good))-106. and p1_slev(i,j) gt 0. then plots, lon1(i,j), lat1(i,j), color=0., psym=8, symsize=1

; color rangeが557 Paから±12 Pa
loadct, 39
plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=3, title='Orb920, topography removed, 557 Pa +/- 12 Pa', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p1_slev(i,j) gt 0 and albedomap1(i,j) ge 0.15 and SZA1(i,j) le 60 and SZA2_mod(i,j) le 60 and EA1(i,j) le 15 and EA2_mod(i,j) le 15then plots, lon1(i,j), lat1(i,j), color=(p1_slev(i,j)-median(p1_slev(good))-12.)/24.*254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p1_slev(i,j) gt median(p1_slev(good))+12. and albedomap1(i,j) ge 0.15 and SZA1(i,j) le 60 and SZA2_mod(i,j) le 60 and EA1(i,j) le 15 and EA2_mod(i,j) le 15then plots, lon1(i,j), lat1(i,j), color=254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p1_slev(i,j) lt median(p1_slev(good))-12. and p1_slev(i,j) gt 0. and albedomap1(i,j) ge 0.15 and SZA1(i,j) le 60 and SZA2_mod(i,j) le 60 and EA1(i,j) le 15 and EA2_mod(i,j) le 15then plots, lon1(i,j), lat1(i,j), color=0., psym=8, symsize=1

; altitude corecttionをした後のヒストグラムの表示
loadct, 39
;good = where(abs(dif) ge 0.001)
dif2 = p1_slev - median(p1_slev(good))
dif3 = p1 - median(p1(good))
hist = histogram(dif2(good), min=-1000, max=1000, binsize=1)
hist3 = histogram(dif3(good), min=-1000, max=1000, binsize=1)
bin = dindgen(2001) - 1000d
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,400], xr=[-30,30], thick=3, charsize=3, title='Histogram: Orb920, topography removed, ref: 557 Pa', xtitle='Pa'

;plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,500], xr=[-230,230], thick=3, charsize=3, title='Histogram: Orb920, topography removed, ref: 557 Pa', xtitle='Pa'
;oplot, bin, hist3, psym=10, color=254

end