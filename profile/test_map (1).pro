; 高度補正をするプログラム

Pro altitude_correction

device, decomposed = 0, retain = 2
loadct, 39

restore, '/Volumes/2014_2015/OMEGA/ORB0920_3.sav'
ind = where_xyz(longi ge 274 and longi le 277 and lati ge 54 and lati le 56, xind=xind, yind=yind)
SZA1 = reform(geocube(*,8,*))*1.e-4
EA1 = reform(geocube(*,9,*))*1.e-4
PA1 = reform(geocube(*,10,*))*1.e-4
Alt1 = reform(alt)

restore, '//Volumes/2014_2015/OMEGA/ORB0931_3.sav'
ind = where_xyz(longi ge 272 and longi le 277 and lati ge 50 and lati le 61, xind=xind, yind=yind)
SZA2 = reform(geocube(*,8,*))*1.e-4
EA2 = reform(geocube(*,9,*))*1.e-4
PA2 = reform(geocube(*,10,*))*1.e-4
Alt2 = reform(alt)
SZA2_mod = SZA2
EA2_mod = EA2

;restore, '/Users/shohei/tmp/pressuremap_ORB0920_3.sav'
restore, '/Users/shohei/tmp/SPmap_ORB0920_3_albedo_update.sav'
lat1 = lati
lon1 = longi
p1 = pressure
albedomap1 = albedomap
TAMAP1 = TAMAP

;restore, '/Users/shohei/tmp/pressuremap_ORB0931_1.sav'
restore, '/Users/shohei/tmp/SPmap_ORB0931_3_albedo_update.sav'
lat2 = lati
lon2 = longi
p2 = pressure
p2_mod = p2
p2_mod(*,*) = -0d/0d
TAMAP2 = TAMAP

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
help, p1, p2
;stop

A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

window, 0, xs=2400, ys=1200
!P.Multi = [0, 4, 2]

loadct, 39
plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='SZA', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(SZA1(i,j))/75.*254., psym=8, symsize=1

plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='EA', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(EA1(i,j))/30.*254., psym=8, symsize=1

plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-622.)/212.*254., psym=8, symsize=1

;plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb931, from 600 Pa to 1200 Pa'
;for i = 0, 127 do for j = 0, 595 do if p2(i,j) gt 0 then plots, lon2(i,j), lat2(i,j), color=(p2(i,j)-600.)/300.*254., psym=8, symsize=1

loadct, 39
plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb931, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p2_mod(i,j)-622.)/212.*254., psym=8, symsize=1

loadct, 39
plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920-Orb931, from -50 Pa to 50 Pa', xtitle='Lon', ytitle='Lat'
loadct, 33
;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+12)/24.*254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and albedomap1(i,j) ge 0.15 and SZA1(i,j) le 60 and SZA2_mod(i,j) le 60 and EA1(i,j) le 15 and EA2_mod(i,j) le 15 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+12)/24.*254., psym=8, symsize=1

loadct, 39
dif = p1-p2_mod
;good = where(abs(dif) ge 0.001)
good = where(abs(dif) ge 0.001 and albedomap1 ge 0.15 and SZA1 le 60 and SZA2_mod le 60 and EA1 le 15 and EA2_mod le 15)
hist = histogram(dif(good), min=-100, max=100, binsize=1)
bin = dindgen(201) - 100d
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,300], xr=[-30,30], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa'

loadct, 39
plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='SZA', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(alt1(i,j)-6d)/5d*254., psym=8, symsize=1

;correction of Topograohical effects --->
g = 3.72d ;m/s
R = 192d ;mks unit
p1_slev = p1
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then p1_slev(i,j) = p1(i,j) * exp( (-alt1(i,j)*1d3) / (R*TAMAP1(i,j)/g) )

plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1_slev(i,j)-622.)/212.*254., psym=8, symsize=1
stop





window, 0, xs=1200, ys=1200
!P.Multi = [0, 2, 2]

plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-500.)/600.*254., psym=8, symsize=0.8

;plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb931, from 600 Pa to 1200 Pa'
;for i = 0, 127 do for j = 0, 595 do if p2(i,j) gt 0 then plots, lon2(i,j), lat2(i,j), color=(p2(i,j)-600.)/300.*254., psym=8, symsize=1

plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb931, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p2_mod(i,j)-500.)/600.*254., psym=8, symsize=0.8

plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920-Orb931, from -50 Pa to 50 Pa', xtitle='Lon', ytitle='Lat'
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+50)/100.*254., psym=8, symsize=0.8

dif = p1-p2_mod
;good = where(abs(dif) ge 0.001)
good = where(abs(dif) ge 0.001 and albedomap1 ge 0.18)
hist = histogram(dif(good), min=-100, max=100, binsize=1)
bin = dindgen(201) - 100d
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,100], xr=[-100,100], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa'

stop
end