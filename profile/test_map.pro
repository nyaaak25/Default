
; ヒストグラムを作成するプログラム

FUNCTION MYGAUSS, X, P
  RETURN, P[0] + GAUSS1(X, P[1:3])
END


Pro test_map
device, decomposed = 0, retain = 2
loadct, 39

; ヒストグラムを2個重ねて表示させる
minlon = 274
maxlon = 277
minlat = 54
maxlat = 56

; EWのヒストグラムを作成
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0920_3.sav'
TAMAP1 = TAMAP
albedomap1 = inputalbedo
Alt1 = altitude
lat1 = lati
lon1 = longi
p1 = pressure
p1 = exp(p1)

; *** calc index and create array ***
; index search
; 1つ目のORB用
ip1 = n_elements(lat1(*,0))
io1 = n_elements(lat1(0,*))

; *** 差し引いたfileのrestore ***
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/mode1_EW_work_920-931.sav'
p1_mod = p2_mod

; deff p2_mod
deff = dblarr(ip1, io1)
deff_1 = p1 - p1_mod

for i = 0, ip1-1 do for j = 0, io1-1 do if p1(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then deff(i,j) = deff_1(i,j)
for i = 0, ip1-1 do for j = 0, io1-1 do if deff(i,j) eq 0 then deff(i,j) = -0d/0d
; *** histogram ***
good = where(abs(deff) ge 0.01)
med_hist = deff(good)

stop

; fitting
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_B2_update.sav'
TAMAP2 = TAMAP
albedomap2 = inputalbedo
Alt2 = altitude
lat2 = lati
lon2 = longi
p2 = pressure

; *** 差し引いたfileのrestore ***
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/mode1_fit_work_920-931.sav'
p_mod = p2_mod
deff2 = p2 - p_mod

; *** histogram ***
good = where(abs(deff) ge 0.01)
med_hist = deff(good)

good2 = where(abs(deff2) ge 0.01)
med_hist2 = deff2(good2)

loadct, 39
hist = histogram(med_hist, min=-100, max=100, binsize=1)
bin = (findgen(n_elements(med_hist))*1) + min(med_hist)

hist2 = histogram(med_hist2, min=-100, max=100, binsize=1)
bin2 = (findgen(n_elements(med_hist2))*1) + min(med_hist2)


window, 2, xs=1200, ys=800
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,250], xr=[-35,35], thick=3, charsize=2, title='Relative error', xtitle='Pa'
;plots, bin2, hist2, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,250], xr=[-35,35], thick=3











stop


















restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0920_3.sav'
ind = where_xyz(longi ge 274 and longi le 277 and lati ge 54 and lati le 56, xind=xind, yind=yind)
SZA1 = reform(geocube(*,8,*))*1.e-4
EA1 = reform(geocube(*,9,*))*1.e-4
PA1 = reform(geocube(*,10,*))*1.e-4
Alt1 = reform(geocube(*,12,*))

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0920_3.sav'
ind = where_xyz(longi ge 272 and longi le 277 and lati ge 50 and lati le 61, xind=xind, yind=yind)
SZA2 = reform(geocube(*,8,*))*1.e-4
EA2 = reform(geocube(*,9,*))*1.e-4
PA2 = reform(geocube(*,10,*))*1.e-4
Alt2 = reform(geocube(*,12,*))
SZA2_mod = SZA2
EA2_mod = EA2

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0920_3.sav'
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_B2_update.sav'
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_albedo_update.sav'
lon1 = longi
lat1 = lati
p1 = exp(pressure)
;albedomap1 = albedomap

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0931_3.sav'
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0931_3_B2_update.sav'
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0931_3_albedo_update.sav'
lat2 = lati
lon2 = longi
p2 = exp(pressure)
p2_mod = p2
p2_mod(*,*) = -0d/0d

;albedomap2 = albedomap
;albedomap2_mod = albedomap
;albedomap2_mod(*,*) = -0d/0d

for i = 0, 127 do begin
  for j = 0, 595 do begin
    a = WHERE_XYZ(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*)) eq min(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*))), XIND=xind, YIND=yind, count)
    if count ge 1 then begin
      if abs(lat1(i,j)-lat2(xind(0),yind(0))) le 0.03 and abs(lon1(i,j)-lon2(xind(0),yind(0))) le 0.03 then begin
        p2_mod(i,j) = p2(xind(0),yind(0))
        ;albedomap2_mod(i,j) = albedomap2(xind(0),yind(0))
        ;SZA2_mod(i,j) = SZA2(xind(0),yind(0))
        ;EA2_mod(i,j) = EA2(xind(0),yind(0))
      endif
    endif
  endfor
endfor
help, p1, p2
stop

A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

; ヒストグラムにガウシアンをfitさせてみる
; ヒストグラムのbinとかを作る
dif = p1-p2_mod
;good = where(abs(dif) ge 0.001)
good = where(abs(dif) ge 0.001 and albedomap1 ge 0.15 and SZA2 le 60)
hist = histogram(dif(good), min=-100, max=100, binsize=1)
bin = dindgen(201) - 100d

; ヒストグラムにガウシアンをfittingさせる

ind = where(bin ge -70 and bin le 70)
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
start = [11d, 1d, 5d, 2000d]
hist1 = hist + 10d
hist2 = hist1[ind]
bin1 = bin[ind]

pi(0).fixed = 1
pi(2).limited(*) = 1
pi(2).limits(1) = 1d
pi(2).limits(1) = 7d


err = hist2*10
result = MPFITFUN('MYGAUSS', bin1, hist2, err, start, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, NPRINT=0, yfit=yfit, /nan)

window, 3, xs=700, ys=700
!P.Multi = [0, 2, 2]
loadct, 39
plot, bin1, hist2, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,300], xr=[-70,70], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa'
oplot, bin1, result(0)+gauss1(bin1, result(1:3)), color=254, thick=5

stop


; アルベドと地表面気圧の差を見る

loadct,39
dif1 = p1-p2_mod
dif2 = albedomap1-albedomap2_mod

window, 3, xs=800, ys=800
plot, findgen(10), xr=[-0.026,0.026],yr=[-15.5,15.5],back=255, color=0, /nodata, charsize=2, title='Albedo / SP',xs=1, ys=1, xtitle='Albedo', ytitle='Pressure (Pa)'
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and dif2(i,j) gt -0.025 and dif2(i,j) lt 0.025 and dif1(i,j)gt -15 and dif1(i,j) lt 15 then plots, dif2(i,j),dif1(i,j), color=0, psym=8, symsize=0.8
;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and albedomap1(i,j) ge 0.15 and SZA2(i,j) le 60 then oplot, dif2,color = 254

stop


; surface pressure image
window, 1, xs=1200, ys=1400
!P.Multi = [0, 2, 2]

loadct, 39
;plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='SZA', xtitle='Lon', ytitle='Lat'
;for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(SZA1(i,j))/75.*254., psym=8, symsize=1

;plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='EA', xtitle='Lon', ytitle='Lat'
;for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(EA(i,j))/30.*254., psym=8, symsize=1

plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-622.)/212.*254., psym=8, symsize=0.8

;plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb931, from 600 Pa to 1200 Pa'
;for i = 0, 127 do for j = 0, 595 do if p2(i,j) gt 0 then plots, lon2(i,j), lat2(i,j), color=(p2(i,j)-600.)/300.*254., psym=8, symsize=1

loadct, 39
plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb931, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p2_mod(i,j)-622.)/212.*254., psym=8, symsize=0.8

loadct, 39
plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920-Orb931, from -50 Pa to 50 Pa', xtitle='Lon', ytitle='Lat'
loadct, 33
;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+12)/24.*254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and albedomap1(i,j) ge 0.15 and SZA2(i,j) le 60 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+12)/24.*254., psym=8, symsize=0.8

loadct, 39
dif = p1-p2_mod
;good = where(abs(dif) ge 0.001)
good = where(abs(dif) ge 0.001 and albedomap1 ge 0.15 and SZA2 le 60)
hist = histogram(dif(good), min=-100, max=100, binsize=1)
bin = dindgen(201) - 100d
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,300], xr=[-30,30], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa'

stop

; Albedo image

window, 1, xs=1200, ys=1400
!P.Multi = [0, 2, 2]

plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920, 0.15 to 0.3', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(albedomap1(i,j)-0.15d)/0.15d*254., psym=8, symsize=0.8

loadct, 39
plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb931, from 0.15 to 0.3', xtitle='Lon', ytitle='Lat'
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(albedomap2_mod(i,j)-0.15d)/0.15d*254., psym=8, symsize=0.8

loadct, 39
plot, findgen(10),  yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920-Orb931', xtitle='Lon', ytitle='Lat'
loadct, 33
;for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-p2_mod(i,j)+12)/24.*254., psym=8, symsize=1
for i = 0, 127 do for j = 0, 595 do if p2_mod(i,j) gt 0 and albedomap1(i,j) ge 0.15 and SZA2(i,j) le 60 then plots, lon1(i,j), lat1(i,j), color=(albedomap1(i,j)-albedomap2_mod(i,j)+0.025d)/0.05d*254., psym=8, symsize=0.8

loadct, 39
dif = albedomap1-albedomap2_mod
;good = where(abs(dif) ge 0.001)
good = where(abs(dif) ge 0.00001 and albedomap1 ge 0.15 and SZA2 le 60)
hist = histogram(dif(good), min=-0.05, max=0.05, binsize=0.001)
bin = dindgen(101)*0.001 - 0.05d
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,600], xr=[-0.05,0.05], thick=3, charsize=2, title='Histogram: Orb920-Orb931'

stop



window, 0, xs=1200, ys=1200
!P.Multi = [0, 2, 2]

plot, findgen(10), yr=[55.5, 58.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2, title='Orb920, from 500 Pa to 1100 Pa', xtitle='Lon', ytitle='Lat'
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