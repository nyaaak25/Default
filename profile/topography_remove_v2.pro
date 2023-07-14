; ================================================================================================
; 2軌道の差を取って、導出精度を検証するプログラム
; create by Akira Kazama
; 
; topography_remove_v2: 2023.04.24 Mon 16:08
; ・これまでバラバラになんとなく作っていたプログラムを統一
; ・lati, longiを指定したら勝手にplotしてくれる仕様を追加
; ================================================================================================



Pro topography_remove_v2

; plotをするためのおまじない
device, decomposed = 0, retain = 2
loadct,39

; ++++++++++++ program start ++++++++++++
;  
; *** restore EW file and store values ***
; EW file include: ret pressure, trans, lati, longi, MOLA, dust, SZA, EA, Albedo, MCD pressure, TA
; 格納されているfile name: pressure, trans, lati, longi, altitude, dustmap, SZA_all, EA_all, inputalbedo, MCDpressure, TAmap
; note: ret pressureをpressureに変換するためには、expを取る必要があり

; ORB0920_3
ORB1 = 'ORB0920_3' ; (TBD), ここでORB名を入れただけで色々変えなくて住むような仕様に変更
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_B2_update.sav'
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0920_3.sav'
TAMAP1 = TAMAP
albedomap1 = inputalbedo
Alt1 = altitude
lat1 = lati
lon1 = longi
p1 = pressure
p1 = exp(p1)

; ORB0931_3
ORB2 = 'ORB0931_3'
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0931_3_B2_update.sav'
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0931_3.sav'
TAMAP2 = TAMAP
albedomap2 = inputalbedo
Alt2 = altitude
lat2 = lati
lon2 = longi
p2 = pressure
p2 = exp(p2)


; *** calc index and create array ***
; index search
; 1つ目のORB用
ip1 = n_elements(lat1(*,0))
io1 = n_elements(lat1(0,*))
; 2つ目のORB用
ip2 = n_elements(lat2(*,0))
io2 = n_elements(lat2(0,*))

; create array
p2_mod = dblarr(ip2, io2)

; *** 2軌道による圧力の差し引きを行う ***
; note: 1度回したら差し引いた圧力が格納されるsav fileを作成し、それを呼び出すようにする。2回目以降はコメントアウト。
; 
; mode1: 全く同じ場所はないので、片方のORBを基準に1番近いpixelを取ってきて差し引きを行う 
; mode2: 空間同士を補間して、空間分解能を向上させる(TBD)
; mode3: 4点/9点平均を行って、S/Nを向上させる(TBD)

; mode1::
;for i = 0, 127 do begin
;  for j = 0, 595 do begin
;    a = WHERE_XYZ(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*)) eq min(abs(lat1(i,j)-lat2(*,*))+abs(lon1(i,j)-lon2(*,*))), XIND=xind, YIND=yind, count)
;    if count ge 1 then begin
;      if abs(lat1(i,j)-lat2(xind(0),yind(0))) le 0.03 and abs(lon1(i,j)-lon2(xind(0),yind(0))) le 0.03 then begin
;        p2_mod(i,j) = p2(xind(0),yind(0))
;      endif
;    endif
;  endfor
;endfor
;
;save, p2_mod, filename='/Users/nyonn/IDLWorkspace/Default/savfile/mode1_fit_work_920-931.sav'

; ++++++++++++ index指定 ++++++++++++
; *** plotする領域を指定 ***
; M論で使用したregion一覧
; RegionⅠ: 58 to 60, RegionⅡ: 56 to 58, RegionⅢ: 54 to 56(クレーター), RegionⅣ: 52 to 54
minlon = 274
maxlon = 277
minlat = 54
maxlat = 56

ind1 = where_xyz(lon1 ge minlon and lon1 le maxlon and lat1 ge minlat and lat1 le maxlat, xind=xind, yind=yind)
ind2 = where_xyz(lon2 ge minlon and lon2 le maxlon and lat2 ge minlat and lat2 le maxlat, xind=xind, yind=yind)

; *** 差し引いたfileのrestore ***
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/mode1_EW_work_920-931.sav'
p2_mod = p2_mod

deff = p1 - p2_mod

; ++++++++++++ altitude correction ++++++++++++
g = 3.72d ;m/s
Gconst = 6.67430d-11
MMars = 6.42d23
RMars = 3.4d6
g = -Gconst*MMars/(-RMars*RMars)
R = 192d ;mks uni

; ORB0920_3
zref1 = mean(alt1(ind1))
p1_slev = dblarr(ip1, io1)
for i = 0, ip1-1 do for j = 0, io1-1 do if p1(i,j) gt 0 then p1_slev(i,j) = p1(i,j) * exp( (alt1(i,j)) / (R*TAMAP1(i,j)/(-Gconst*MMars/(-1d*(RMars+alt1(i,J))*(RMars+alt1(i,J)) ))))

; ORB0931_3
zref2 = mean(alt2(ind2))
p2_slev = dblarr(ip2, io2)
for i = 0, ip2-1 do for j = 0, io2-1 do if p2(i,j) gt 0 then p2_slev(i,j) = p2(i,j) * exp( (alt2(i,j)*1d3) / (R*TAMAP2(i,j)/(-Gconst*MMars/(-1d*(RMars+alt2(i,J)*1d3)*(RMars+alt2(i,J)*1d3) ))))

; ++++ data selection ++++
; ORB0920
for i = 0, ip1-1 do for j = 0, io1-1 do if p1(i,j) eq 1 then p1(i,j) = -0d/0d
for i = 0, ip2-1 do for j = 0, io2-1 do if p2(i,j) eq 1 then p2(i,j) = -0d/0d

; deff p2_mod
for i = 0, ip2-1 do for j = 0, io2-1 do if p2_mod(i,j) eq 0 then p2_mod(i,j) = -0d/0d

; 高度補正
; EW
for i = 0, ip1-1 do for j = 0, io1-1 do if p1_slev(i,j) eq 0 then p1_slev(i,j) = -0d/0d
for i = 0, ip1-1 do for j = 0, io1-1 do if p1_slev(i,j) gt 10000 then p1_slev(i,j) = -0d/0d

; 高度補正
; fitting
for i = 0, ip2-1 do for j = 0, io2-1 do if p2_slev(i,j) eq 0 then p2_slev(i,j) = -0d/0d
for i = 0, ip2-1 do for j = 0, io2-1 do if p2_slev(i,j) gt 10000 then p2_slev(i,j) = -0d/0d


; ++++++++++++ plot start ++++++++++++
; *** color pallet用 ***
dev = min(p1(ind1))
dev2 = max(p1(ind1)) - min(p1(ind1))
dev3 = min(p2)
dev4 = max(p2) - min(p2)
dev5 = min(p1_slev)
dev6 = max(p1_slev) - min(p1_slev)
dev7 = min(p2_slev)
dev8 = max(p2_slev) - min(p2_slev)

; 差何Paまでみたかのplot
maxran = 20
minran = -20

dev9 = minran
dev10 = maxran - minran


; *** プロットの枠組みを指定 ***
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

; *** plot ***
window, 0, xs=1200, ys=800
; 何分割で何個表示させたいかを決める
!P.Multi = [0, 2, 2]
  
; *** 1つめのORBのret pressureをplot *** 
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='ORB920_3 SPmap', xtitle='Lon', ytitle='Lat',position=[0.15,0.55,0.5,0.9]
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, Range=[max(p1(ind1)),min(p1(ind1))],position=[0.058,0.57,0.065,0.87],title='Pa', TCHARSIZE=15,/vertical
loadct, 16
for i = 0, ip1-1 do for j = 0, io1-1 do if p1(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-dev)/dev2 *254., psym=8, symsize=1

; *** 2つめのORBのret pressureをplot *** 
loadct, 39
plot, findgen(10),xs=1, ys=1, yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='ORB931_3 SPmap', xtitle='Lon', ytitle='Lat',position=[0.6,0.55,0.95,0.9]
loadct, 16
for i = 0, ip2-1 do for j = 0, io2-1 do if p2(i,j) ge 0 and lat2(i,j) gt minlat and lat2(i,j) lt maxlat and lon2(i,j) gt minlon and lon2(i,j) lt maxlon then plots, lon2(i,j), lat2(i,j), color=(p2(i,j)-dev)/dev2 *254., psym=8, symsize=1

; *** 2つのORBを差し引いた図のplot *** 
loadct, 39
plot, findgen(10),  xs=1, ys=1,yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='ORB920 - ORB931', xtitle='Lon', ytitle='Lat',position=[0.15,0.1,0.5,0.4]
cgLOADCT, 33
cgColorbar, Divisions=4, Minor=5, Range=[-20,20],position=[0.058,0.07,0.065,0.37],title='Pressure (Pa)', TCHARSIZE=15,/vertical
loadct, 33
for i = 0, ip1 - 1 do for j = 0, io1 - 1 do if p1(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(deff(i,j)-dev9)/dev10 *254., psym=8, symsize=1

; *** histogram ***
good = where(abs(deff) lt 100)
med_hist = deff(good)
loadct, 39
hist = histogram(med_hist, min=-100, max=100, binsize=1)
bin = (findgen(n_elements(med_hist))*1) + min(med_hist)
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,10000], xr=[-50,50], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa',position=[0.6,0.1,0.95,0.4]

stop

; *** 高度補正をしたplot ***
; 同じ color scale
loadct, 39
plot, findgen(10),xs=1, ys=1, yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='altitude correction', xtitle='Lon', ytitle='Lat',position=[0.6,0.55,0.95,0.9]
loadct, 16
for i = 0, ip1-1 do for j = 0, io1-1 do if p1_slev(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p1_slev(i,j)-dev5)/dev6 *254., psym=8, symsize=1

; *** 高度補正をしたplot ***
; ±24Pa の color scale
med = p1_slev - median(p1_slev)

loadct, 39
plot, findgen(10),xs=1, ys=1, yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='median +- 30', xtitle='Lon', ytitle='Lat',position=[0.15,0.1,0.5,0.4]
cgLOADCT, 33
cgColorbar, Divisions=4, Minor=5, Range=[-30,30],position=[0.058,0.07,0.065,0.37],title='Pressure (Pa)', TCHARSIZE=15,/vertical
loadct, 33
for i = 0, ip1-1 do for j = 0, io1-1 do if p1_slev(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(med(i,j)-30.)/60.*254., psym=8, symsize=1

; *** histogram ***
good = where(abs(deff) lt 100)
med_hist = deff(good)
loadct, 39
hist = histogram(med_hist, min=-30, max=30, binsize=1)
bin = (findgen(n_elements(med_hist))*1) + min(med_hist)
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,250], xr=[-35,35], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa',position=[0.6,0.1,0.95,0.4]

stop




;loadct, 39
;good = where(abs(dif) ge 0.001)
;dif2 = p1_slev - median(p1_slev)
;good = where(abs(dif2) ge 0.001 )
;dif3 = p1 - median(p1(good))
;hist = histogram(dif2(good), min=-1000, max=1000, binsize=1)
;hist3 = histogram(dif3(good), min=-1000, max=1000, binsize=1)
;bin = dindgen(2001) - 1000d
;plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,1500], xr=[-70,70], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa',position=[0.6,0.1,0.95,0.4]

stop
; ヒストグラムのプロット
loadct, 39
dif = p1(ind) -p2_mod(ind)
good = where(abs(dif) ge 0.001 );and albedomap1 ge 0.18)
;hist = histogram(dif(good), min=-100, max=100, binsize=1)
hist = histogram(dif, min=-100, max=100, binsize=1)
bin = dindgen(201) - 100d
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,250], xr=[-40,40], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Pa',position=[0.6,0.1,0.95,0.4]
stop  
  
;  good = where(abs(dif) ge 0.00001 and p1 ge 0.15 and SZA1 le 60)
;  hist = histogram(dif(good), min=-0.05, max=0.05, binsize=0.001)
;  bin = dindgen(101)*0.001 - 0.05d
;  plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,600], xr=[-0.03,0.03], thick=3, charsize=2, title='Histogram: Orb920-Orb931', xtitle='Albedo',position=[0.6,0.1,0.95,0.4]

  stop
  snapshot = TVRD(True=1)
  ;Write_JPEG, 'Albedo_920-931.jpg', snapshot, True=1, Quality=100


end