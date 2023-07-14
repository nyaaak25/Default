; ================================================================================================
; 2軌道の差を取って、導出精度を検証するプログラム
; create by Akira Kazama
; 
; topography_remove_v2: 2023.04.24 Mon 16:08
; ・これまでバラバラになんとなく作っていたプログラムを統一
; ・lati, longiを指定したら勝手にplotしてくれる仕様を追加
; 
;  topography_remove_nyon: 2023.04.25 Tue 18:18
; ・EW法とfitting法の差分図を作成するためのプログラム
; ================================================================================================



Pro topography_remove_nyon

; plotをするためのおまじない
device, decomposed = 0, retain = 2
loadct,39

; ++++++++++++ program start ++++++++++++
;  
; *** restore EW file and store values ***
; EW file include: ret pressure, trans, lati, longi, MOLA, dust, SZA, EA, Albedo, MCD pressure, TA
; 格納されているfile name: pressure, trans, lati, longi, altitude, dustmap, SZA_all, EA_all, inputalbedo, MCDpressure, TAmap
; note: ret pressureをpressureに変換するためには、expを取る必要があり

; EW method
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0931_3.sav'
TAMAP1 = TAMAP
albedomap1 = inputalbedo
Alt1 = altitude
lat1 = lati
lon1 = longi
p1 = pressure
p1 = exp(p1)

; fitting
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0931_3_B2_update.sav'
TAMAP2 = TAMAP
albedomap2 = inputalbedo
Alt2 = altitude
lat2 = lati
lon2 = longi
p2 = pressure

; *** calc index and create array ***
; index search
; 1つ目のORB用
ip1 = n_elements(lat1(*,0))
io1 = n_elements(lat1(0,*))

; 2つ目のORB用
ip2 = n_elements(lat2(*,0))
io2 = n_elements(lat2(0,*))

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
;save, p2_mod, filename='/Users/nyonn/IDLWorkspace/Default/savfile/mode1_EW_work_920-931.sav'

; ++++++++++++ index指定 ++++++++++++
; *** plotする領域を指定 ***
; M論で使用したregion一覧
; RegionⅠ: 58 to 60, RegionⅡ: 56 to 58, RegionⅢ: 54 to 56(クレーター), RegionⅣ: 52 to 54
minlon = 274
maxlon = 277
minlat = 54
maxlat = 56

ind = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat and lati le maxlat, xind=xind, yind=yind)

; ++++++++++++ altitude correction ++++++++++++
g = 3.72d ;m/s
Gconst = 6.67430d-11
MMars = 6.42d23
RMars = 3.4d6
g = -Gconst*MMars/(-RMars*RMars)
R = 192d ;mks uni

; EW result
zref1 = mean(alt1(ind))
p1_slev = dblarr(ip1, io1)
for i = 0, ip1-1 do for j = 0, io1-1 do if p1(i,j) gt 0 then p1_slev(i,j) = p1(i,j) * exp( (alt1(i,j)) / (R*TAMAP1(i,j)/(-Gconst*MMars/(-1d*(RMars+alt1(i,J))*(RMars+alt1(i,J)) ))))

; fitting result
zref2 = mean(alt2(ind))
p2_slev = dblarr(ip2, io2)
for i = 0, ip2-1 do for j = 0, io2-1 do if p2(i,j) gt 0 then p2_slev(i,j) = p2(i,j) * exp( (alt2(i,j)) / (R*TAMAP2(i,j)/(-Gconst*MMars/(-1d*(RMars+alt2(i,J))*(RMars+alt2(i,J)) ))))

; ++++ data selection ++++
; EW
for i = 0, ip1-1 do for j = 0, io1-1 do if p1(i,j) eq 1 then p1(i,j) = -0d/0d
; fitting
for i = 0, ip2-1 do for j = 0, io2-1 do if p2(i,j) eq 0 then p2(i,j) = -0d/0d

; 高度補正
; EW
for i = 0, ip1-1 do for j = 0, io1-1 do if p1_slev(i,j) eq 0 then p1_slev(i,j) = -0d/0d
for i = 0, ip1-1 do for j = 0, io1-1 do if p1_slev(i,j) gt 10000 then p1_slev(i,j) = -0d/0d

; 高度補正
; fitting
for i = 0, ip2-1 do for j = 0, io2-1 do if p2_slev(i,j) eq 0 then p2_slev(i,j) = -0d/0d
for i = 0, ip2-1 do for j = 0, io2-1 do if p2_slev(i,j) gt 10000 then p2_slev(i,j) = -0d/0d

; *** fittingとEWの差をみる ***
deff = p1 - p2
deff_alt = p1_slev - p2_slev

; ++++++++++++ plot start ++++++++++++
; 
; *** プロットの枠組みを指定 ***
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

; *** plot ***
;  +++ retrieval result and altitude correction result +++
window, 0, xs=1200, ys=800
; 何分割で何個表示させたいかを決める
!P.Multi = [0, 2, 2]

; *** color pallet用 ***
dev = min(p1)
dev2 = max(p1) - min(p1)
dev3 = min(p2)
dev4 = max(p2) - min(p2)
dev5 = min(p1_slev)
dev6 = max(p1_slev) - min(p1_slev)
dev7 = min(p2_slev)
dev8 = max(p2_slev) - min(p2_slev)
  
; *** EW ret pressureをplot *** 
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='EW: retrieval', xtitle='Lon', ytitle='Lat',position=[0.15,0.55,0.5,0.9]
cgLOADCT, 16

; fittingと同じcolor bar使用
cgColorbar, Divisions=4, Minor=5, Range=[max(p2),min(p2)],position=[0.058,0.57,0.065,0.87],title='Pa', TCHARSIZE=15,/vertical
loadct, 16
for i = 0, ip1-1 do for j = 0, io1-1 do if p1(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p1(i,j) - dev3)/dev4 *254., psym=8, symsize=1

; 別のcolor barでplot
;cgColorbar, Divisions=4, Minor=5, Range=[max(p1),min(p1)],position=[0.058,0.57,0.065,0.87],title='Pa', TCHARSIZE=15,/vertical
;loadct, 16
;for i = 0, ip1-1 do for j = 0, io1-1 do if p1(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-median(p1) - dev)/dev2 *254., psym=8, symsize=1


; *** fittingのret pressureをplot *** 
loadct, 39
plot, findgen(10),xs=1, ys=1, yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='fitting: retrieval', xtitle='Lon', ytitle='Lat',position=[0.65,0.55,0.99,0.9]
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, Range=[max(p2),min(p2)],position=[0.558,0.57,0.565,0.87],title='Pa', TCHARSIZE=15,/vertical
loadct, 16
for i = 0, ip2-1 do for j = 0, io2-1 do if p2(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p2(i,j)- dev3)/dev4 *254., psym=8, symsize=1


; *** EW: 高度補正をしたplot ***
loadct, 39
plot, findgen(10),xs=1, ys=1, yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='EW: altitude correction', xtitle='Lon', ytitle='Lat',position=[0.15,0.1,0.5,0.4]
cgLOADCT, 16
; fittingと同じcolor bar使用
cgColorbar, Divisions=4, Minor=5, Range=[max(p2_slev),min(p2_slev)],position=[0.058,0.07,0.065,0.37],title='Pressure (Pa)', TCHARSIZE=15,/vertical
loadct, 16
for i = 0, ip1-1 do for j = 0, io1-1 do if p1_slev(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p1_slev(i,j)-dev7)/dev8 *254., psym=8, symsize=1
; 別のcolor bar使用
;cgColorbar, Divisions=4, Minor=5, Range=[max(p1_slev),min(p1_slev)],position=[0.058,0.07,0.065,0.37],title='Pressure (Pa)', TCHARSIZE=15,/vertical
;loadct, 16
;for i = 0, ip1-1 do for j = 0, io1-1 do if p1_slev(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p1_slev(i,j)-median(p1_slev)-dev5)/dev6 *254., psym=8, symsize=1


; *** fitting: 高度補正をしたplot ***
loadct, 39
plot, findgen(10),xs=1, ys=1, yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='fitting: altitude correction', xtitle='Lon', ytitle='Lat',position=[0.65,0.1,0.99,0.4]
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, Range=[max(p2_slev),min(p2_slev)],position=[0.558,0.07,0.565,0.37],title='Pressure (Pa)', TCHARSIZE=15,/vertical
loadct, 16
for i = 0, ip2-1 do for j = 0, io2-1 do if p2_slev(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p2_slev(i,j)-dev7)/dev8 *254., psym=8, symsize=1


;  +++ comparison of EW to fitting +++
window, 1, xs=1200, ys=800
; 何分割で何個表示させたいかを決める
!P.Multi = [0, 2, 2]

maxran = 20
minran = -20

; *** color pallet ***
dev9 = minran
dev10 = maxran - minran

; *** fitting - EW: retrieval ***
loadct, 39
plot, findgen(10),  xs=1, ys=1,yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='EW-fitting: ret', xtitle='Lon', ytitle='Lat',position=[0.15,0.55,0.5,0.9]
cgLOADCT, 33
;cgColorbar, Divisions=4, Minor=5, Range=[-20,20],position=[0.058,0.57,0.065,0.87],title='Pressure (Pa)', TCHARSIZE=15,/vertical
cgColorbar, Divisions=4, Minor=5, Range=[minran,maxran],position=[0.058,0.57,0.065,0.87],title='Pressure (Pa)', TCHARSIZE=15,/vertical
loadct, 33
;for i = 0, ip1 - 1 do for j = 0, io1 - 1 do if deff(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(deff(i,j)- dev9)/dev10 *254., psym=8, symsize=1
for i = 0, ip1 - 1 do for j = 0, io1 - 1 do if p1_slev(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(deff(i,j)- dev9)/dev10 *254., psym=8, symsize=1

; *** fitting - EW: altitude remove ***
loadct, 39
plot, findgen(10),xs=1, ys=1, yr=[minlat - 0.2, maxlat + 0.2], xr=[minlon - 0.2, maxlon + 0.2], back=255, color=0, /nodata, charsize=2, title='EW-fitting: alt', xtitle='Lon', ytitle='Lat',position=[0.65,0.55,0.99,0.9]
cgLOADCT, 33
cgColorbar, Divisions=4, Minor=5, Range=[minran,maxran],position=[0.558,0.57,0.565,0.87],title='Pa', TCHARSIZE=15,/vertical
loadct, 33
for i = 0, ip1-1 do for j = 0, io1-1 do if p1_slev(i,j) gt 0 and lat1(i,j) gt minlat and lat1(i,j) lt maxlat and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(deff_alt(i,j)-dev9)/dev10 *254., psym=8, symsize=1

; *** histgram EW-fitting: retrieval ***
loadct, 39
good = where(abs(deff) ge 0.001)
hist = histogram(deff(good), binsize=1)
bin = (findgen(n_elements(deff(good)))*1) + min(deff(good))
plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,500], xr=[-30,30], thick=3, charsize=2, title='Histogram: retrieval', xtitle='Pa',position=[0.15,0.1,0.5,0.4]


; *** histgram EW-fitting: altituede remove ***
loadct, 39
good1 = where(abs(deff_alt) ge 0.001 )
hist1 = histogram(deff_alt(good1), binsize=1)
bin1 = (findgen(n_elements(deff_alt(good1)))*1) + min(deff_alt(good1))
plot, bin1, hist1, psym=10, back=255, color=0, xs=1, ys=1, yr=[0,500], xr=[-20,20], thick=3, charsize=2, title='Histogram: altitude remove', xtitle='Pa',position=[0.65,0.1,0.99,0.4]
stop  

snapshot = TVRD(True=1)
;Write_JPEG, 'Albedo_920-931.jpg', snapshot, True=1, Quality=100
stop


end