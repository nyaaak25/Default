; QL作成のための練習
; EW_work fileを読み込んで、QLをより詳しくplotしたいときのためのプログラム

pro data_look

device, decomposed = 0, retain = 2
loadct, 39

; 外れ値を探す
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0931_3.sav'

lat1 = lati
minlati = min(lati)
maxlati = max(lati)
lon1 = longi
minlon = min(longi)
maxlon = max(longi)

ip = n_elements(LATI(*,0))
io = n_elements(LATI(0,*))

; retrieval pressure
p1 = pressure
p1 = exp(p1)
for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 1 then p1(i,j) = !VALUES.F_NAN
for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 870 then p1(i,j) = !VALUES.F_NAN
minp1 = min(p1)
maxp1 = max(p1)
medip1 = median(p1)

; MCD pressure
p2 = MCDpressure
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 0d then p2(i,j) = !VALUES.F_NAN
minp2 = min(p2)
maxp2 = max(p2)

; derivetion pressure
p3 = p1 - p2
medip3 = median(p3)

; 何%変動までみたいかによってcolor rangeを変更
; 現在は5%変動までをみられるように設定
dev_range = 0.05*medip1
minp3 = medip3 - dev_range
maxp3 = medip3 + dev_range

; dust opacity
p4 = dustmap
for i = 0, ip-1 do for j = 0, io-1 do if p4(i,j) eq 0d then p4(i,j) = !VALUES.F_NAN
minp4 = min(p4)
maxp4 = max(p4)

; MOLA altitude
p5 = altitude
minp5 = min(p5)
maxp5 = max(p5)

; albedo input
p6 = inputalbedo
for i = 0, ip-1 do for j = 0, io-1 do if p6(i,j) eq 0d then p6(i,j) = !VALUES.F_NAN
minp6 = min(p6)
maxp6 = max(p6)

; プロットの枠組み
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

window, 0, xs=1800, ys=1600
!P.Multi = [0, 3, 2]

; Retrieved pressure plot
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.55,0.31,0.92], /nodata, charsize=3, title='Retrieved pressure', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.04,0.58,0.046,0.86],title='Pa', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1

; MCD prediction pressure
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.55,0.65,0.92], /nodata, charsize=3, title='MCD pressure', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp2,maxp2],position=[0.38,0.58,0.386,0.86],title='Pa', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p2(i,j)-minp2)/(maxp2-minp2)*254., psym=8, symsize=1

; retrieved - MCD
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.55,0.99,0.92], /nodata, charsize=3, title='Retrieval - MCD', xtitle='Lon', ytitle='Lat'
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp3,maxp3],position=[0.72,0.58,0.726,0.86],title='Pa', TCHARSIZE=10, /vertical
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p3(i,j)-minp3)/(maxp3-minp3)*254., psym=8, symsize=1

; input Albedo map
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.08,0.31,0.45], /nodata, charsize=3, title='Albedo map', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp6,maxp6],position=[0.04,0.12,0.046,0.4],title='Albedo', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p6(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p6(i,j)-minp6)/(maxp6-minp6)*254., psym=8, symsize=1

; dust opacity map
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.08,0.65,0.45], /nodata, charsize=3, title='dust opacity', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp4,maxp4],position=[0.38,0.12,0.386,0.4],title='Dust opacity', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p4(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p4(i,j)-minp4)/(maxp4-minp4)*254., psym=8, symsize=1

; MOLA altitude map
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.45], /nodata, charsize=3, title='MOLA altitude', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp5,maxp5],position=[0.72,0.12,0.726,0.4],title='m', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p5(i,j)-minp5)/(maxp5-minp5)*254., psym=8, symsize=1

loadct, 39
xyouts,0.06,0.975,'ORB0931_3',charsize=2.5,color=0,/normal
xyouts,0.2,0.975,'Time : 2004/10/10 10:14',color=0,/normal,charsize=2
xyouts,0.45,0.975,'Ls : 98.936',color=0,/normal,charsize=2
xyouts,0.57,0.975,'Local Time : 16.0292',color=0,/normal,charsize=2
stop
snapshot = TVRD(True=1)
Write_JPEG, 'SP_ORB0931_3.jpg', snapshot, True=1, Quality=75


stop
; パラパラ漫画を作るための練習
minp1 = 800
maxp1 = 1400

window, 0, xs=4000, ys=800
!P.Multi = [0, 2, 2]

; 左上
loadct,39
plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360],  position=[0.1,0.6,0.49,0.9], back=255, color=0, /nodata, charsize=2, title='Retrieved pressure', xtitle='Longitude [deg]', ytitle='Latitude [deg]'
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=1, Range=[minp1,maxp1],position=[0.04,0.62,0.044,0.86], title='Pa', TCHARSIZE=1, /vertical

;左下
plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360],  position=[0.1,0.15,0.49,0.45], back=255, color=0, /nodata, charsize=2, title='Retrieved pressure', xtitle='Longitude [deg]', ytitle='Latitude [deg]'
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=1, Range=[minp1,maxp1],position=[0.04,0.17,0.044,0.41], title='Pa', TCHARSIZE=1, /vertical

;右上 
plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360],  position=[0.61,0.6,0.99,0.9], back=255, color=0, /nodata, charsize=2, title='Retrieved pressure', xtitle='Longitude [deg]', ytitle='Latitude [deg]'
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=1, Range=[minp1,maxp1],position=[0.55,0.62,0.554,0.86], title='Pa', TCHARSIZE=1, /vertical

;右下
plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360],  position=[0.61,0.15,0.99,0.45], back=255, color=0, /nodata, charsize=2, title='Retrieved pressure', xtitle='Longitude [deg]', ytitle='Latitude [deg]'
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=1, Range=[minp1,maxp1],position=[0.55,0.17,0.554,0.41], title='Pa', TCHARSIZE=1, /vertical



stop

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0351_3.sav'
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/EW_work_ORB0920_3.sav'
ip = n_elements(LATI(*,0))
io = n_elements(LATI(0,*))

lat1 = lati
minlati = min(lati)
;minlati = 18
maxlati = max(lati)
;maxlati = 27
lon1 = longi
minlon = min(longi)
;minlon = 201
maxlon = max(longi)
;maxlon = 204

lati_ind = where_xyz(lat1 ge minlati and lat1 le maxlati and lon1 ge minlon and lon1 le maxlon,xind=xind, yind=yind)


; retrieval pressure
p1 = pressure
pp1 = dblarr(ip,io)
p1 = exp(p1)
for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 1 then p1(i,j) = !VALUES.F_NAN

pp1(lati_ind) = p1(lati_ind)
for i = 0, ip-1 do for j = 0, io-1 do if pp1(i,j) eq 0 then pp1(i,j) = !VALUES.F_NAN
minp1 = min(pp1)
maxp1 = max(pp1)
medip1 = median(pp1)
dev_range = 0.01*medip1
minp1 = minp1 - dev_range
maxp1 = maxp1 + dev_range

; MCD pressure
p2 = MCDpressure
pp2 = dblarr(ip,io)
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 0 then p2(i,j) = !VALUES.F_NAN

pp2(lati_ind) = p2(lati_ind)
for i = 0, ip-1 do for j = 0, io-1 do if pp2(i,j) eq 0 then pp2(i,j) = !VALUES.F_NAN
minp2 = min(pp2)
maxp2 = max(pp2)

; derivetion pressure
p3 = pp1 - pp2
medip3 = median(p3)

; 何%変動までみたいかによってcolor rangeを変更
dev_range = 0.05*medip1
minp3 = medip3 - dev_range
maxp3 = medip3 + dev_range

; dust opacity
p4 = dustmap
for i = 0, ip-1 do for j = 0, io-1 do if p4(i,j) eq 0 then p4(i,j) = !VALUES.F_NAN
minp4 = min(p4)
maxp4 = max(p4)

; MOLA altitude
p5 = altitude
pp5 = dblarr(ip,io)

pp5(lati_ind) = p5(lati_ind)
for i = 0, ip-1 do for j = 0, io-1 do if pp5(i,j) eq 0 then pp5(i,j) = !VALUES.F_NAN

minp5 = min(pp5)
maxp5 = max(pp5)

;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0920_3.sav'

; albedo input
p6 = inputalbedo

; altitude
;Alt1 = reform(alt)
;p6 = altitude
minp6 = min(p6)
maxp6 = max(p6)

; derivetion pressure
p7 = pp1 - pp5
medip7 = median(p7)

; 何%変動までみたいかによってcolor rangeを変更
dev_range1 = 0.05*medip1
minp7 = medip7 - dev_range1
maxp7 = medip7 + dev_range1


; プロットの枠組み
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

window, 0, xs=1800, ys=1600
!P.Multi = [0, 3, 2]

; Retrieved pressure plot
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.55,0.31,0.92], /nodata, charsize=3, title='Retrieved pressure', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.04,0.58,0.046,0.86],title='Pa', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 0 and lat1(i,j) gt minlati and lat1(i,j) lt maxlati and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(pp1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1

; MCD prediction pressure
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.55,0.65,0.92], /nodata, charsize=3, title='MCD pressure', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp2,maxp2],position=[0.38,0.58,0.386,0.86],title='Pa', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 and lat1(i,j) gt minlati and lat1(i,j) lt maxlati and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(pp2(i,j)-minp2)/(maxp2-minp2)*254., psym=8, symsize=1

; retrieved - MCD
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.55,0.99,0.92], /nodata, charsize=3, title='Retrieval - MCD', xtitle='Lon', ytitle='Lat'
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp3,maxp3],position=[0.72,0.58,0.726,0.86],title='Pa', TCHARSIZE=10, /vertical
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 and lat1(i,j) gt minlati and lat1(i,j) lt maxlati and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p3(i,j)-minp3)/(maxp3-minp3)*254., psym=8, symsize=1

; input Albedo map
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.08,0.31,0.45], /nodata, charsize=3, title='Albedo map', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp6,maxp6],position=[0.04,0.12,0.046,0.4],title='Albedo', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p6(i,j) gt 0 and lat1(i,j) gt minlati and lat1(i,j) lt maxlati and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p6(i,j)-minp6)/(maxp6-minp6)*254., psym=8, symsize=1

; dust opacity map
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.08,0.65,0.45], /nodata, charsize=3, title='dust opacity', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp4,maxp4],position=[0.38,0.12,0.386,0.4],title='Dust opacity', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p4(i,j) gt 0 and lat1(i,j) gt minlati and lat1(i,j) lt maxlati and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p4(i,j)-minp4)/(maxp4-minp4)*254., psym=8, symsize=1

; MOLA altitude map
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.45], /nodata, charsize=3, title='MOLA altitude', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp5,maxp5],position=[0.72,0.12,0.726,0.4],title='m', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if lat1(i,j) gt minlati and lat1(i,j) lt maxlati and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(pp5(i,j)-minp5)/(maxp5-minp5)*254., psym=8, symsize=1

; 疑似地形除去
;loadct,39
;plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.45], /nodata, charsize=3, title='MOLA altitude', xtitle='Lon', ytitle='Lat'
;cgLOADCT, 39
;cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp7,maxp7],position=[0.72,0.12,0.726,0.4],title='m', TCHARSIZE=10, /vertical
;loadct, 39
;for i = 0, ip-1 do for j = 0, io-1 do if lat1(i,j) gt minlati and lat1(i,j) lt maxlati and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then plots, lon1(i,j), lat1(i,j), color=(p7(i,j)-minp7)/(maxp7-minp7)*254., psym=8, symsize=1


loadct, 39
xyouts,0.06,0.975,'ORB0351_3',charsize=2.5,color=0,/normal


snapshot = TVRD(True=1)
Write_JPEG, 'ORB0351_5per.jpg', snapshot, True=1, Quality=75

stop
; プロットの枠組み
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

window, 0, xs=2000, ys=1600
!P.Multi = [0, 3, 3]

; 上
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.68,0.31,0.92], /nodata, charsize=3, title='Retrieved pressure', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.04,0.58,0.046,0.86],title='Pa', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1

loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.55,0.65,0.92], /nodata, charsize=3, title='MCD pressure', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp2,maxp2],position=[0.38,0.58,0.386,0.86],title='Pa', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p2(i,j)-minp2)/(maxp2-minp2)*254., psym=8, symsize=1

loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.55,0.99,0.92], /nodata, charsize=3, title='Retrieval - MCD', xtitle='Lon', ytitle='Lat'
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp3,maxp3],position=[0.72,0.58,0.726,0.86],title='Pa', TCHARSIZE=10, /vertical 
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p3(i,j)-minp3)/(maxp3-minp3)*254., psym=8, symsize=1

; 中
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.41,0.31,0.65], /nodata, charsize=3, title='Albedo map', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp6,maxp6],position=[0.04,0.12,0.046,0.4],title='Pa', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p6(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p6(i,j)-minp6)/(maxp6-minp6)*254., psym=8, symsize=1

loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.08,0.65,0.45], /nodata, charsize=3, title='dust opacity', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp6,maxp6],position=[0.38,0.12,0.386,0.4],title='Dust opacity', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p4(i,j)-minp4)/(maxp4-minp4)*254., psym=8, symsize=1

loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.45], /nodata, charsize=3, title='MOLA altitude', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp5,maxp5],position=[0.72,0.12,0.726,0.4],title='m', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p5(i,j)-minp5)/(maxp5-minp5)*254., psym=8, symsize=1

; 下
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.14,0.31,0.38], /nodata, charsize=3, title='MOLA altitude', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp5,maxp5],position=[0.72,0.12,0.726,0.4],title='m', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p5(i,j)-minp5)/(maxp5-minp5)*254., psym=8, symsize=1

loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.45], /nodata, charsize=3, title='MOLA altitude', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp5,maxp5],position=[0.72,0.12,0.726,0.4],title='m', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p5(i,j)-minp5)/(maxp5-minp5)*254., psym=8, symsize=1

loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.45], /nodata, charsize=3, title='MOLA altitude', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp5,maxp5],position=[0.72,0.12,0.726,0.4],title='m', TCHARSIZE=10, /vertical
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p5(i,j)-minp5)/(maxp5-minp5)*254., psym=8, symsize=1

loadct, 39
xyouts,0.06,0.975,'ORB0920',charsize=2.5,color=0,/normal
xyouts,0.2,0.975,'Time : 2002/3/4 12:41',color=0,/normal,charsize=2
xyouts,0.45,0.975,'Ls : 67', color=0,/normal,charsize=2
xyouts,0.57,0.975,'Local Time : 13.3',color=0,/normal,charsize=2

stop

end