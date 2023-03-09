pro QL_create

device, decomposed = 0, retain = 2
loadct,39

restore, '/work1/LUT/SP/table/absorption/work_CO2_ORB0920_3_b2.sav'
lat1 = lati
minlati = min(lati)
maxlati = max(lati)
lon1 = longi
minlon = min(longi)
maxlon = max(longi)
p1 = pressure
p1 = exp(p1)

for i = 0, 127 do for j = 0, 595 do if p1(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN
minp1 = min(p1)
maxp1 = max(p1)


; プロットの枠組み
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

; plot
window, 0, xs=1200, ys=800
!P.Multi = [0, 3, 1]


loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.08,0.31,0.92], /nodata, charsize=3, title='Retrieved pressure', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.04,0.12,0.046,0.9],title='Pa', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1

loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.08,0.65,0.92], /nodata, charsize=3, title='MCD pressure', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.38,0.12,0.386,0.9],title='Pa', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1

loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.92], /nodata, charsize=3, title='MCD pressure', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.72,0.12,0.726,0.9],title='Pa', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, 127 do for j = 0, 595 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1

loadct, 39
xyouts,0.06,0.975,'ORB0920_3',charsize=2.5,color=0,/normal
xyouts,0.2,0.975,'Time : 2004/12/13 12:11',color=0,/normal,charsize=2
xyouts,0.45,0.975,'Ls : 54.111',color=0,/normal,charsize=2
xyouts,0.57,0.975,'Local Time : 12.222',color=0,/normal,charsize=2


stop
end
