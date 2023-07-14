;------------------------------------------------------------------------------------------------------------------------
;
; 熱潮汐波の解析をするためのプログラム
;
; create by Akira Kazama
; 
; moving_average　 ::2023.6.13 Tue 12:00:00
; 移動平均を取るためのプログラム
;
;------------------------------------------------------------------------------------------------------------------------

pro altitude_correction_test
device, decomposed = 0, retain = 2


; mcdのデータをみてみる！どんなものが格納されているかを確認
; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'
path2 = '/work1/LUT/SP/mcddata/'

; mcdのデータをrestore: lat0, long,30
restore, path2 + 'lon0-lat30_mcddata.sav'
mcd_ps1 = mcd_ps
mcd_ta1 = mcd_ta
MOLA1 = MOLA  ; 1370.00 m
sea_ps1 = dblarr(49,25) ; 49はLT方向、25はLs方向。それぞれ0.5LT, 15°刻み

; mcdのデータをrestore: lat150, long,30
restore, path2 + 'lon150-lat30_mcddata.sav'
mcd_ps2 = mcd_ps
mcd_ta2 = mcd_ta
MOLA2 = MOLA  ;-2533.75 m
sea_ps2 = dblarr(49,25)

; 緯度方向も追加したデータ
restore, path2 + 'lon150_mcddata.sav'
mcd_ps3 = mcd_ps
mcd_ta3 = mcd_ta
MOLA3 = MOLA

sea_ps3 = dblarr(49,25,19)

; crearte LT arr
LT_arr = findgen(49, start = 0)
LT_arr = LT_arr * 0.5

; crearte LS arr
LS_arr = findgen(25, start = 0)
LS_arr = LS_arr * 15

; crearte Lat arr
Lat_arr = findgen(19, Increment = 10,start = -90)
;ind = where(lat_arr gt -70 and lat_arr lt 70)
ind = where(lat_arr ge -30 and lat_arr le 30)

; ser-leveld conversion
g = 3.72 ;m/s2
Gconst = 6.67430d-11
MMars = 6.42d23
RMars = 3.4d6
g = -Gconst*MMars/(-RMars*RMars)
R = 192d ;J/K/kg

for i = 0, 49-1 do for j = 0, 25-1 do  sea_ps1(i,j) = mcd_ps1(i,j) * exp( (MOLA1) / (R*mcd_ta1(i,j)/(-Gconst*MMars/(-1d*(RMars+MOLA1)*(RMars+MOLA1) ))))
for i = 0, 49-1 do for j = 0, 25-1 do  sea_ps2(i,j) = mcd_ps2(i,j) * exp( (MOLA2) / (R*mcd_ta2(i,j)/(-Gconst*MMars/(-1d*(RMars+MOLA2)*(RMars+MOLA2) ))))
for i = 0, 49-1 do for j = 0, 25-1 do for k = 0, 19-1 do sea_ps3(i,j,k) = mcd_ps3(i,j,k) * exp( (MOLA3(k)) / (R*mcd_ta3(i,j,k)/(-Gconst*MMars/(-1d*(RMars+MOLA3(k))*(RMars+MOLA3(k)) ))))

; local time 12時のものをもってくる
mcd_ps_12_1 = sea_ps1(24,*)
mcd_ps_12_2 = sea_ps2(24,*)
mcd_ps_12_3 = sea_ps3(24,*,*)

pure_ps = dblarr(49,25)
pure_ps2 = dblarr(49,25)
pure_ps3 = dblarr(49,25,19)

; LT12を引き算してみる
;for i=0, 49-1 do pure_ps(i,*) = sea_ps1(i,*) - mcd_ps_12

; 割り算をしてみる
for i=0, 49-1 do pure_ps(i,*) = ((sea_ps1(i,*) - mcd_ps_12_1) / mcd_ps_12_1)*100
for i=0, 49-1 do pure_ps2(i,*) = ((sea_ps2(i,*) - mcd_ps_12_2) / mcd_ps_12_2)*100

for i=0, 49-1 do pure_ps3(i,*,*) = ((sea_ps3(i,*,*) - mcd_ps_12_3) / mcd_ps_12_3)*100

; OMEGAのリトリーバル範囲（±70°）に限定する
Lat_arr = Lat_arr(ind)
pure_ps3 = pure_ps3(*,*,ind)
arr_index = n_elements(Lat_arr)

minp6 = 0
maxp6 = 360

minp7 = -90
maxp7 = 90

; plotの枠組み
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

window, 0, xs=1500, ys=1000
loadct,39
;plot, findgen(100),xs=1, ys=1, yr=[400, 650], xr=[0, 360], back=255, color=0, /nodata, charsize=3, title='tidal wave detection [MCD]', ytitle='mean pressure [Pa]'
;plot, findgen(100),xs=1, ys=1, yr=[400, 550], xr=[0, 24], back=255, color=0, /nodata, charsize=3, title='tidal wave detection [MCD]', ytitle='mean pressure [Pa]'
plot, findgen(100),xs=1, ys=1, yr=[-7,7], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidal wave detection Ls 270', xtitle='Local time', ytitle='mean pressure [%]'
cgLOADCT, 39
;cgColorbar, Divisions=4, Minor=5, charsize=2, position=[0.063,0.12,0.067,0.86], Range=[minp6,maxp6],title='Ls [deg]', TCHARSIZE=10, /vertical
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp7,maxp7],title='Latitude [deg]', position=[0.063,0.12,0.067,0.86], TCHARSIZE=10, /vertical 

;for i =0, 25-1 do begin
  for j =0, arr_index-1 do begin
    loadct, 39
    ;plots, LT_arr, mcd_ps1(*,0), color=0, psym=8, symsize=2
    ;plots, LT_arr, mcd_ps2(*,0), color=254, psym=8, symsize=2
    ;plots, LS_arr(i), sea_ps1(24,i), color=0, psym=8, symsize=2
    ;plots, LT_arr(10:35), pure_ps(10:35,i), color=(LS_arr(i)/(maxp6-minp6))*254., psym=8, symsize=2
    ;plots, LT_arr(*), pure_ps3(*,i,j), color=(LS_arr(i)/(maxp6-minp6))*254., psym=8, symsize=2
    plots, LT_arr(15:33), pure_ps3(10:33,18,j), color=((Lat_arr(j)+maxp7)/(maxp7-minp7))*254., psym=8, symsize=2
    ;plots, LT_arr(10:35), sea_ps3(10:35,0,j), color=(Lat_arr(j)/(maxp7-minp7))*254., psym=8, symsize=2
  endfor
;endfor


stop

snapshot = TVRD(True=1)
Write_JPEG, path_ql + 'Latitude_ls0.jpg', snapshot, True=1, Quality=75

stop


; MOLA情報がOMEGAデータで大丈夫かどうかの確認をしてみる
; まずはOMEHAの高度を使用してMCDのPSの高度補正をしたマップを作成してみる

; ===================== resrore working file =====================
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'
path2 = '/work1/LUT/SP/mcddata/'

; 取り敢えず1ORBで確かめてみる
file = '/work1/LUT/SP/EWwork/EW_work_ORB1050_4.sav'
restore, file

; ===== altitude correction ======
g = 3.72 ;m/s2
Gconst = 6.67430d-11
MMars = 6.42d23
RMars = 3.4d6
g = -Gconst*MMars/(-RMars*RMars)
R = 192d ;J/K/kg

; n_indを格納
ip = n_elements(LATI(*,0))
io = n_elements(LATI(0,*))

; longitude, latitude
lat1 = lati
lon1 = longi

; ls
SOLAR_LONGITUDE = SOLAR_LONGITUDE
local_time = local_time

; MCD surface pressure
p1_mcd = mcdpressure 

; bad pixelにNanを格納
for i = 0, ip-1 do for j = 0, io-1 do if p1_mcd(i,j) eq 0 then p1_mcd(i,j) = !VALUES.F_NAN

; =============== mcd pressure altutude corecction using OMEGA data ===============
p1 = dblarr(ip,io)
p5 = altitude
TAmap = TAMAP

for i = 0, ip-1 do for j = 0, io-1 do if p1_mcd(i,j) gt 0 then p1(i,j) = p1_mcd(i,j) * exp( (p5(i,j)) / (R*TAmap(i,j)/(-Gconst*MMars/(-1d*(RMars+p5(i,J))*(RMars+p5(i,J)) ))))
for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN


; =============== mcd pressure altutude corecction using MCD data ===============
MOLA_mcd = dblarr(ip,io)
p2 = dblarr(ip,io)

for l = 0, ip -1 do begin ;loop for slit scan
    for k = 0, io -1 do begin ;test getting surface feature

      ; =========== MCD version 6.1 ===========
        lat = lati(l,k)
        lon = longi(l,k)
        Ls = SOLAR_LONGITUDE
        Loct = local_time
        hrkey = 1 ; set high resolution mode on (hrkey=0 to set high resolution off)
        zkey = 3    ; specify that xz is the altitude above surface (m)
        xz = 0. ; (array of altitude: step=2000 m start: 5m)
        datekey = 1       ; <integer> type of input date (1=Mars date)
        xdate = Ls       ; <double precision> date (IF datekey = 1 : Value of Ls [deg.])
        dset ='/work1/LUT/MCDv6.1/MCD_6.1_queleclim/data/'
        scena = 27         ; <integer> scenario (1 = Climatology ave solar)
        perturkey = 1     ; <integer>  perturbation type (1= none)
        seedin   = 7.0    ; <real>
        gwlength = 0.0    ; <real>  for small scale (ie: gravity wave) perturbations;
        extvarkeys = LONARR(101) ; add 1 element because indexes in IDL start at 0
        for i0 = 0, 100 do extvarkeys[i0] = 1 ; <integer> array output type (extvar(i) = 0 : don’t compute)
        meanvar = FLTARR(6) ; add 1 element because indexes in IDL start at 0
        for i0 = 0, 5 do meanvar[i0] = 0.0 ; <real> mean unperturbed values (array of 5)
        extvar = FLTARR(101) ; add 1 element because indexes in IDL start at 0
        for i0 = 0, 100 do extvar[i0] = 0.0 ; <real>  extra variables (array of 100)
        pres = 0.0        ; <real> atmospheric pressure (Pa)
        dens = 0.0        ; <real> atmospheric density (kg/m^3)
        temp = 0.0        ; <real> atmospheric temperature (K)
        zonwind = 0.0     ; <real> zonal wind component (East-West)
        merwind = 0.0     ; <real> meridional wind component (North-South)
        seedout = 0.0     ; <real> current value of the seed of the random number generator
        ierr = 0          ; <integer> error flag (0 = OK)

        cd, '/work1/LUT/MCDv6.1/MCD_6.1_queleclim/mcd/interfaces/idl/'
        a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
        dset,scena,perturkey,seedin,gwlength,extvarkeys, $
        pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)

        MOLA_altitude = extvar(4)

        ; +++++++++ save MCD data +++++++++++
        MOLA_mcd(l,k) = MOLA_altitude

    endfor
endfor

for i = 0, ip-1 do for j = 0, io-1 do if p1_mcd(i,j) gt 0 then p2(i,j) = p1_mcd(i,j) * exp( (MOLA_mcd(i,j)) / (R*TAmap(i,j)/(-Gconst*MMars/(-1d*(RMars+MOLA_mcd(i,J))*(RMars+MOLA_mcd(i,J)) ))))
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 0 then p2(i,j) = !VALUES.F_NAN


; ================= plot 準備 ====================
lat1 = lati
minlati = min(lati)
maxlati = max(lati)
lon1 = longi
minlon = min(longi)
maxlon = max(longi)

; using OMEGA data
minp1 = min(p1)
maxp1 = max(p1)

; using mcd data
minp2 = min(p2)
maxp2 = max(p2)

; OMEGAとmcd dataの差
p3 = p1 - p2
minp3 = min(p3)
maxp3 = max(p3)

; MOLA map using OMEGA
minp5 = min(p5)
maxp5 = max(p5)
medip5 = median(p5)

; MOLA map using OMEGA
p4 = MOLA_mcd
minp4 = min(p4)
maxp4 = max(p4)

p6 = p5 -p4
minp6 = min(p6)
maxp6 = max(p6)


; ==================== plot　開始 ==============================
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

window, 22, xs=1800, ys=1600
!P.Multi = [0, 3, 2]

; altitude correction using OMEGA data
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.55,0.31,0.92], /nodata, charsize=3, title='using OMEGA data', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.04,0.58,0.046,0.86],title='Pa', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1

; MCD prediction pressure
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.55,0.65,0.92], /nodata, charsize=3, title='using MCD data', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
;cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp2,maxp2],position=[0.38,0.58,0.386,0.86],title='Pa', TCHARSIZE=10, /vertical 
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.38,0.58,0.386,0.86],title='Pa', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p2(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1

; retrieved - MCD
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.55,0.99,0.92], /nodata, charsize=3, title='OMEGA - MCD', xtitle='Lon', ytitle='Lat'
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp3,maxp3],position=[0.72,0.58,0.726,0.86],title='Pa', TCHARSIZE=10, /vertical 
for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p3(i,j)-minp3)/(maxp3-minp3)*254., psym=8, symsize=1

; input MOLA map using OMEGA
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.08,0.31,0.45], /nodata, charsize=3, title='using OMEGA', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp5,maxp5],position=[0.04,0.12,0.046,0.4],title='m', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p5(i,j)-minp5)/(maxp5-minp5)*254., psym=8, symsize=1

; dust opacity map
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.08,0.65,0.45], /nodata, charsize=3, title='using mcd data', xtitle='Lon', ytitle='Lat'
cgLOADCT, 16
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp4,maxp4],position=[0.38,0.12,0.386,0.4],title='m', TCHARSIZE=10, /vertical 
loadct, 16
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p4(i,j)-minp4)/(maxp4-minp4)*254., psym=8, symsize=1

; MOLA altitude map
loadct,39
plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.45], /nodata, charsize=3, title='OMEGA-mcd', xtitle='Lon', ytitle='Lat'
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp6,maxp6],position=[0.72,0.12,0.726,0.4], TCHARSIZE=10, /vertical 
loadct, 39
for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p6(i,j)-minp6)/(maxp6-minp6)*254., psym=8, symsize=1

print, mean(p1, /NAN)
print, mean(p2, /NAN)
print, stddev(p1, /NAN)
print, stddev(p2, /NAN)

stop


end