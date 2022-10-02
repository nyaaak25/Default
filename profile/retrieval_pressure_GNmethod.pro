; まずは1点でfittingをかけてみる
Pro retrieval_pressure_GNmethod
; fitting test

path = '/data2/omega/sav/'
path2 = '/work1/LUT/SP/table/absorption/'

;　観測したいデータを選択
restore, path+'ORB0363_3.sav'
restore, path + 'specmars.sav'

; CO2 absorption emission line   
CO2=where(wvl gt 1.81 and wvl lt 2.19)
wvl=wvl[CO2]

jdat=jdat(*,CO2,*)
specmars = specmars(CO2)

; 観測地点の選択
; ind=32015, xind=15, yind=1000, lati: 48.431797 S, longi: 60.808998 E [Forget+, retrievalすると1036 Pa] ORB0030_1
; ind=where_xyz(longi ge 60.79 and longi le 60.81 and lati ge -48.44 and lati le -48.43,xind=xind,yind=yind)

; ind=76094,xind=62, yind=594, lati:22.705700 N , longi:311.76300 E (48.237 W)  [Forget+, retrievalすると852 Pa] ORB0363_3
ind=where_xyz(longi ge 311.73 and longi le 311.78 and lati ge 22.7 and lati le 22.72,xind=xind,yind=yind)

; ind=37135, xind=15, yind=1160, lati:7.764°S, longi:24.980°E　[Forget+, retrievalすると470 Pa] ORB1201_3
; ind=where_xyz(longi ge 24.979 and longi le 24.981 and lati ge -7.765 and lati le -7.763,xind=xind,yind=yind)

; ind = 68287, xind=63, yind=533, lati:51.068897 N, longi:276.50281 E, 青木さんの結果を再現 ORB0931_3
; ind=where_xyz(longi ge 276.50 and longi le 276.51 and lati ge 51.05 and lati le 51.1,xind=xind,yind=yind)
; ind = where_xyz(longi eq 275.82278 and lati eq 52.079899,xind=xind,yind=yind)


nanserch=where(jdat ge 0 and jdat le 0.0001)
jdat(nanserch)=!VALUES.F_NAN

;    ; ----- MCD ------
lat = median(lati(xind,yind),/double)
lon = median(longi(xind,yind),/double)
Ls = SOLAR_LONGITUDE
Loct = 12. + (lon - SUB_SOLAR_LONGITUDE)/15. ;! TBC!!! !
dust = 2  ; our best guess MY24 scenario, with solar average conditions
hrkey = 1 ; set high resolution mode on (hrkey=0 to set high resolution off)
zkey = 3    ; specify that xz is the altitude above surface (m)
xz = 0. ; (array of altitude: step=2000 m start: 5m)
datekey = 1       ; <integer> type of input date (1=Mars date)
xdate = ls        ; <double precision> date (IF datekey = 1 : Value of Ls [deg.])
dset ='/work1/LUT/MCD/MCD5.3/data/' ;‘MCD_DATA/’  ; <character*50> data set
scena = 1         ; <integer> scenario (1 = Climatology ave solar)
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
;    cd, ‘/Users/Shohei/Documents/Programs/MCD_ver4.3/mcd/idl/’
;   mcd_idl, Ls, Loct, lat, lon, dust, zkey, hrkey, xz, meanvarz, extvarz
cd, '/work1/LUT/MCD/MCD5.3/mcd/idl/'
a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
  dset,scena,perturkey,seedin,gwlength,extvarkeys, $
  pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
scaleH = extvar(13)
Z1 = scaleH * 0.1d
Z2 = scaleH * 4d
dust_opacity = extvar(36)*1.2090217E+00/4.6232791E+00
ice_opacity = 0.d;*4.2410289E-01/***  ;!TBD!
SP = PRES
xz=[Z1] ; (array of altitude with step=2000 m)
a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
  dset,scena,perturkey,seedin,gwlength,extvarkeys, $
  pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
TA = temp
xz=[Z2] ; (array of altitude with step=2000 m)
a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
  dset,scena,perturkey,seedin,gwlength,extvarkeys, $
  pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
TB = temp
;;-----------------

x = [wvl(0),wvl(1),wvl(2),wvl(23),wvl(24),wvl(25)]
Y = [jdat(xind,0,yind), jdat(xind,1,yind), jdat(xind,2,yind), jdat(xind,23,yind), jdat(xind,24,yind), jdat(xind,25,yind)]
coef = linfit(X,Y)
cont = coef(0) + coef(1)*wvl

SZA = reform(geocube(xind,8,yind))*1.e-4
EA = reform(geocube(xind,9,yind))*1.e-4
PA = reform(geocube(xind,10,yind))*1.e-4

Albedo_input = jdat(xind,0,yind)/specmars(0) / cos(geocube(xind,8,yind)*1e-4*!DTOR)

width = 1.0 - jdat(xind,*,yind)/cont
Intensity = width

for wave = 0,26 do begin

  if wave eq 0 then LMS0 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 1 then LMS1 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 2 then LMS2 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 3 then LMS3 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 4 then LMS4 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 5 then LMS5 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 6 then LMS6 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 7 then LMS7 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 8 then LMS8 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 9 then LMS9 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 10 then LMS10 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 11 then LMS11 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 12 then LMS12 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 13 then LMS13 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 14 then LMS14 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 15 then LMS15 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 16 then LMS16 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 17 then LMS17 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 18 then LMS18 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 19 then LMS19 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 20 then LMS20 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 21 then LMS21 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 22 then LMS22 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 23 then LMS23 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 24 then LMS24 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 25 then LMS25 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 26 then LMS26 = ARS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)

endfor

; GN methodを実行
initial = alog(400) ; 初期値の設定、400 Pa ~ 800 Paくらいで様子を見る
pressure_ID = GNmethod(LMS_calc, initial)

print(exp(pressure_ID))



end

