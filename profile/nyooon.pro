function gauss, x, myu, sigma

gauss_func = 1/(sigma*sqrt(2*!pi))*exp(-(x-myu)^2/2*sigma^2)

return, gauss_func

end

Pro nyooon

; file 検索方法
Result = file_test(output_name_rad)
if result eq 0 then begin



; binary fileを読む
ind = 0l
h = 0d
hi = fltarr(31)
t = 0d
ti = fltarr(31)
p = 0d
pi = fltarr(31)
Kw = fltarr(31,101101)
Kw2 = fltarr(31,101101)

file1 = '/Users/nyonn/IDLWorkspace/Default/savfile/LUTable_T1_285_T2_200_PRS1500.k'
file2 = '/Users/nyonn/IDLWorkspace/Default/savfile/CO2_SP15_TA5_TB3.k'

openr, 1, file1
; for i=0,31 -1  do readu, 1, ind, h, t, p, Kw(*,i)
for i =0, 30 do readu, 1, ind, hi(i), ti(i), pi(i), Kw(*,i)
close,1

openr, lun, file1, /get_lun
fs = fstat(lun)

len = fs.size / n_bytes_in_data_structure
for i = 0L, len - 1 do begin
  readu, lun, var


restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0030_1.sav'


; ガウシアンの定義
x = reform(wvl(0:127))
myu2 = 1.9d
myu3 = 2.3d
sigma2 = 0.5d / (2* sqrt(2*alog(2)))
sigma3 = 0.56d / (2* sqrt(2*alog(2)))

y2 = gauss1(x, [1.9d, 0.5d / (2 * sqrt(2*alog(2)))])
y3 = gauss1(x, [1.9d, 0.56d / (2 * sqrt(2*alog(2)))])

y = y2 + y3

gauss_2 = gauss(x,myu2,sigma2)
gauss_3 = gauss(x,myu3,sigma3)

gaussian = gauss_2 + gauss_3




Pressure_grid = dblarr(15)
Pressure_grid(0) = alog(50d)
Pressure_grid(1) = alog(150d)
Pressure_grid(2) = alog(180d)
Pressure_grid(3) = alog(215d)
Pressure_grid(4) = alog(257d)
Pressure_grid(5) = alog(308d)
Pressure_grid(6) = alog(369d)
Pressure_grid(7) = alog(442d)
Pressure_grid(8) = alog(529d)
Pressure_grid(9) = alog(633d)
Pressure_grid(10) = alog(758d)
Pressure_grid(11) = alog(907d)
Pressure_grid(12) = alog(1096d)
Pressure_grid(13) = alog(1300d)
Pressure_grid(14) = alog(1500d)

min_ind = 5
total_LMS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.43, 0.14, 0.9, 1.0, 0.54]

before_pres = (Pressure_grid(min_ind - 1) + Pressure_grid(min_ind))/2
after_pres = (Pressure_grid(min_ind + 1) + Pressure_grid(min_ind))/2

before_interpol = interpol(total_LMS, Pressure_grid, before_pres)
after_interpol = interpol(total_LMS, Pressure_grid, after_pres)

print, "before:", before_interpol
print, "after:", after_interpol

result_compare = before_interpol < after_interpol

print, "result:", result_compare

repeat begin
if before_interpol eq result_compare then begin

  loop_ind_1 = (Pressure_grid(min_ind - 1) + before_pres)/2
  loop_ind_2 = (Pressure_grid(min_ind) + before_pres)/2
  
  print, "loop1:", loop_ind_1
  print, "loop2:", loop_ind_2
  
  before_interpol = interpol(total_LMS, Pressure_grid, loop_ind_1)
  after_interpol = interpol(total_LMS, Pressure_grid, loop_ind_2)
  
  result_compare = before_interpol < after_interpol
  print, "loop_result:", result_compare

endif else  begin

  loop_ind_1 = (Pressure_grid(min_ind + 1) + after_pres)/2
  loop_ind_2 = (Pressure_grid(min_ind) + after_pres)/2
    
  before_interpol = interpol(total_LMS, Pressure_grid, loop_ind_1)
  after_interpol = interpol(total_LMS, Pressure_grid, loop_ind_2)
  
  result_compare = before_interpol < after_interpol
  print, "loop_result:", result_compare

endelse
endrep until result_compare le 0.0001

;aa = [1,2,3]
;
;for i =0, 3 do begin
;  if i eq 0 then I1 = aa + 1
;  if i eq 1 then I2 = aa+ 2
;  if i eq 3 then I3 = aa +3
;endfor
;total_I = I1+I2+I3
;print, I1
;print, I2
;print, total_I

; IDLの思考整理.pro file
; 試したいことを色々試せるfile

; Lsとind場所を知りたいときに使うコード
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0363_0.sav'
;Ls=strmid(SOLAR_LONGITUDE,5,7)

; Densityで振ってみることを試してみる
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0030_1.sav'
ind=where_xyz(longi ge 60.79 and longi le 60.81 and lati ge -48.44 and lati le -48.43,xind=xind,yind=yind)
; ind=32015, lati: 48.431797 S, longi: 60.808998 E [Forget+, retrievalすると1036 Pa]
; 
; restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0363_3.sav'
; ind=where_xyz(longi ge 311.73 and longi le 311.78 and lati ge 22.7 and lati le 22.72,xind=xind,yind=yind)
; ind=76094,xind=62, yind=594, lati:22.705700 N , longi:311.76300 E (48.237 W)  [Forget+, retrievalすると852 Pa]

; restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB1201_3.sav'
; ind=where_xyz(longi ge 24.979 and longi le 24.981 and lati ge -7.765 and lati le -7.763,xind=xind,yind=yind)
; ind=37135, xind=15, yind=1160, lati:7.764°S, longi:24.980°E　[Forget+, retrievalすると470 Pa]

; restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0931_3.sav'
; ind=where_xyz(longi ge 276.50 and longi le 276.51 and lati ge 51.05 and lati le 51.1,xind=xind,yind=yind)
; ind = 68287, xind=63, yind=533, lati:51.068897 N, longi:276.50281 E, 青木さんの結果を再現

openr,2,'/Users/nyonn/IDLWorkspace/Default/profile/specsol_0403.dat'
specmars=0B
;;spcmarsに入れるOMEGAの波長分
;格納する場所：specmars
specmars=fltarr(352)
readf,2,specmars
close,2
specmars = specmars/dmars/dmars

wvl_ind = where(wvl gt 0.35 and wvl lt 1.0)
wvl=wvl(wvl_ind)
specmars = specmars(wvl_ind)

io=n_elements(jdat(*,1,1))
ip=n_elements(jdat(1,1,*))
nwvl=n_elements(wvl)

flux=dblarr(io,nwvl,ip)
for i=0,io-1 do begin
  for o=0,ip-1 do begin
    flux(i,*,o)=jdat(i,wvl_ind,o)/specmars
  endfor
endfor

radiance=dblarr(n_elements(ind),nwvl)
for i=0,nwvl-1 do begin
  newflux=reform(flux(*,i,*))
  radiance(*,i)=double(newflux(ind))
endfor

nanserch=where(radiance ge -1 and radiance le 0.0001)
radiance(nanserch)=!VALUES.F_NAN

plot, wvl, radiance,xstyle=1



; curve fittingを試す！
path1 = '/Users/nyonn/IDLWorkspace/Default/LUT/SP1_TA1_TB1_SZA1_EA1_PA1_Dust1_WaterI1_SurfaceA1_rad.dat'
path2 = '/Users/nyonn/IDLWorkspace/Default/LUT/SP1_TA1_TB1_SZA1_EA3_PA1_Dust5_WaterI2_SurfaceA6_rad.dat'

wn=dblarr(27)
rad1=dblarr(27)
rad2=dblarr(27)


;read files
openr,lun,path1,/get_lun
for i = 0, 27-1 do begin
  readf,lun,a,b
  wn(i)=a
  rad1(i)=b
endfor
free_lun,lun

openr,lun,path2,/get_lun
for i = 0, 27-1 do begin
  readf,lun,a,b
  rad2(i)=b
endfor
free_lun,lun

wn = (1/wn)*10000
wn = reverse(wn)

rad1 = reverse(rad1)
rad2 = reverse(rad2)


stop

; ここのモデル部分を気圧のパラメータを振れば変わるものを入れれたら良い。
; 圧力の関数にすれば良いから、、P=index_0のときはrad0 → 多次元補間されたものを入れればいい
expr = rad2
result = mpfitexpr(expr, wn, rad1,rad2)

stop

; nスキャンは何個？？
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0030_1.sav'
;openr,2,'/Users/nyonn/IDLWorkspace/Default/profile/specsol_0403.dat'
;specmars=0B
;specmars=fltarr(352)
;readf,2,specmars
;close,2
;
;specmars = specmars/dmars/dmars
;
;CO2=where(wvl gt 1.81 and wvl lt 2.19)
;wvl=wvl[CO2]
;jdat=jdat(*,CO2,*)
;specmars = specmars(CO2)
;
;span_lati = max(lati(0,*)) - min(lati(0,*))
;nlay_span = ceil(span_lati/1.875)
;ip = n_elements(LATI(*,0))
;
;for n = 0, 1 do begin
;  points = where(lati(0,*) ge min(lati(0,*)) + float(n)*1.875 and lati(0,*) lt min(lati(0,*)) + float(n+1)*1.875, count)
;  count_nscan = count
;  
;  for i = 0, 5 do begin ;loop for slit scan
;    for j = 0, 10 do begin
;      Albedo_input = jdat(j,0,points(i))/ specmars(0) / cos(geocube(j,8,points(i))*1e-4*!DTOR)
;      print,n,i,j
;      print,'Albedo:', Albedo_input
;      print, 'specmars:', specmars(0)
;      print, 'cos:', cos(geocube(j,8,points(i))*1e-4*!DTOR)
;    endfor
;  endfor
;endfor
;

;; 水和鉱物がある場所のスペクトルを見たい！
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0232_0.sav'
;;specsol_0403:太陽輝度情報
openr,2,'/Users/nyonn/IDLWorkspace/Default/profile/specsol_0403.dat'
specmars=0B
;;spcmarsに入れるOMEGAの波長分
;格納する場所：specmars
specmars=fltarr(352)
readf,2,specmars
close,2

specmars = specmars/dmars/dmars
Ls=strmid(SOLAR_LONGITUDE,5,7)

wvl_ind = where(wvl gt 1.0 and wvl lt 2.6)
wvl=wvl(wvl_ind)
specmars = specmars(wvl_ind)

io=n_elements(jdat(*,1,1))
ip=n_elements(jdat(1,1,*))
nwvl=n_elements(wvl)

flux=dblarr(io,nwvl,ip)
for i=0,io-1 do begin
  for o=0,ip-1 do begin
    flux(i,*,o)=jdat(i,wvl_ind,o)/specmars
  endfor
endfor

; 大シルチス台地
ind=where_xyz(longi ge 68.7 and longi le 69.3 and lati ge 7.8 and lati le 8.4,xind=xind,yind=yind)
radiance=dblarr(n_elements(ind),nwvl)

for i=0,nwvl-1 do begin
  newflux=reform(flux(*,i,*))
  radiance(*,i)=double(newflux(ind))
endfor

nanserch=where(radiance ge 0 and radiance le 0.0001)
radiance(nanserch)=!VALUES.F_NAN

plot, wvl, radiance,xrange=[1.1,2.7],xstyle=1
; plot, wvl, radiance,yrange=[0.13,0.19]

; for文を回さない
;for ISZA = 1, 2 do begin ;4) SZA
;  if ISZA eq 1 then SZA = 00.0
;  if ISZA eq 2 then SZA = 15.0
;  
;  for IPA = 2,2 do begin ;6) Phase angle
;    if IPA eq 1 then PA = 00.0
;    if IPA eq 2 then PA = 45.0
;    if IPA eq 3 then PA = 90.0
;      
;      print, SZA + PA
;      
;  endfor
;endfor


; 吸収
;wn=dblarr(27)
;rad1=dblarr(27)
;
;file1 = '/Users/nyonn/IDLWorkspace/Default/LUT/SP1_TA1_TB1_SZA1_EA1_PA1_Dust1_WaterI1_SurfaceA1_rad.dat'
;openr,lun,file1,/get_lun
;for i = 0, 27-1 do begin
;  readf,lun,a,b
;  wn(i)=a
;  rad1(i)=b
;endfor
;
;wn = (1/wn)*10000
;wn = reverse(wn)
;band=where(wn gt 1.89 and wn lt 2.05)
;
;rad1 = reverse(rad1)
;
;x = [wn(0), wn(3), wn(5), wn(23), wn(24), wn(25)]
;y1 = [rad1(0), rad1(3), rad1(5), rad1(23), rad1(24), rad1(25)]
;coef1 = linfit(x,y1)
;
;cont1 = coef1(0) + coef1(1)*wn
;width1 = 1.0 - rad1/cont1
;y_calc = total(width1[band])


path =  '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0920_3.sav'
restore, path

; ARSで計算した場所の波長帯を取ってきている
CO2=where(wvl gt 1.81 and wvl lt 2.19)
wvl=wvl[CO2]

nwvl=n_elements(wvl)

io=n_elements(jdat(*,1,1))
ip=n_elements(jdat(1,1,*))

ip_1 = n_elements(LATI(*,0))
io_1 = n_elements(LATI(0,*))

obs_spec = reform(jdat(*,0,*))

; CO2の吸収線の中
band=where(wvl gt 1.94 and wvl lt 2.09)

;OMEGAの38番目と48番目の素子が死んでた
nanserch=where(jdat ge 0 and jdat le 0.0001)
jdat(nanserch)=!VALUES.F_NAN

jdat = jdat(*,CO2,*)

; ======== CO2吸収量 mapping =============
for i=0, io-1 do begin
  for j=0, ip-1 do begin
;for i=0, 2 do begin
;  for j=0, 2 do begin
    x = [wvl(0),wvl(1),wvl(2),wvl(23),wvl(24),wvl(25)]
    y = [jdat(i,0,j), jdat(i,1,j), jdat(i,2,j), jdat(i,23,j), jdat(i,24,j), jdat(i,25,j)]
    coef = linfit(x,y)
    cont = coef(0) + coef(1)*wvl
    width = 1- jdat(i,*,j)/cont
    
    SZA = reform(geocube(i,8,j))*1.e-4
    EA = reform(geocube(i,9,j))*1.e-4
    PA = reform(geocube(i,10,j))*1.e-4
  
    obs_spec(i,j) = total(width[*,band],/nan)
  endfor
endfor

;x = [wvl(0),wvl(3),wvl(5),wvl(23),wvl(24),wvl(25)]
;y = [jdat(0,0,0), jdat(0,3,0), jdat(0,5,0), jdat(0,23,0), jdat(0,24,0), jdat(0,25,0)]
;coef = linfit(x,y)
;cont = coef(0) + coef(1)*wvl
;width = 1- jdat(0,*,0)/cont
;
;obs_spec = total(width[*,band],/nan)


print, obs_spec


;path =  '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0931_3.sav'
;restore, path
;
;; ARSで計算した場所の波長帯を取ってきている
;CO2=where(wvl gt 1.81 and wvl lt 2.19)
;wvl=wvl[CO2]
;
;nwvl=n_elements(wvl)
;
;io=n_elements(jdat(*,1,1))
;ip=n_elements(jdat(1,1,*))
;
;; flux=dblarr(io,nwvl,ip)
;flux = reform(jdat(*,0,*))
;
;; 全ての緯度・経度を持ってきて、plotさせる
;maxlongi = max(longi)
;minlongi = min(longi)
;maxlati = max(lati)
;minlati = min(lati)
;
;ind=where_xyz(longi ge minlongi and longi le maxlongi and lati ge minlati and lati le maxlati,xind=xind,yind=yind)
;
;; CO2の吸収線の中
;band=where(wn gt 1.89 and wn lt 2.05)
;
;;where_xyzで探してきたindを入れると、その場所のlongi,latiを取り出してくることができる。今は全範囲。
;longi1=longi(ind)
;lati1=lati(ind)
;
;;データ番号と波長のデータセットを作成する
;radiance=dblarr(n_elements(ind),nwvl)
;
;;reformは次元を落としてくれるよ
;for i=0,nwvl-1 do begin
;  newflux1=reform(flux(*,i,*))
;  radiance(*,i)=double(newflux1(ind))
;endfor
;
;;OMEGAの38番目と48番目の素子が死んでた
;nanserch=where(radiance ge 0 and radiance le 0.0001)
;radiance(nanserch)=!VALUES.F_NAN
;
;; ======== CO2吸収量 mapping =============
;obs_spec = dblarr(n_elements(ind))
;
;for i=0, 76288-1 do begin
;  x = [wvl(0),wvl(3),wvl(5),wvl(23),wvl(24),wvl(25)]
;  y = [radiance(i,0), radiance(i,3), radiance(i,5), radiance(i,23), radiance(i,24), radiance(i,25)]
;  coef = linfit(x,y)
;  cont = coef(0) + coef(1)*wvl
;  width = 1- radiance(i,*)/cont
;
;  obs_spec(i) = total(width[band],/nan)
;endfor


;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;ip = n_elements(lati(*,0))
;io = n_elements(lati(0,*))
;    
;trans = reform(jdat(*,0,*))
;pressure = trans
;  
;CO2=where(wvl gt 1.81 and wvl lt 2.19)
;wvl=wvl[CO2]
;jdat=jdat(*,CO2,*)
;
;band=where(wvl gt 1.9 and wvl lt 2.1)
;
;x = [wvl(0),wvl(3),wvl(5),wvl(23),wvl(24),wvl(25)]
;
;;latitude grid
;span_lati = max(lati(0,*)) - min(lati(0,*))
;nlay_span = ceil(span_lati/1.875)
;
;; MCDのグリッドが3.75°だから、1.875°ごとにMCDの読み出しを行う
;; 基本的にはOMEGAは緯度方向にスキャンを行う
;
;for i = 0, io-1 do begin ;loop for slit scan
;  for j = 0, ip-1 do begin  
;
;    Y = [jdat(j,0,i), jdat(j,3,i), jdat(j,5,i), jdat(j,23,i), jdat(j,24,i), jdat(j,25,i)]
;    coef = linfit(X,Y)
;    cont = coef(0) + coef(1)*wvl
;
;    SZA = reform(geocube(j,8,i))*1.e-4
;    EA = reform(geocube(j,9,i))*1.e-4
;    PA = reform(geocube(j,10,i))*1.e-4
;    
;    ; Albedo_input = jdat(j,0,points(i))/specmars(0) / cos(geocube(j,8,points(i))*1e-4*!DTOR)
;    width = 1.0 - jdat(j,*,i)/cont
;    trans(j,i) = total(width[band], /nan)
;    ; pressure(j,points(i)) = ret_pressure(trans(j,points(i)), TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input)
;
;  endfor
;endfor
;
;
;loadct, 39
;Set_Plot, 'Z'
;Device, Set_Resolution=[1000,800], Set_Pixel_Depth=24, Decomposed=0
;
;!p.multi=[0,2,1]
;
;min_pressure = min(trans(*,*),/nan) ;0.
;max_pressure = max(trans(*,points),/nan) ;0.1; max(trans(*,where(flag(points) eq 0)),/nan)*1.1
;rn_pressure = max_pressure - min_pressure
;plot,longi(*,*),lati(*,*),xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
;  title='Ls:'+STRCOMPRESS(Solar_longitude)+' Max:'+string(max_pressure,format='(f4.2)' ),charsize=1.
;for i = 0, io-1 do begin
;  for j = 0, ip-1 do begin
;    color = (trans(j,i)-min_pressure)/rn_pressure*254.
;    if color gt 254 then color = 254
;    if color lt 0 then color = 0
;    plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
;  endfor
;endfor
;
;!Y.OMargin = [0,0]
;!p.multi=0
;snapshot = TVRD(True=1)
;Write_JPEG, 'test.jpg', snapshot, True=1, Quality=75
;
;stop

end