Pro nyooon

; IDLの思考整理.pro file
; 試したいことを色々試せるfile


for ISZA = 1, 2 do begin ;4) SZA
  if ISZA eq 1 then SZA = 00.0
  if ISZA eq 2 then SZA = 15.0
  
  for IPA = 2,2 do begin ;6) Phase angle
    if IPA eq 1 then PA = 00.0
    if IPA eq 2 then PA = 45.0
    if IPA eq 3 then PA = 90.0
      
      print, SZA + PA
      
  endfor
endfor


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