;------------------------------------------------------------------------------------------------------------------------
Pro Cal_LUtable_obsspec_ver2
; +++
; to read outputs from ARS and calculate equivalent width (=absorption depth) for each spectrum
; create by Shohei Aoki
; 
; edit by Akira Kazama
; 
; cal_lutable_obsspec　 ::2022.6.1 Wed 11:43:00
; OMEGAで観測されたスペクトルの吸収深さ量を計算するプログラム
; 
; cal_lutable_obsspec_ver2  ::2022.06.07 Tue 11:10:00
; I/F -> 放射輝度から吸収量を計算するプログラム
; 
; +++
;------------------------------------------------------------------------------------------------------------------------

; ============ file情報 ==========================
; 観測データ
path = '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0931_3.sav'
restore, path

; ============== CO2吸収量の計算を行う ==============================

; CO2の吸収帯を持ってくる
CO2=where(wvl gt 1.81 and wvl lt 2.19)
wvl=wvl[CO2]

nwvl=n_elements(wvl)
io=n_elements(jdat(*,1,1))
ip=n_elements(jdat(1,1,*))
flux=dblarr(io,nwvl,ip)

for i=0,io-1 do begin
  for o=0,ip-1 do begin
    flux(i,*,o)=jdat(i,co2,o)
  endfor
endfor

; 全ての緯度・経度を持ってきて、plotさせる
maxlongi = max(longi)
minlongi = min(longi)
maxlati = max(lati)
minlati = min(lati)

ind=where_xyz(longi ge minlongi and longi le maxlongi and lati ge minlati and lati le maxlati,xind=xind,yind=yind)

; CO2の吸収線の中
band=where(wvl gt 1.9 and wvl lt 2.1)

;where_xyzで探してきたindを入れると、その場所のlongi,latiを取り出してくることができる。今は全範囲。
longi1=longi(ind)
lati1=lati(ind)

;データ番号と波長のデータセットを作成する
radiance=dblarr(n_elements(ind),nwvl)

;reformは次元を落としてくれるよ
for i=0,nwvl-1 do begin
  newflux1=reform(flux(*,i,*))
  radiance(*,i)=double(newflux1(ind))
endfor

;OMEGAの38番目と48番目の素子が死んでた
nanserch=where(radiance ge 0 and radiance le 0.0001)
radiance(nanserch)=!VALUES.F_NAN

; ======== CO2吸収量 mapping =============
obs_spec = dblarr(n_elements(ind))

for i=0, 75007 do begin
  x = [wvl(0),wvl(3),wvl(5),wvl(23),wvl(24),wvl(25)]
  y = [radiance(i,0), radiance(i,3), radiance(i,5), radiance(i,23), radiance(i,24), radiance(i,25)]
  coef = linfit(x,y)
  cont = coef(0) + coef(1)*wvl
  width = 1- radiance(i,*)/cont
  
  obs_spec(i) = total(width[band],/nan)  
endfor

;save,obs_spec,$
;  filename='/Users/nyonn/IDLWorkspace/Default/savfile/Table_SP_obs_calc_orb0931_3.sav'
;stop
END