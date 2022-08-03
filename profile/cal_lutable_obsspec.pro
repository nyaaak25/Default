;------------------------------------------------------------------------------------------------------------------------
Pro cal_lutable_obsspec

;
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
; cal_lutable_obsspec  :: 2022.07.11 Mon 11:32:00
; version3  二次元に配列を落とさないで、jdat(*,CO2,*)の形でtransを返すプログラム
; 
; +++
;------------------------------------------------------------------------------------------------------------------------

; ============ file情報 ==========================
; 観測データ
path = '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0313_4.sav'
restore, path

; ============== CO2吸収量の計算を行う ==============================

; ARSで計算した場所の波長帯を取ってきている
CO2=where(wvl gt 1.81 and wvl lt 2.19)
wvl=wvl[CO2]

nwvl=n_elements(wvl)

io=n_elements(jdat(*,1,1))
ip=n_elements(jdat(1,1,*))

obs_spec = reform(jdat(*,0,*))

; CO2の吸収線の中
band=where(wvl gt 1.85 and wvl lt 2.10)

;OMEGAの38番目と48番目の素子が死んでた
nanserch=where(jdat ge 0 and jdat le 0.0001)
jdat(nanserch)=!VALUES.F_NAN

jdat = jdat(*,CO2,*)

for i=0, io-1 do begin
  for j=0, ip-1 do begin
    x = [wvl(0),wvl(1),wvl(2),wvl(23),wvl(24),wvl(25)]
    y = [jdat(i,0,j), jdat(i,1,j), jdat(i,2,j), jdat(i,23,j), jdat(i,24,j), jdat(i,25,j)]
    coef = linfit(x,y)
    cont = coef(0) + coef(1)*wvl
    width = 1- jdat(i,*,j)/cont
  
    obs_spec(i,j) = total(width[*,band],/nan)
  endfor
endfor

save,obs_spec,$
  filename='/Users/nyonn/IDLWorkspace/Default/savfile/Table_SP_obs_calc_orb0313_4.sav'
stop
END