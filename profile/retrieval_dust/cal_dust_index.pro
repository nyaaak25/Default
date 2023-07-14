;------------------------------------------------------------------------------------------------------------------------
;
; ダストのインデックスを作成するためのプログラム
; create by Akira Kazama
; 
; dust_index　 ::2023.7.13 Thu 15:55:00
; 緯度平均を取って、ダスト量を導出するためのプログラム
;
;
;------------------------------------------------------------------------------------------------------------------------


pro cal_dust_index

path = '/work1/LUT/dust/table/output/'

dust_index_file = dblarr(31)

; pathの指定とspecmarsの読み込み
pathfile = '/data2/omega/sav/'
restore, pathfile + 'specmars.sav'
specmars = specmars

for i = 0, 30 do begin

  ; dustのインデックスを作成するためのファイルを読み込む
  if i lt 10 then str_i = strupcase(strmid(i,0,2,/REVERSE_OFFSET))
  if i ge 10 then str_i = strupcase(strmid(i,1,3,/REVERSE_OFFSET))
  print, str_i
  file = path + 'loc1_dust' + str_i + '_albedo.dat'

  ; 毎回毎回初期化をする
  rad1 = dblarr(11)
  wn = dblarr(11)

  ; fileを読み込んで、放射輝度からI/Fを計算する
  openr,lun,file,/get_lun
  for k = 0, 11-1 do begin
    readf,lun,a,b
    wn(k)=a
    rad1(k)=b
  endfor
  free_lun,lun

  wn = reverse(wn)
  wav = 1/wn
  wn = (1/wn)*10000

  ; wn(5): 2.777 um                
  rad1 = reverse(rad1)
  rad1 = (rad1/wav^2)*1e-7  ; 波数から波長換算 (erg = 10-7W, cm^2 = 10-4m^2, μm = 10^4cm)
  rad1 = rad1(5) / specmars(140)

  dust_index_file(i) = rad1
  
endfor

save, dust_index_file, filename = path + 'dust_index_cal.sav'

end