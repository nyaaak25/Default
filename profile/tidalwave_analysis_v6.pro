;------------------------------------------------------------------------------------------------------------------------
;
; 熱潮汐波の解析をするためのプログラム
;
; create by Akira Kazama
; 
; tidalwave_analysis　 ::2023.5.13 Sat 07:42:00
; original seasonal variotion plot tool
;
; tidalwave_analysis_v1　 ::2023.5.30 Tue 17:32:00
; focus on specfic latitude and longitude
; 色々な機能を追加（ダスト/SZA別にplot, MCDとの差をみる）
; TBD! error barを作成、MY3年分を並べてplotする

; tidalwave_analysis_v2　 ::2023.6.13 Tue 11:32:00
; MY27-29を混ぜて熱潮汐LTを稼いでplotをする
; retrievalされたデータも含めて解析を行う [2023/6/20追加]
;
; tidalwave_analysis_v3　 ::2023.6.13 Tue 16:59:00
; 色々なメモ代わりに使えるプログラム

; tidalwave_analysis_v4　 ::2023.6.28 Wed 15:21:00
; Localtime方向の平均値+標準偏差を取るようなplotを作成する
;
; tidalwave_analysis_v5　 ::2023.6.28 Wed 15:28:00
; 標準偏差付きのデータを作成する

; tidalwave_analysis_v5_1　 ::2023.7.4 Tue 15:07:00
; 熱潮汐波の検出をする
; 平均方法：各スペクトルで中央値を探す
; データクオリティ2/3以上のものを使用する, 1σに入らないデータは除く
; データセレクション：高度±10 km, リムデータを除く, Albedo: 0.2-0.5
;
; tidalwave_analysis_v5_2　 ::2023.7.5 Wed 9:32:00
; 熱潮汐波の検出をする
; function化 + Ls/ local time indexを付加
; 平均方法：各スペクトルで中央値を探す
; データクオリティ2/3以上のものを使用する, 1σに入らないデータは除く
; データセレクション：高度±10 km, リムデータを除く
; !TBD! albedo, water ice, dust index use

; tidalwave_analysis_v6　 ::2023.7.5 Wed 16:10:00
; 熱潮汐波の検出をする
; v5_2で作成したファイルを読み込んで、画像を作るプログラム
;
;------------------------------------------------------------------------------------------------------------------------


pro tidalwave_analysis_v6

device, decomposed = 0, retain = 2

; ================ path + Ls/ localtime index serach ================ 
;path_data = '/work1/LUT/SP/data_analysis/'
path_data = '/work1/LUT/SP/data_analysis_index/'


; ====== plotの枠組み ========
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

;; plotの枠組みを作成  
window, 0, xs=1500, ys=1000
loadct,39
;plot, findgen(100),xs=1, ys=1, yr=[-30, 30], xr=[0, 24], back=255, color=0, /nodata, charsize=3, title='tidel wave detection', xtitle='Local time', ytitle='mean pressure [%]'
plot, findgen(100),xs=1, ys=1, yr=[-15, 15], xr=[0, 24], back=255, color=0, /nodata, charsize=3, title='tidel wave detection', xtitle='Local time', ytitle='mean pressure [%]'


; +++++++++ Ls / Local time / latitude 指定 ++++++++++++++
; Lsの定義
min_Ls = 0
max_Ls = 90
str_min_ls = strupcase(strmid(min_ls,1,4,/REVERSE_OFFSET))
str_max_ls = strupcase(strmid(max_ls,1,4,/REVERSE_OFFSET))

; 2つの季節を混ぜてみる
min_Ls2 = 90
max_Ls2 = 180
str_min_ls2 = strupcase(strmid(min_ls2,1,4,/REVERSE_OFFSET))
str_max_ls2 = strupcase(strmid(max_ls2,1,4,/REVERSE_OFFSET))

; local timeの定義
dev_local = 1d  ; Local time 何分刻みにしたいかをここで指定

for i = 0, 24d, dev_local do begin
    str_i = strupcase(strmid(i,1,4,/REVERSE_OFFSET))

    ; 緯度を定義
    min_latitude = 30
    max_latitude = 70
    str_min = strupcase(strmid(min_latitude,1,4,/REVERSE_OFFSET))
    str_max = strupcase(strmid(max_latitude,1,4,/REVERSE_OFFSET))

    ; fileができていないところのskipをする
    ; _lon: longitudeを100度から180度に限定しているファイル

    files = file_search(path_data + 'local_' + '*.sav')
    file_find = where(files eq path_data + 'local_' + str_i + '_lat' + str_min + '_' + str_max + '_' + 'ls_' + str_min_ls + '_' + str_max_ls + '_lon.sav')
    if file_find(0) eq -1 then goto, skip2
    if file_find(0) ge 0 then restore, path_data + 'local_' + str_i + '_lat' + str_min + '_' + str_max + '_' + 'ls_' + str_min_ls + '_' + str_max_ls + '_lon.sav'

    ; 1つ目の季節をplotしてみる
    ret_data1 = spectrum(*,0)
    mcd_data1 = spectrum(*,1)

    skip2:

    ; 2つの季節を混ぜて考えてみる
    file_find2 = where(files eq path_data + 'local_' + str_i + '_lat' + str_min + '_' + str_max + '_' + 'ls_' + str_min_ls2 + '_' + str_max_ls2 + '_lon.sav')
    if file_find(0) eq -1 and file_find2(0) eq -1 then goto, skip1
    if file_find2(0) eq -1 then goto, skip3
    if file_find2(0) ge 0 then restore, path_data + 'local_' + str_i + '_lat' + str_min + '_' + str_max + '_' + 'ls_' + str_min_ls2 + '_' + str_max_ls2 + '_lon.sav'

    ; ２つ目の季節をplotしてみる
    ret_data2 = spectrum(*,0)
    mcd_data2 = spectrum(*,1)

    skip3:

    if file_find2(0) eq -1 and file_find(0) ge 0 then begin
        ret_data = ret_data1
        mcd_data = mcd_data1
    endif

    if file_find(0) eq -1 and file_find2(0) ge 0 then begin
        ret_data = ret_data2
        mcd_data = mcd_data2
    endif

    if file_find(0) ge 0 and file_find2(0) ge 0 then begin
        ret_data = [ret_data1, ret_data2]
        mcd_data = [mcd_data1, mcd_data2]
    endif


    ; 熱潮汐の検出のための解析をここで行う
    ; OMEGAのリトリーバルデータ
    medi_ret = median(ret_data)
    std_ret = stddev(ret_data,/nan)

    ret_ran = [medi_ret - std_ret, medi_ret + std_ret]

    ; mcd データ
    medi_mcd = median(mcd_data)
    std_mcd = stddev(mcd_data,/nan)

    mcd_ran = [medi_mcd - std_mcd, medi_mcd + std_mcd]

    plots, i, medi_ret, color=0, psym=8, symsize=2
    plots, i, ret_ran, color=0, THICK=3
    
    plots, i, medi_mcd, color=254, psym=8, symsize=2
    ;plots, i, mcd_data, color=110, psym=8, symsize=0.5
    plots, i, mcd_ran, color=254, THICK=3

    ind_hist = where(abs(ret_data) gt 0)
    if ind_hist(0) eq -1 then goto, skip1

    ;window, 1
    ;data = ret_data
    ;n = 0.1d
    ;loadct, 39 
    ;hist = histogram(data, binsize=n, /nan) ;data：ヒストグラムを作りたいデータarray、binsize：どのくらいの細かさか
    ;bin = (findgen(n_elements(data))*n) + min(data,/nan) 
    ;plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr= [0,n_elements(data)/50], xr=[min(data,/nan)-1d, + max(data,/nan) + 1d], title='Histogram'

    ;print, n_elements(data)
    ;print, min(data,/nan)
    ;print, max(data,/nan)
    ;print, mean(data,/nan)
    ;print, medi_ret
    ;stop

    skip1:
  
endfor


stop

end