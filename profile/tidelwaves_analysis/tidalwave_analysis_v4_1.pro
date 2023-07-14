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

; tidalwave_analysis_v3　 ::2023.6.13 Tue 16:59:00
; 他index別にみる
;
; tidalwave_analysis_v4　 ::2023.6.28 Wed 15:21:00
; Localtime方向の平均値+標準偏差を取るようなplotを作成する
;
; tidalwave_analysis_v4_1　 ::2023.7.3 Mon 11:45:00
; Localtime方向の平均値+標準偏差を取るようなplotを作成する
; ひとつひとつのスペクトルで作成してみる
; 
;
;------------------------------------------------------------------------------------------------------------------------

pro tidalwave_analysis_v4_1

device, decomposed = 0, retain = 2

; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'

; restore file
; choose latitude
minlat = 0
maxlat = 30
str_min = strupcase(strmid(minlat,1,4,/REVERSE_OFFSET))
str_max = strupcase(strmid(maxlat,1,4,/REVERSE_OFFSET))

; file search for EW file
files=file_search(path_work+'*.sav',count=count)

; mcd data save
restore, path_work + 'moving_latitude_' + str_min + '_' + str_max + '.sav'
localtime = lt_ind
ps = mcd_ps_ind
ps12 = mcd_ps_ind_12


; === retireval data save ===
; _1new is 1 sigma remove
; _1new_DQ is  data quality less than 3 is not used.
; ===========================

;restore, path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new_DQ.sav'
restore, path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new_std_DQ.sav'
Ls = ls_ind
ret_ps = ps_ind

; ====== plotの枠組み ========
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

;; plotの枠組みを作成  
window, 0, xs=1500, ys=1000
loadct,39
plot, findgen(100),xs=1, ys=1, yr=[400, 1300], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='Local time', ytitle='pressure [Pa]'
;plot, findgen(100),xs=1, ys=1, yr=[-30, 30], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='Local time', ytitle='mean pressure [%]'
;plot, findgen(100),xs=1, ys=1, yr=[-10, 10], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='Local time', ytitle='mean pressure [%]'

div = ((ret_ps - ps12) /ps12)*100d
div_mcd = ((ps - ps12) /ps12)*100d

; ===== LT方向に積分を行ってみる =========

;!! new version !!
; 使うLsを限定する
ind_ls = where(Ls ge 180d and Ls le 360d)
ind_ls2 = where(Ls ge 0d and Ls le 90d)
div(ind_ls) = -0d/0d
div_mcd(ind_ls) = -0d/0d
div(ind_ls2) = -0d/0d

; どのくらいの分解能で平均値を取るかをここで指定
dev = 1d

; Loop start
;for i= 0d, 24d, dev do begin
for i= 16d, 16d, dev do begin
    ind = where(localtime gt i and localtime lt i + dev)
    ind_2 = where(abs(div(ind)) gt 0, count)

    sum_cont = 0d
    
    for loop = 0, count-1 do begin
        restore, files(ind(ind_2(loop)))
        fileorbit = strupcase(strmid(files(ind(ind_2(loop))),12,9,/REVERSE_OFFSET))

        minlon = min(longi)
        maxlon = max(longi)
        ind_xyz = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat and lati le maxlat, xind=xind, yind=yind)

        ip = n_elements(LATI(*,0))
        io = n_elements(LATI(0,*))

        ; p1 を定義する
        p1 = p1_slev

        ; bad pixelにNanを格納
        for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN

        ; bad dataを弾くように設定
        pm = 1  ; 何σに入るかをここで定義する！
        std = stddev(p1,/Nan)
        medi = median(p1)

        min_range = medi - (pm * std)
        max_range = medi + (pm * std)

        new_p1 = dblarr(ip,io)
        ind2 = where_xyz(p1 ge min_range and p1 le max_range, xind=xind, yind=yind)
        new_p1(ind2) = p1(ind2)
        for i = 0, ip-1 do for j = 0, io-1 do if new_p1(i,j) eq 0 then new_p1(i,j) = !VALUES.F_NAN

        local = localtime(ind(ind_2(loop)))
        for i = 0, ip-1 do for j = 0, io-1 do if new_p1(i,j) gt 0 then plots, local, new_p1(i,j), color=0, psym=6, symsize=0.5

        ; その中で平均を取る
        various_p1 = mean(new_p1(ind_xyz), /Nan)
        print, 'mean: ', various_p1

        cont_ind = where(new_p1(ind_xyz) gt 0)
        cont_abs = n_elements(cont_ind)
        print, 'sample; ', cont_abs

        sum_cont = sum_cont + cont_abs
        print, 'sum: ', sum_cont
        print, 'ORB:', fileorbit

        ;window, 1
        ;data = new_p1(ind_xyz)
        ;n = 1d
        ;loadct, 39 
        ;hist = histogram(data, binsize=n, /nan) ;data：ヒストグラムを作りたいデータarray、binsize：どのくらいの細かさか
        ;bin = (findgen(n_elements(data))*n) + min(data,/nan) 
        ;plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr= [0,10], xr=[min(data,/nan)-10d, + max(data,/nan) + 10d], title='Histogram'

    endfor
endfor


stop

end