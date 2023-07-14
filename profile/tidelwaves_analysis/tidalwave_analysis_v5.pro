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

; tidalwave_analysis_v5　 ::2023.6.28 Wed 15:28:00
; 標準偏差付きのデータを作成する
; 
;
;------------------------------------------------------------------------------------------------------------------------

pro tidalwave_analysis_v5

device, decomposed = 0, retain = 2

; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'
path2 = '/work1/LUT/SP/mcddata/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

; restore file
; choose latitude
minlat = 0
maxlat = 30
str_min = strupcase(strmid(minlat,1,4,/REVERSE_OFFSET))
str_max = strupcase(strmid(maxlat,1,4,/REVERSE_OFFSET))

; mcd data save
restore, path_work + 'moving_latitude_' + str_min + '_' + str_max + '.sav'
;Ls = ls_ind
localtime = lt_ind
ps = mcd_ps_ind
ps12 = mcd_ps_ind_12

; === retireval data save ===
; _1new is 1 sigma remove
; _1new_std is indluding std information
; _1new_DQ is  data quality less than 3 is not used.
; _1new_std_DQ is  data quality less than 3/2 is not used.
; ===========================

;restore, path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new_DQ.sav'
;restore, path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new.sav'
restore, path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new_std_DQ.sav'

Ls = ls_ind
ret_ps = ps_ind
ret_std = ps_std
albedo = albedo_ind
dust = dust_ind

; create index
ip = n_elements(Ls)

; ====== plotの枠組み ========
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

;; plotの枠組みを作成  
window, 0, xs=1500, ys=1000
loadct,39
plot, findgen(100),xs=1, ys=1, yr=[-105, 105], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='LT', ytitle='mean pressure [%]'
;plot, findgen(100),xs=1, ys=1, yr=[-20, 20], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='LT', ytitle='mean pressure [%]'

div = ((ret_ps - ps12) /ps12)*100d

; 平均値内のそもそもの標準偏差でどんだけのばらつきがあるかをみてみる
min_div = ret_ps - ret_std
max_div = ret_ps + ret_std
min_ran = ((min_div - ps12) /ps12)*100d
max_ran = ((max_div - ps12) /ps12)*100d

; dust/albedo conditionによって色を変える
minc1 = min(albedo)
maxc1 = max(albedo)

;minc1 = min(dust)
;maxc1 = max(dust)

; Color bar をつける
cgLOADCT, 39
cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minc1,maxc1],title='Albedo', position=[0.063,0.12,0.067,0.86], TCHARSIZE=10, /vertical 


; ===== LT 12:00のMCD圧力データで引き算を行う =========
for i = 0, ip-1 do begin
    if i ge 3694 then Ls(i) = Ls(i) + 360d
    if i ge 6222 then Ls(i) = Ls(i) + 360d

    ; 標準偏差を書く、±5%の絶対値の誤差が乗る可能性があるので、それがどれだけ聞くかをみてみる
    ; 一つの観測結果に±5%の変動が乗っていると思って、error barを作成
    ; ------->
    ;min_ret = ret_ps(i) - ret_ps(i) * 0.05d
    ;max_ret = ret_ps(i) + ret_ps(i) * 0.05d
    ;min_ran = ((min_ret - ps12(i)) /ps12(i))*100d
    ;max_ran = ((max_ret - ps12(i)) /ps12(i))*100d

    ;ran = [min_ran, max_ran]
    ; <-------
    
    ;if Ps(i) gt 0d and i ge 0 and i le 3393 and Ls(i) gt 0d and Ls(i) lt 180d then begin
    if Ps(i) gt 0d and i ge 0 and i le 3393 and Ls(i) gt 90d and Ls(i) lt 180d then begin
    ;if Ps(i) gt 0d and i ge 0 and i le 3393 and Ls(i) gt 90d and Ls(i) lt 180d and ret_std(i) lt 15 then begin
        ran = [min_ran(i), max_ran(i)]
        loadct, 39
        plots, localtime(i), ran, color=(albedo(i)/(maxc1-minc1))*254, THICK=3
        plots, localtime(i), div(i), color=(albedo(i)/(maxc1-minc1))*254, psym=8, symsize=2
        ;plots, localtime(i), div(i), color=0, psym=8, symsize=2
    endif

    ;if Ps(i) gt 0d and i ge 3694 and i le 6221 and Ls(i) gt 360d and Ls(i) lt 540d then begin
    if Ps(i) gt 0d and i ge 3694 and i le 6221 and Ls(i) gt 450d and Ls(i) lt 540d then begin
    ;if Ps(i) gt 0d and i ge 3694 and i le 6221 and Ls(i) gt 450d and Ls(i) lt 540d and ret_std(i) lt 15 then begin
        ran = [min_ran(i), max_ran(i)]
        loadct, 39
        plots, localtime(i), ran, color=(albedo(i)/(maxc1-minc1))*254, THICK=3
        plots, localtime(i), div(i), color=(albedo(i)/(maxc1-minc1))*254, psym=8, symsize=2
        ;plots, localtime(i), div(i), color=60, psym=8, symsize=2
    endif

    ;if Ps(i) gt 0d and i ge 6222 and Ls(i) gt 720d and Ls(i) lt 900d then begin
    if Ps(i) gt 0d and i ge 6222 and Ls(i) gt 810d and Ls(i) lt 900d then begin
    ;if Ps(i) gt 0d and i ge 6222 and Ls(i) gt 810d and Ls(i) lt 900d and ret_std(i) lt 15 then begin
        ran = [min_ran(i), max_ran(i)]
        loadct, 39
        plots, localtime(i), ran, color=(albedo(i)/(maxc1-minc1))*254, THICK=3
        plots, localtime(i), div(i), color=(albedo(i)/(maxc1-minc1))*254, psym=8, symsize=2
        ;plots, localtime(i), div(i), color=254, psym=8, symsize=2
    endif

endfor
stop

; image save
snapshot = TVRD(True=1)
Write_JPEG, path_ql + 'focuson_Ls90-180_5err_DQ_sourth.jpg', snapshot, True=1, Quality=75


stop

end