;------------------------------------------------------------------------------------------------------------------------
;
; 熱潮汐波の解析をするためのプログラム
;
; create by Akira Kazama
; 
; tidalwave_analysis　 ::2023.5.13 Sat 07:42:00
; tidalwave_analysis　 ::2023.5.13 Sat 07:42:00
; original seasonal variotion plot tool
;
; tidalwave_analysis_v1　 ::2023.5.30 Tue 17:32:00
; focus on specfic latitude and longitude
; 色々な機能を追加（ダスト/SZA別にplot, MCDとの差をみる）
; TBD! error barを作成、MY3年分を並べてplotする
;
; tidalwave_analysis_v2　 ::2023.6.13 Tue 11:32:00
; MY27-29を混ぜて熱潮汐LTを稼いでplotをする
; 移動平均を取って行う
; retrievalされたデータも含めて解析を行う [2023/6/20追加]
;
;------------------------------------------------------------------------------------------------------------------------

pro tidalwave_analysis_v2

device, decomposed = 0, retain = 2

; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'
path2 = '/work1/LUT/SP/mcddata/'

; restore file
; choose latitude
minlat = -30
maxlat = 0
str_min = strupcase(strmid(minlat,1,4,/REVERSE_OFFSET))
str_max = strupcase(strmid(maxlat,1,4,/REVERSE_OFFSET))

; mcd data save
restore, path_work + 'moving_latitude_' + str_min + '_' + str_max + '_new.sav'
;Ls = ls_ind
localtime = lt_ind
ps = mcd_ps_ind
ps12 = mcd_ps_ind_12

; retireval data save
restore, path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new.sav'
Ls = ls_ind
ret_ps = ps_ind

; create index
ip = n_elements(Ls)
div_season_ps = dblarr(ip)

; Ls color bar
minp6 = 0
maxp6 = 360

; ====== plotの枠組み ========
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

;; plotの枠組みを作成  
window, 0, xs=1500, ys=1000
loadct,39
;plot, findgen(100),xs=1, ys=1, yr=[-105, 105], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='Ls [Deg]', ytitle='mean pressure [%]'
plot, findgen(100),xs=1, ys=1, yr=[-8, 8], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='LT', ytitle='mean pressure [%]'
;cgLOADCT, 39
;cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp6,maxp6],title='Ls [deg]', position=[0.063,0.12,0.067,0.86], TCHARSIZE=10, /vertical 


; ====== 一つの季節できちんと熱絵潮汐がみられているかcheckのためのplot =======
; -------->
;plot, findgen(100),xs=1, ys=1, yr=[290, 1000], xr=[0, 1080], back=255, color=0, /nodata, charsize=3, title='tidal wave detection', xtitle='Ls [Deg]', ytitle='mean pressure [Pa]'
;plot, findgen(100),xs=1, ys=1, yr=[290, 1000], xr=[0, 24], back=255, color=0, /nodata, charsize=3, title='tidal wave detection', xtitle='Ls [Deg]', ytitle='mean pressure [Pa]'
;for i = 0, ip-1 do begin

;    if Ps(i) gt 0 then begin
;        loadct, 39
;        plots, ls(i), ps(i), color=0, psym=8, symsize=2 
;    endif

;endfor
; <--------

; ===== LT 12:00のMCD圧力データで引き算を行う =========
div = ((ret_ps - ps12) /ps12)*100d
;div = ((ps - ps12) /ps12)*100d


for i = 0, ip-1 do begin
;for i = 0, 3694 do begin
    ;if i ge 3694 then Ls(i) = Ls(i) + 360d
    ;if i ge 6222 then Ls(i) = Ls(i) + 360d

    if Ps(i) gt 0d and i ge 0 and i le 3393 and div(i) lt 8 and div(i) gt -8  then begin
    ;if Ps(i) gt 0d and i ge 0 and i le 3393 then begin
        loadct, 39
        ;plots, localtime(i), div(i), color=(LS(i)/(maxp6-minp6))*254, psym=8, symsize=2
        plots, localtime(i), div(i), color=0, psym=8, symsize=2
        ;plots, ls(i), ps(i), color=0, psym=8, symsize=2
        ;plots, ls(i), p_mean, color=254, psym=8, symsize=2
    endif

    if Ps(i) gt 0d and i ge 3694 and i le 6221 and div(i) lt 8 and div(i) gt -8 then begin
    ;if Ps(i) gt 0d and i ge 3694 and i le 6221 then begin
    ;if Ps(i) gt 0d and i ge 3694 and i le 6221 and Ls(i) gt 0d and Ls(i) lt 180d then begin
        ;oadct, 39
        ;plots, localtime(i), div(i), color=(LS(i)/(maxp6-minp6))*254, psym=8, symsize=2
        plots, localtime(i), div(i), color=60, psym=8, symsize=2
        ;plots, ls(i), ps(i), color=0, psym=8, symsize=2
        ;plots, ls(i), p_mean, color=254, psym=8, symsize=2
    endif

    if Ps(i) gt 0d and i ge 6222 and div(i) lt 8 and div(i) gt -8 then begin
    ;if Ps(i) gt 0d and i ge 6222 then begin
        ;loadct, 39
        ;plots, localtime(i), div(i), color=(LS(i)/(maxp6-minp6))*254, psym=8, symsize=2
        plots, localtime(i), div(i), color=254, psym=8, symsize=2
        ;plots, ls(i), ps(i), color=0, psym=8, symsize=2
        ;plots, ls(i), p_mean, color=254, psym=8, symsize=2
    endif

endfor
stop



; ====== 移動平均をとって熱潮汐波の検出を試みる ======
; Ls bin sizeを決める（移動平均の間隔を設定）
bin_size = 5

; ==== 移動平均を取る ======
for i = 0, ip-1 do begin

    ; ±bin_sizeの範囲を探してくる
    ls_ind = where(Ls(*) gt Ls(i) - bin_size and Ls(*) lt Ls(i) + bin_size)
    
    ; その中でのPsを探して、平均値を取る
    p_ind = Ps(Ls_ind)
    p_mean = mean(p_ind, /NaN)

    ; 平均値を引く
    p_div = Ps(i) - p_mean
    div_season_ps(i) = p_div

    maxp6 = 1080d
    minp6 = 0d

    if i ge 3694 then Ls(i) = Ls(i) + 360d
    if i ge 6222 then Ls(i) = Ls(i) + 360d

    ;if Ps(i) gt 0d and Ls(i) gt 90 and Ls(i) lt 180 then begin
    if Ps(i) gt 0d then begin
        loadct, 39
        plots, localtime(i), p_div, color=(LS(i)/(maxp6-minp6))*254, psym=8, symsize=2
        ;plots, localtime(i), Ps(i), color=0, psym=8, symsize=2
        ;plots, ls(i), ps(i), color=0, psym=8, symsize=2
        ;plots, ls(i), p_mean, color=254, psym=8, symsize=2
    endif

endfor

stop



snapshot = TVRD(True=1)
Write_JPEG, path_ql + 'MYmix_lat0-30_moving.jpg', snapshot, True=1, Quality=75


stop

end