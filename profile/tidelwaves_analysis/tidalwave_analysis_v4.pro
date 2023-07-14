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
; TBD! MY3年分を並べてplotする

; tidalwave_analysis_v2　 ::2023.6.13 Tue 11:32:00
; MY27-29を混ぜて熱潮汐LTを稼いでplotをする
; retrievalされたデータも含めて解析を行う [2023/6/20追加]

; tidalwave_analysis_v3　 ::2023.6.13 Tue 16:59:00
; 他index別にみる
;
; tidalwave_analysis_v4　 ::2023.6.28 Wed 15:21:00
; Localtime方向の平均値+標準偏差を取るようなplotを作成する
; 重み付け平均をするための工夫をおこなう
; 
;
;------------------------------------------------------------------------------------------------------------------------


;+
; PURPOSE:
;  Calculates the weighted mean of a set of measurements.
;
; CATEGORY:
;  Statistics
;
; CALLING SEQUENCE:
;  result = wmean(val, dval, [/nan, error = error])
; 
; INPUTS:
;  val : A vector of data values
;  dval: The error on each data value
;
; KEYWORD PARAMETERS:
;  nan: If set, treat NANs as missing data
;  error: Set to a named variable to hold the error on the weighted
;         mean
; 
; OUTPUTS:
;  The weighted mean
;
; MODIFICATION HISTORY
;  April 2009 Written by Chris Beaumont
;  April 23 2009 fixed integer truncation bug
;- 

function wmean, val, dval, nan = nan, error = error
compile_opt idl2
on_error, 2

;- check inputs
if n_params() ne 2 then begin
   print, 'wmean calling sequence:'
   print,' result = wmean(val, dval, [/nan, error = error])'
   return, !values.f_nan
endif
sz = n_elements(val)
if n_elements(dval) ne sz then $
   message, 'val and dval must contain the same number of elements'


;- the calculation
;- require at least one finite value
if keyword_set(nan) then begin
   hit = where(finite(val) and finite(dval), ct)
   if ct eq 0 then begin
      error = !values.f_nan
      return, !values.f_nan
   endif
endif

dw = total(1D / dval^2, nan = keyword_set(nan))
sum = total(1D * val / dval^2, nan = keyword_set(nan)) / dw
error = sqrt(1D / dw)

return, sum

end




pro tidalwave_analysis_v4


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
; _1new_std_DQ is ata quality less than 3/2 is not used.
; ===========================

restore, path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new_std_DQ.sav'
;restore, path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new.sav'
Ls = ls_ind
ret_ps = ps_ind
ret_std = ps_std

; ====== plotの枠組み ========
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

;; plotの枠組みを作成  
window, 0, xs=1500, ys=1000
loadct,39
plot, findgen(100),xs=1, ys=1, yr=[-40, 40], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='Local time', ytitle='mean pressure [%]'
;plot, findgen(100),xs=1, ys=1, yr=[-10, 10], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='Local time', ytitle='mean pressure [%]'

div = ((ret_ps - ps12) /ps12)*100d
div_mcd = ((ps - ps12) /ps12)*100d

lt12_std = (ret_std / ps12)*100d

min_div = ret_ps - ret_std
max_div = ret_ps + ret_std
min_ran = ((min_div - ps12) /ps12)*100d
max_ran = ((max_div - ps12) /ps12)*100d

; ===== LT方向に積分を行ってみる =========

;!! new version !!
; 使うLsを限定する
ind_ls = where(Ls ge 180d and Ls le 360d)
;ind_ls2 = where(Ls ge 0d and Ls le 90d)
div(ind_ls) = -0d/0d
div_mcd(ind_ls) = -0d/0d
;div(ind_ls2) = -0d/0d

; どのくらいの分解能で平均値を取るかをここで指定
dev = 0.5d

; Loop start
for i= 16d, 24d, dev do begin


    ind = where(localtime gt i and localtime lt i + dev)
    ave1 = mean(div(ind), /nan)
    ave_mcd = mean(div_mcd(ind),/nan)

    ; err 分散を記入
    std = stddev(div(ind),/nan)
    std_mcd = stddev(div_mcd(ind),/nan)

    ; 重み付け平均
   ave = wmean(div(ind), lt12_std(ind), /nan)

   std_ran = [ave-std, ave+std]
   std_ran_mcd = [ave_mcd-std_mcd, ave_mcd+std_mcd]

   plots, i , div(ind), color= 100, psym=8, symsize=1
    plots, i, ave_mcd, color=60, psym=8, symsize=2
    plots, i, ave, color=0, psym=8, symsize=2
    ;plots, i, ave1, color=254, psym=8, symsize=2
    ;plots, i, std_ran, color=0, THICK=3

    ;plots, i, std_ran_mcd, color=254, THICK=3

    ind2 = where(abs(div(ind)) gt 0,cont)
    cont_abs = cont

    sample = strmid(cont_abs,1,7)

    print, 'sample number: ', i, cont
    
    ; ==== Histogram ====
    ;if cont_abs gt 0 then begin 
    ;    window, 1
    ;    data = div(ind)
    ;    n = 1d
    ;    loadct, 39 
    ;    hist = histogram(data, binsize=n) ;data：ヒストグラムを作りたいデータarray、binsize：どのくらいの細かさか
    ;    bin = (findgen(n_elements(data))*n) + min(data) 
    ;    plot, bin, hist, psym=10, back=255, color=0, xs=1, ys=1, yr= [0,cont_abs/3d], xr=[min(data)-10d, + max(data) + 10d], title='Histogram'
    ;endif
    ;stop
    stop

endfor

; !! old verdion !!
;ind7 = where(localtime gt 6d and localtime lt 7d)
;ind8 = where(localtime gt 7d and localtime lt 8d)
;ind9 = where(localtime gt 8d and localtime lt 9d)
;ind10 = where(localtime gt 9d and localtime lt 10d)
;ind11 = where(localtime gt 10d and localtime lt 11d)
;ind12 = where(localtime gt 11d and localtime lt 12d)
;ind13 = where(localtime gt 12d and localtime lt 13d)
;ind14 = where(localtime gt 13d and localtime lt 14d)
;ind15 = where(localtime gt 14d and localtime lt 15d)
;ind16 = where(localtime gt 15d and localtime lt 16d)
;ind17 = where(localtime gt 16d and localtime lt 17d)
;
;mean7 = mean(div(ind7), /nan)
;mean8 = mean(div(ind8), /nan)
;mean9 = mean(div(ind9), /nan)
;mean10 = mean(div(ind10), /nan)
;mean11 = mean(div(ind11), /nan)
;mean12 = mean(div(ind12), /nan)
;mean13 = mean(div(ind13), /nan)
;mean14 = mean(div(ind14), /nan)
;mean15 = mean(div(ind15), /nan)
;mean16 = mean(div(ind16), /nan)
;mean17 = mean(div(ind17), /nan)

;new_loc = [7d, 8d, 9d, 10d, 11d, 12d, 13d, 14d, 15d, 16d, 17d]
;new_div = [mean7, mean8, mean9, mean10, mean11, mean12, mean13, mean14, mean15, mean16, mean17]
;plots, new_loc, new_div, color=0, psym=8, symsize=2


stop

end