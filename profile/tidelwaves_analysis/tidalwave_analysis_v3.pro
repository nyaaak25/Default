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
; 色々なメモ代わりに使えるプログラム
; 
;
;------------------------------------------------------------------------------------------------------------------------

pro tidalwave_analysis_v3

device, decomposed = 0, retain = 2

; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'
path2 = '/work1/LUT/SP/mcddata/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB0181_0.sav') ; MY27 Ls 0
loop_b = fix(loop_b(0))

; restore file
; choose latitude
minlat = -30
maxlat = 0
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
; _1new_DQ is  data quality less than 3 is not used.
; ===========================

;restore, path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new_DQ.sav'
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
plot, findgen(100),xs=1, ys=1, yr=[-50, 50], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='Ls [Deg]', ytitle='mean pressure [%]'
;plot, findgen(100),xs=1, ys=1, yr=[-15, 15], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='Ls [Deg]', ytitle='mean pressure [%]'
;plot, findgen(100),xs=1, ys=1, yr=[-8, 8], xr=[0, 24], back=255, color=0, /nodata, charsize=3, position=[0.18,0.1,0.98,0.9], title='tidel wave detection', xtitle='LT', ytitle='mean pressure [%]'
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

div = ((ret_ps - ps12) /ps12)*100d

; ===== LT 12:00のMCD圧力データで引き算を行う =========
for i = 0, ip-1 do begin
    if i ge 3694 then Ls(i) = Ls(i) + 360d
    if i ge 6222 then Ls(i) = Ls(i) + 360d

    ; 標準偏差を書く、±5%の絶対値の誤差が乗る可能性があるので、それがどれだけ聞くかをみてみる
    min_ret = ret_ps(i) - ret_ps(i) * 0.05d
    max_ret = ret_ps(i) + ret_ps(i) * 0.05d
    min_ran = ((min_ret - ps12(i)) /ps12(i))*100d
    max_ran = ((max_ret - ps12(i)) /ps12(i))*100d

    ran = [min_ran, max_ran]

    ;if Ps(i) gt 0d and i ge 0 and i le 3393 and Ls(i) gt 0d and Ls(i) lt 180d then begin
    if Ps(i) gt 0d and i ge 0 and i le 3393 and Ls(i) gt 90d and Ls(i) lt 180d then begin
        loadct, 39
        ;plots, localtime(i), div(i), color=(LS(i)/(maxp6-minp6))*254, psym=8, symsize=2
        plots, localtime(i), div(i), color=0, psym=8, symsize=2
        ;plots, localtime(i), ran, color=0, THICK=3
        ;plots, ls(i), ps(i), color=0, psym=8, symsize=2
        ;plots, ls(i), p_mean, color=254, psym=8, symsize=2
    endif

    ;if Ps(i) gt 0d and i ge 3694 and i le 6221 and Ls(i) gt 360d and Ls(i) lt 540d then begin
    if Ps(i) gt 0d and i ge 3694 and i le 6221 and Ls(i) gt 450d and Ls(i) lt 540d then begin
        loadct, 39
        ;plots, localtime(i), div(i), color=(LS(i)/(maxp6-minp6))*254, psym=8, symsize=2
        plots, localtime(i), div(i), color=60, psym=8, symsize=2
        ;plots, localtime(i), ran, color=60, THICK=3
        ;plots, ls(i), ps(i), color=0, psym=8, symsize=2
        ;plots, ls(i), p_mean, color=254, psym=8, symsize=2
    endif

    ;if Ps(i) gt 0d and i ge 6222 and Ls(i) gt 720d and Ls(i) lt 900d then begin
    if Ps(i) gt 0d and i ge 6222 and Ls(i) gt 810d and Ls(i) lt 900d then begin
        loadct, 39
        ;plots, localtime(i), div(i), color=(LS(i)/(maxp6-minp6))*254, psym=8, symsize=2
        plots, localtime(i), div(i), color=254, psym=8, symsize=2
        ;plots, localtime(i), ran, color=254, THICK=3
        ;plots, ls(i), ps(i), color=0, psym=8, symsize=2
        ;plots, ls(i), p_mean, color=254, psym=8, symsize=2
    endif

endfor
stop




; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

; MY27
loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB0181_0.sav') ; MY27 Ls 0
loop_b = fix(loop_b(0))
;loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB2595_0.sav') ; MY27 Ls 360
;loop_a = fix(loop_a(0))

; MY 28 Ls plot
;loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB2607_0.sav') ; MY28 Ls 0
;loop_b = fix(loop_b(0))
;loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB5055_4.sav') ; MY28 Ls 360
;loop_a = fix(loop_a(0))

; MY 29 Ls plot
;loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB5062_0.sav') ; MY29 Ls 0
;loop_b = fix(loop_b(0))
loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB7454_4.sav') ; MY29 Ls 360
loop_a = fix(loop_a(0))

; create index
Ls_ind = dblarr(loop_a - loop_b)
mean_p1 = dblarr(loop_a - loop_b)

; ====== 任意の季節を選択 ======
; Ls
minls = 0
maxls = 360

; ===== curiosity data ======
; ----->
g = 3.72 ;m/s2
Gconst = 6.67430d-11
MMars = 6.42d23
RMars = 3.4d6
g = -Gconst*MMars/(-RMars*RMars)
R = 192d ;J/K/kg

cs_slev = dblarr(25)
curiosity_ls = dblarr(25)
curiosity_ps = dblarr(25)
curiosity_ta = dblarr(25)

file_curiosity = '/work1/LUT/Soft/curiosity.dat'
openr,lun,file_curiosity,/get_lun
for i = 0, 25-1 do begin
    readf,lun,a,b
    curiosity_ls(i)=a
    curiosity_ps(i)=b
endfor

file_curiosity_ta = '/work1/LUT/Soft/curiosity_ta.dat'
openr,lun,file_curiosity_ta,/get_lun
for i = 0, 25-1 do begin
    readf,lun,a
    curiosity_ta(i)=a
endfor

cs_altitude = -4488.0483d

for i = 0, 25-1 do cs_slev(i) = curiosity_ps(i) * exp( (cs_altitude) / (R*curiosity_ta(i)/(-Gconst*MMars/(-1d*(RMars+cs_altitude)*(RMars+cs_altitude) ))))

; <-----


; ============ plot ここから =================

; プロットの枠組み
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

;; plotの枠組みを作成
;window, 0, xs=1500, ys=1000
window, 0, xs=750, ys=500
loadct,39
plot, findgen(360),xs=1, ys=1, yr=[290, 1000], xr=[minls, maxls], back=255, color=0, /nodata, charsize=3, title='MY29 ret - mcd: Ls ' +strmid(minls,5,6)+ ' to ' +strmid(maxls,5,6), xtitle='Ls [Deg]', ytitle='mean pressure [Pa]'

; ============= file loopを始める ====================
for loop = loop_b, loop_a do begin
    
    ; fileのrestore
    file = files(loop)
    restore, file

    ; 軌道番号を文字として保存
    fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))

    ; n_indを格納
    ip = n_elements(LATI(*,0))
    io = n_elements(LATI(0,*))

    ; QLが作成されなかったfileはskipするように選択
    skip_ind = where(pressure eq 0)
    n_skip = n_elements(skip_ind)
    skip_ratio = n_skip *100d / (ip * io)
    if skip_ratio gt 97d  then goto, skip1

    ; longitude, latitude
    lat1 = lati
    lon1 = longi

    ; ls
    ls = SOLAR_LONGITUDE

    ; retrieval surface pressure [Done altitude correction]
    p2 = p1_slev

    ; bad pixelにNanを格納
    for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 0 then p2(i,j) = !VALUES.F_NAN

    ; MCD surface pressure
    p1_mcd = mcdpressure 

    ; bad pixelにNanを格納
    for i = 0, ip-1 do for j = 0, io-1 do if p1_mcd(i,j) eq 0 then p1_mcd(i,j) = !VALUES.F_NAN

    ; mcd pressure altutude corecction
    ip = n_elements(lat1(*,0))
    io = n_elements(lat1(0,*))
    p1 = dblarr(ip,io)
    p5 = altitude
    TAmap = TAMAP

    for i = 0, ip-1 do for j = 0, io-1 do if p1_mcd(i,j) gt 0 then p1(i,j) = p1_mcd(i,j) * exp( (p5(i,j)) / (R*TAmap(i,j)/(-Gconst*MMars/(-1d*(RMars+p5(i,J))*(RMars+p5(i,J)) ))))


; =================== choose define value ========================
    ; ====== 任意のSZAを選択 ======
    ; restore dust information
    SZA = SZA_all
    for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then SZA(i,j) = !VALUES.F_NAN

    maxSZA = max(SZA)
    define_minSZA = 55
    define_maxSZA = 60

    ; ====== 任意のalbedoを選択 ======
    ; restore albedo information
    albedo = inputalbedo
    for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then albedo(i,j) = !VALUES.F_NAN

    minalbedo = min(albedo)
    if minalbedo lt  0.1 then print, fileorbit
    if minalbedo lt 0.1 then stop
    maxalbedo = max(albedo)
    define_minalbedo = 0.2
    define_maxalbedo = 0.3

    ; ====== 任意のdust選択 ======
    ; restore dust information
    dust = dustmap

    for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then dust(i,j) = !VALUES.F_NAN
    for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN

    maxdust = max(dust)
    define_mindust = 1.3
    define_maxdust = 1.4

    ; ====== 任意の緯度選択 ======
    ; lat: -70 to -30、-30 to 0、0 to 30、30 to 70
    minlon = min(lon1)
    maxlon = max(lon1)
    minlat = -70
    maxlat = -30

    ; 上記の任意の緯度経度内の圧力を持ってくる
    ind = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat and lati le maxlat, xind=xind, yind=yind)
    if ind(0) eq -1 then goto, skip1

    ; その中で平均を取る
    various_p1 = mean(p1(ind), /Nan)
    various_p2 = mean(p2(ind), /Nan)
    div_ps = various_p1 - various_p2

    ; (TBD) error barをつける（±2.5%の範囲）
    ;yerr = various_p2 * 0.025


; ==================== 実際の plotを作成 ==============================    
    ; plot

    ;; とりあえず全部をなにも考えずにplot用
    if various_p2 gt 0  and ls gt minls and ls lt maxls then begin

    ; ダスト別/SZA別にplotするとき用の指定
    ;if various_p2 gt 0 and ls gt minls and ls lt maxls and maxdust gt define_mindust and maxdust lt define_maxdust then begin
    ;if various_p2 gt 0 and ls gt minls and ls lt maxls and maxSZA gt define_minSZA and maxSZA lt define_maxSZA then begin
    ;if various_p2 gt 0 and ls gt minls and ls lt maxls and maxalbedo gt define_minalbedo and maxalbedo lt define_maxalbedo then begin

    ; MCDとの差が大きいところ（下） or　小さいところを見るとき（上）に使うもの
    ;if various_p2 gt 0 and various_p2 lt 1000 and ls gt minls and ls lt maxls and div_ps gt -50 and div_ps lt 50 then begin
    ;if various_p2 gt 0 and various_p2 lt 1000 and ls gt minls and ls lt maxls and div_ps gt 200 then begin

        loadct, 39
        plots, ls, various_p1, color=0, psym=8, symsize=2
        plots, ls, various_p2, color=254, psym=8, symsize=2

        ; ======== print option ============
        ;print, 'pressure: ', various_p1
        ;print, 'div pressure: ', div_ps
        ;print, 'ls: ', ls
        ;print, 'local time: ', mean(local_time)
        ;print, 'dust', maxdust
        print, 'ORB number: ', fileorbit
        ;print, 'orbloop: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1


    endif


    skip1:

endfor
stop

; ===== curiosity seosonal value plot ======
loadct,39
plots, curiosity_ls, cs_slev, color=205, thick = 6
plots, curiosity_ls, cs_slev, color=205, psym = 8, symsize=2

; image save
snapshot = TVRD(True=1)
Write_JPEG, path_ql + 'focuson_Ls90-180_5err_DQ_sourth.jpg', snapshot, True=1, Quality=75


stop

end