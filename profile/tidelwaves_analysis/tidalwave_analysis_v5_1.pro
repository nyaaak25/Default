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
;
;------------------------------------------------------------------------------------------------------------------------

function data_selection, index

; TBD
aa = 1

return, aa

end



pro tidalwave_analysis_v5_1

; function化をして、Lsとlocal timeをいれたら、そのLocal time/Lsにあたる軌道だけをもってくる
; inputはloop index, outputは何%変動かどうか
; さらにそこから↓のしきい値を決めて、データセレクションをするfunctionへとぶちこむ

; local_time / Ls indexをまずは作成する


device, decomposed = 0, retain = 2

; ================ path file & file search ================ 
path_work = '/work1/LUT/SP/EWwork/'
path_mcd = '/work1/LUT/SP/mcd12/'
path_sav = '/data2/omega/sav/'
path_index = '/work1/LUT/SP/EW_INDEX/'
path_common = '/work1/LUT/Common/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

; MY 27 Ls plot
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

; =============== index 作成 =================
ls_ind = dblarr(loop_a - loop_b +1)
ps_ind = dblarr(loop_a - loop_b +1)
lt_ind = dblarr(loop_a - loop_b +1)
ps_std = dblarr(loop_a - loop_b +1)
dust_ind = dblarr(loop_a - loop_b +1)
albedo_ind = dblarr(loop_a - loop_b +1)

; =============== loop start =================
for loop = loop_b, loop_b do begin

    ; +++++ restore file +++++
    file = files(loop)
    restore, file

    ; n_indを格納
    ip = n_elements(LATI(*,0))
    io = n_elements(LATI(0,*))

    ; 軌道番号を文字として保存
    fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))
    input_file = strupcase(strmid(file,9,4,/REVERSE_OFFSET))
    qub_file = strupcase(strmid(file,9,2,/REVERSE_OFFSET))

    ; QLが作成されなかったデータは弾く
    skip_ind = where(pressure eq 0)
    n_skip = n_elements(skip_ind)
    skip_ratio = n_skip *100d / (ip * io)
    if skip_ratio gt 97d then print, 'warning: skip_ratio'
    if skip_ratio gt 97d then goto, skip1
    
    ; ============== data quality 判定 ===================
    ; fileの開く場所の絶対pathを指定する
    if qub_file ge 00 and qub_file le 25 then path_qub = '/data2/omega/PSA/MEX-M-OMEGA-2-EDR-FLIGHT-V1.0/DATA/ORB' + qub_file + '/' + fileorbit + '.QUB'
    if qub_file ge 26 and qub_file le 48 then path_qub = '/data2/omega/PSA/MEX-M-OMEGA-2-EDR-FLIGHT-EXT1-V1.0/DATA/ORB' + qub_file + '/' + fileorbit + '.QUB'
    if qub_file ge 49 and qub_file le 76 then path_qub = '/data2/omega/PSA/MEX-M-OMEGA-2-EDR-FLIGHT-EXT2-V1.0/DATA/ORB' + qub_file + '/' + fileorbit + '.QUB'
    
    openr,2,path_qub,ERROR=errflag
    close,2
    readcube, path_qub, idat, sdat0, sdat1,info
    data_quality=fix(info(5)) 

    ; data qualityを判定して、必要ないものはskipする
    if input_file le 1900 then begin
        if data_quality lt 2 then print, 'data_quality is Bad'
        if data_quality lt 2 then goto, skip1
    endif

    if input_file ge 1900 then begin
        if data_quality le 3 then print, 'data_quality is Bad'
        if data_quality le 3 then goto, skip1
    endif

    ; ========== limb data の判定をする==============
    restore, path_sav + fileorbit + '.sav'
    tgh = reform(geocube(*,12,*))/1000. ;km
    ind_limb = where_xyz(tgh le 65.535,xind=xind, yind=yind)
    if ind_limb(0) eq -1 then print, 'This is all limb data'
    if ind_limb(0) eq -1 then goto, skip1

    ; =========== 任意の緯度選択 =============
    minlon = min(longi)
    maxlon = max(longi)
    minlat = 30
    maxlat = 35
    str_min = strupcase(strmid(minlat,1,4,/REVERSE_OFFSET))
    str_max = strupcase(strmid(maxlat,1,4,/REVERSE_OFFSET))
    
    ; 上記の任意の緯度経度内の圧力を持ってくる
    ind = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat and lati le maxlat, xind=xind, yind=yind)
    if ind(0) eq -1 then print, 'index error: goto skip'
    if ind(0) eq -1 then goto, skip1

    ; loopを回す場所を限定的にする
    ip1 = min(xind)
    ip2 = max(xind)
    io1 = min(yind)
    io2 = max(yind)
    
    ; indexが0だったときに間違ったloopを回さないようにするためのおまじない
    if ip2 eq 0 then ip2 = ip2 + 1
    if io2 eq 0 then io2 = io2 + 1

    ; =========== 必要なデータをここで格納する ===========
    ; retrieval surface pressure [Done altitude correction]
    p1 = p1_slev

    ; 緯度経度がある場所だけを格納させる
    new_p1 = dblarr(ip,io)
    new_p1(ind) = p1(ind)
    for i = 0, ip-1 do for j = 0, io-1 do if new_p1(i,j) eq 0 then new_p1(i,j) = !VALUES.F_NAN

    ; altitude
    alt = altitude
    
    ; 高度が±10 kmの場合はNanが入るように設定する
    for i = ip1, ip2-1 do for j = io1, io2- 1 do if alt(i,j) gt 10000d then new_p1(i,j) = !VALUES.F_NAN
    for i = ip1, ip2-1 do for j = io1, io2- 1 do if alt(i,j) lt -10000d then new_p1(i,j) = !VALUES.F_NAN

    ; restore MCD LT 12:00 reference
    restore, path_mcd + fileorbit + '.sav'

    ; mcd surface pressure [Done altitude correction]
    mcd_alt = mcd_ps_alt

    ; mcd LT12 surface pressure [Done altitude correction]
    mcd12_alt = mcd_ps12_alt

    ; solar longitude
    LS = ls

    ; local time
    localtime = local_time

    ; ========== 1 sigmaから外れるデータを外す ============
    pm = 1 ; 何σに入るかをここで定義する
    std = stddev(new_p1,/Nan)
    medi = median(new_p1)

    min_range = medi - (pm * std)
    max_range = medi + (pm * std)

    new_p2 = dblarr(ip,io)
    ind2 = where_xyz(new_p1 ge min_range and new_p1 le max_range, xind=xind, yind=yind)
    new_p2(ind2) = new_p1(ind2)

    for i = 0, ip-1 do for j = 0, io-1 do if new_p2(i,j) eq 0 then mcd_alt(i,j) = !VALUES.F_NAN
    for i = 0, ip-1 do for j = 0, io-1 do if new_p2(i,j) eq 0 then new_p2(i,j) = !VALUES.F_NAN

    ret_div = ((new_p2 - mcd12_alt) /mcd12_alt)*100d
    mcd_div = ((mcd_alt - mcd12_alt) /mcd12_alt)*100d


    stop

    for i = ip1, ip2-1 do begin
        for j = io1, io2-1 do begin


        endfor
    endfor

    stop

    skip1:




endfor
stop


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
plot, findgen(100),xs=1, ys=1, yr=[min(lati), max(lati)], xr=[min(longi), max(longi)], back=255, color=0, /nodata, charsize=3
for i = 0, ip-1 do for j = 0, io-1 do if new_p2(i,j) gt 0 then plots, longi(i,j), lati(i,j), color=(new_p2(i,j)-min(new_p2))/(max(new_p2)-min(new_p2))*254., psym=8, symsize=1

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