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
; 
;
;------------------------------------------------------------------------------------------------------------------------

function data_selection, index, min_latitude, max_latitude

loop_index = index
n_loop_index = n_elements(loop_index)

; ================ path file & file search ================ 
path_work = '/work1/LUT/SP/EWwork/'
path_mcd = '/work1/LUT/SP/mcd12/'
path_sav = '/data2/omega/sav/'
path_index = '/work1/LUT/SP/EW_INDEX/'
path_common = '/work1/LUT/Common/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

; 値の初期化
sum_good = 0
aaa = dblarr(10,10)

for loop = 0, n_loop_index -1 do begin

    ; +++++ restore file and restore file information +++++++
    file = files(loop_index(loop))
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
    ;if skip_ratio gt 97d then print, 'warning: skip_ratio'
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
        ;if data_quality lt 2 then print, 'data_quality is Bad'
        if data_quality lt 2 then goto, skip1
    endif

    if input_file ge 1900 then begin
        ;if data_quality le 3 then print, 'data_quality is Bad'
        if data_quality le 3 then goto, skip1
    endif

    ; ========== limb data の判定をする==============
    restore, path_sav + fileorbit + '.sav'
    tgh = reform(geocube(*,12,*))/1000. ;km
    ind_limb = where_xyz(tgh le 65.535,xind=xind, yind=yind)
    if ind_limb(0) eq -1 then print, 'This is all limb data'
    if ind_limb(0) eq -1 then goto, skip1

    ; =========== 任意の緯度選択 =============
    ;minlon = min(longi)
    ;maxlon = max(longi)
    minlon = 100
    maxlon = 180
    minlat = min_latitude
    maxlat = max_latitude
    str_min = strupcase(strmid(minlat,1,4,/REVERSE_OFFSET))
    str_max = strupcase(strmid(maxlat,1,4,/REVERSE_OFFSET))
    
    ; 上記の任意の緯度経度内の圧力を持ってくる
    ind = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat and lati le maxlat, xind=xind, yind=yind)
    ;if ind(0) eq -1 then print, 'index error: goto skip'
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

    ; limb データは弾くようにする
    for i = ip1, ip2-1 do for j = io1, io2- 1 do if tgh(i,j) ge 65.535 then new_p1(i,j) = !VALUES.F_NAN

    ; +++++ dust / water ice / ice cloud / Albedo index +++++
    ; 格納されているINDEX: rad_2777 (dust) , R_surface_CO2ice (surface ice) , ice_index (ice cloud) , albedo (albedo) , O2, CO, CO2
    restore, path_index + 'INDEX_' + fileorbit + '.sav'

    ; dust index
    ; 0.4以上はダストがあるとみなして、使用しない
    dust_index = rad_2777
    for i = ip1, ip2-1 do for j = io1, io2- 1 do if dust_index(i,j) ge 0.4d then new_p1(i,j) = !VALUES.F_NAN

    ; ice cloud index
    ; 0.7以下は氷があるとみなして、pixelを除外する
    ice_index = ice_index
    for i = ip1, ip2-1 do for j = io1, io2- 1 do if ice_index(i,j) le 0.7d then new_p1(i,j) = !VALUES.F_NAN

    ; albedo index
    ; albedo 0.1から0.5までを採用
    albedo_index = albedo
    for i = ip1, ip2-1 do for j = io1, io2- 1 do if albedo_index(i,j) ge 0.5d then new_p1(i,j) = !VALUES.F_NAN
    for i = ip1, ip2-1 do for j = io1, io2- 1 do if albedo_index(i,j) le 0.1d then new_p1(i,j) = !VALUES.F_NAN

    ; restore MCD LT 12:00 reference
    restore, path_mcd + fileorbit + '.sav'

    ; mcd surface pressure [Done altitude correction]
    mcd_alt = mcd_ps_alt

    ; mcd LT12 surface pressure [Done altitude correction]
    mcd12_alt = mcd_ps12_alt

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
    
    ; 作成したスペクトルをここに格納する
    good_ind = where(abs(ret_div) gt 0)
    n_good = n_elements(good_ind)

    ; 足し合わせるminimumをここで定義しておく
    min_sum = sum_good
    if min_sum gt 0 then new_ret_div_old = new_ret_div
    if min_sum gt 0 then new_mcd_div_old = new_mcd_div

    ; 総スペクトル数を計算する
    sum_good = sum_good + n_good

    ; retrieval spectrum
    new_ret_div = dblarr(sum_good)
    if min_sum gt 0 then new_ret_div(0:min_sum-1) =  new_ret_div_old
    new_ret_div(min_sum:sum_good-1) = ret_div(good_ind)

    ; mcd data
    new_mcd_div = dblarr(sum_good)
    if min_sum gt 0 then new_mcd_div(0:min_sum-1) =  new_mcd_div_old
    new_mcd_div(min_sum:sum_good-1) = mcd_div(good_ind)

    skip1:

endfor

print, 'sample number: ', sum_good
print, 'file orbit: ', fileorbit
if sum_good eq 0 then return, aaa

matrix = reform([new_ret_div,new_mcd_div],sum_good,2)
return, matrix

end



pro tidalwave_analysis_v5_2

device, decomposed = 0, retain = 2

; ================ path + Ls/ localtime index serach ================ 
path_common = '/work1/LUT/Common/'
path_data = '/work1/LUT/SP/data_analysis_index/'

file = path_common + 'Ls_Localtime_index.sav'
restore, file

; ====== plotの枠組み ========
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

;; plotの枠組みを作成  
;window, 0, xs=1500, ys=1000
;loadct,39
;plot, findgen(100),xs=1, ys=1, yr=[-30, 30], xr=[0, 24], back=255, color=0, /nodata, charsize=3, title='tidel wave detection', xtitle='Local time', ytitle='mean pressure [%]'

; Ls/ local time indexからどの軌道番号を使用するかを選択する
Ls_ind = ls_ind
lt_ind = lt_ind

; +++++++++ Ls / Local time / latitude 指定 ++++++++++++++
; Lsの定義
min_Ls = 270
max_Ls = 360
str_min_ls = strupcase(strmid(min_ls,1,4,/REVERSE_OFFSET))
str_max_ls = strupcase(strmid(max_ls,1,4,/REVERSE_OFFSET))

; local timeの定義
;start_local = 0d
dev_local = 1d  ; Local time 何分刻みにしたいかをここで指定

for i = 0, 24d, dev_local do begin
    str_i = strupcase(strmid(i,1,4,/REVERSE_OFFSET))
    index = where(ls_ind ge min_Ls and ls_ind le max_Ls and lt_ind gt i and lt_ind lt i + dev_local)
    if index(0) eq -1 then print, 'skip', i
    if index(0) eq -1 then goto, skip2

    ; 緯度を定義
    min_latitude = 30
    max_latitude = 70
    str_min = strupcase(strmid(min_latitude,1,4,/REVERSE_OFFSET))
    str_max = strupcase(strmid(max_latitude,1,4,/REVERSE_OFFSET))

    ; 熱潮汐のデータ解析に必要なデータの選択を行う
    spectrum = data_selection(index, min_latitude, max_latitude)
    if spectrum(0,0) eq 0 then print, 'skip', i
    if spectrum(0,0) eq 0 then goto, skip2

    save, spectrum, filename = path_data + 'local_' + str_i + '_lat' + str_min + '_' + str_max + '_' + 'ls_' + str_min_ls + '_' + str_max_ls + '_lon.sav'

    ; 熱潮汐の検出のための解析をここで行う
    ; OMEGAのリトリーバルデータ
    ret_data = spectrum(*,0)
    medi_ret = median(ret_data)
    std_ret = stddev(ret_data,/nan)

    ret_ran = [medi_ret - std_ret, medi_ret + std_ret]

    ; mcd データ
    mcd_data = spectrum(*,1)
    medi_mcd = median(mcd_data)
    std_mcd = stddev(mcd_data,/nan)

    mcd_ran = [medi_mcd - std_mcd, medi_mcd + std_mcd]

    ;plots, i, medi_ret, color=0, psym=8, symsize=2
    ;plots, i, ret_ran, color=0, THICK=3
    
    ;plots, i, medi_mcd, color=254, psym=8, symsize=2
    ;plots, i, mcd_ran, color=254, THICK=3

    skip2:

    
endfor


stop

end