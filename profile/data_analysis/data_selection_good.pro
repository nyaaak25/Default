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
; albedo, water ice, dust index use
;
; tidalwave_analysis_v5_3　 ::2023.7.6 Thu 12:37:00
; 熱潮汐波の検出をする
; function化 + Ls/ local time indexを付加
; 平均方法：各スペクトルで中央値を探す
; データクオリティ2/3以上のものを使用する, 1σに入らないデータは除く
; データセレクション：高度±10 km, リムデータを除く
; albedo, water ice, dust index use
; 
; +++ Ls/LTをきちんと判断する
;
;------------------------------------------------------------------------------------------------------------------------

pro data_selection_good

device, decomposed = 0, retain = 2

; ================ path + Ls/ localtime index serach ================ 
path_work = '/work1/LUT/SP/EWwork/'
pathwork = '/work1/LUT/SP/EWwork_good/'
path_sav = '/data2/omega/sav/'
path_index = '/work1/LUT/SP/EW_INDEX/'

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

for loop = loop_b, loop_a -1 do begin

    ; +++++ restore file and restore file information +++++++
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
    if skip_ratio gt 97d then print, 'warning: skip_ratio', fileorbit
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
        if data_quality lt 2 then print, 'data_quality is Bad', fileorbit
        if data_quality lt 2 then goto, skip1
    endif

    if input_file ge 1900 then begin
        if data_quality le 3 then print, 'data_quality is Bad', fileorbit
        if data_quality le 3 then goto, skip1
    endif

    ; ========== limb data の判定をする==============
    restore, path_sav + fileorbit + '.sav'
    tgh = reform(geocube(*,12,*))/1000. ;km
    ind_limb = where_xyz(tgh le 65.535,xind=xind, yind=yind)
    if ind_limb(0) eq -1 then print, 'This is all limb data', fileorbit
    if ind_limb(0) eq -1 then goto, skip1

    save, lati, longi, altitude, trans, pressure, InputAlbedo, dustmap, TAmap, MCDpressure, SZA_all, EA_all, PA_all, local_time, SOLAR_LONGITUDE, p1_slev, bad_index, bad_ratio, SZA_ratio, EA_ratio, dust_ratio, filename = pathwork + 'GOOD_' + fileorbit + '.sav'


    skip1:

endfor

stop

end