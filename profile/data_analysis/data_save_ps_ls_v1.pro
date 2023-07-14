;------------------------------------------------------------------------------------------------------------------------
;
; 熱潮汐波の解析をするためのプログラム
;
; create by Akira Kazama
; 
; data_save_ps_ls　 ::2023.6.13 Tue 12:00:00
; MCD pressureをsaveするプログラム
;
; data_save_ps_ls_v1　 ::2023.6.20 Tue 10:01:00
; 高度補正を行った平均化圧力をsave fileに格納するプログラム
;
;------------------------------------------------------------------------------------------------------------------------


pro data_save_ps_ls_v1

device, decomposed = 0, retain = 2

; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'

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

; ==== index 作成 =====
ls_ind = dblarr(loop_a - loop_b +1)
ps_ind = dblarr(loop_a - loop_b +1)
lt_ind = dblarr(loop_a - loop_b +1)
ps_std = dblarr(loop_a - loop_b +1)
dust_ind = dblarr(loop_a - loop_b +1)
albedo_ind = dblarr(loop_a - loop_b +1)


; ===== fileのrestore loopを始める =======
for loop = loop_b, loop_a do begin
    
    ; initial data
    std = -0d/0d
    ls = -0d/0d
    various_p1 = -0d/0d
    localtime = -0d/0d
    medi_dust = -0d/0d
    medi_albedo = -0d/0d
    
    ; fileのrestore
    file = files(loop)
    restore, file

    ; n_indを格納
    ip = n_elements(LATI(*,0))
    io = n_elements(LATI(0,*))


    ; 軌道番号を文字として保存
    fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))
    input_file = strupcase(strmid(file,9,4,/REVERSE_OFFSET))
    qub_file = strupcase(strmid(file,9,2,/REVERSE_OFFSET))

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

    ; ===================================================

    ; longitude, latitude
    lat1 = lati
    lon1 = longi

    ; ls
    ls = SOLAR_LONGITUDE

    ; ====== 任意の緯度選択 ======
    minlon = min(lon1)
    maxlon = max(lon1)
    minlat = 0
    maxlat = 30
    str_min = strupcase(strmid(minlat,1,4,/REVERSE_OFFSET))
    str_max = strupcase(strmid(maxlat,1,4,/REVERSE_OFFSET))
    
    ; 上記の任意の緯度経度内の圧力を持ってくる
    ind = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat and lati le maxlat, xind=xind, yind=yind)
    if ind(0) eq -1 then print, 'index error: goto skip'
    if ind(0) eq -1 then goto, skip1

    ; create index
    mcd_ps_12 = dblarr(ip,io)
    mcd_ta_12 = dblarr(ip,io)
    
    ; loopを回す場所を限定的にする
    ip1 = min(xind)
    ip2 = max(xind)
    io1 = min(yind)
    io2 = max(yind)
    
    ; indexが0だったときに間違ったloopを回さないようにするためのおまじない
    if ip2 eq 0 then ip2 = ip2 + 1
    if io2 eq 0 then io2 = io2 + 1

    ; albedo map index
    albedo = inputalbedo

    ; dust map index
    dust = dustmap

    ; retrieval surface pressure [Done altitude correction]
    p1 = p1_slev

    ; bad pixelにNanを格納
    for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then dust(i,j) = !VALUES.F_NAN
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

    ; その中で平均を取る
    various_p1 = mean(new_p1(ind), /Nan)
    std = stddev(new_p1(ind), /Nan)
    localtime = median(local_time)
    medi_dust = median(dust)
    medi_albedo = median(albedo)

    skip1:

    ;if fileorbit ge 0181 and fileorbit le 2595 then ls_ind(loop) = ls
    ;if fileorbit ge 2607 and fileorbit le 5055 then ls_ind(loop) = ls + 360d
    ;if fileorbit ge 5062 and fileorbit le 7454 then ls_ind(loop) = ls + 720d

    ls_ind(loop) = ls
    ps_ind(loop) = various_p1
    ps_std(loop) = std
    lt_ind(loop) = localtime
    dust_ind(loop) = medi_dust
    albedo_ind(loop) = medi_albedo

    print, "loop:", (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1


endfor
save, ls_ind, ps_ind, lt_ind, ps_std, dust_ind, albedo_ind, filename = path_work + 'ret_latitude_' + str_min + '_' + str_max + '_1new_std_DQ.sav'
stop


end