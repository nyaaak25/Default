;------------------------------------------------------------------------------------------------------------------------
;
; 熱潮汐波の解析をするためのプログラム
;
; create by Akira Kazama
; 
; tidalwave_analysis　 ::2023.5.13 Sat 07:42:00
;
;------------------------------------------------------------------------------------------------------------------------


pro tidalwave_analysis

device, decomposed = 0, retain = 2

; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'
path_sav = '/data2/omega/sav/'
path_index = '/work1/LUT/SP/EW_INDEX/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB0920_0.sav') ; MY27 Ls 0
loop_b = fix(loop_b(0))
loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB0938_6.sav') ; MY27 Ls 360
loop_a = fix(loop_a(0))

; MY 27 Ls plot
;loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB0181_0.sav') ; MY27 Ls 0
;loop_b = fix(loop_b(0))
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
;loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB7454_4.sav') ; MY29 Ls 360
;loop_a = fix(loop_a(0))

; create index
Ls_ind = dblarr(loop_a - loop_b)
mean_p1 = dblarr(loop_a - loop_b)

; ====== 任意の季節を選択 ======
; Ls
minls = 0
maxls = 360

; ===== curiosity data ======
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

;file_curiosity = '/work1/LUT/Soft/curiosity.dat'
;openr,lun,file_curiosity,/get_lun
;for i = 0, 25-1 do begin
;    readf,lun,a,b
;    curiosity_ls(i)=a
;    curiosity_ps(i)=b
;endfor

;file_curiosity_ta = '/work1/LUT/Soft/curiosity_ta.dat'
;openr,lun,file_curiosity_ta,/get_lun
;for i = 0, 25-1 do begin
;    readf,lun,a
;    curiosity_ta(i)=a
;endfor


;cs_altitude = -4488.0483d

;for i = 0, 25-1 do cs_slev(i) = curiosity_ps(i) * exp( (cs_altitude) / (R*curiosity_ta(i)/(-Gconst*MMars/(-1d*(RMars+cs_altitude)*(RMars+cs_altitude) ))))


; プロットの枠組み
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

;; plotの枠組みを作成
window, 0, xs=1500, ys=1000
loadct,39
;plot, findgen(100),xs=1, ys=1, yr=[290, 1000], xr=[minls, maxls], back=255, color=0, /nodata, charsize=3, title='MY29 Retrieved pressure Ls ' +strmid(minls,5,6)+ ' to ' +strmid(maxls,5,6), xtitle='Ls [Deg]', ytitle='mean pressure [Pa]'
plot, findgen(100),xs=1, ys=1, yr=[400, 700], xr=[96, 102], back=255, color=0, /nodata, charsize=3, xtitle='Ls [Deg]', ytitle='mean pressure [Pa]'

; fileのrestore loopを始める
for loop = loop_b, loop_a do begin
    
    ; fileのrestore
    file = files(loop)
    restore, file

    ; 軌道番号を文字として保存
    fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))
    input_file = strupcase(strmid(file,9,4,/REVERSE_OFFSET))
    qub_file = strupcase(strmid(file,9,2,/REVERSE_OFFSET))

    ; n_indを格納
    ip = n_elements(LATI(*,0))
    io = n_elements(LATI(0,*))

    skip_ind = where(pressure eq 0)
    n_skip = n_elements(skip_ind)
    skip_ratio = n_skip *100d / (ip * io)
    if skip_ratio gt 97d then goto, skip12


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
        if data_quality lt 2 then goto, skip12
    endif

    if input_file ge 1900 then begin
        ;if data_quality le 3 then print, 'data_quality is Bad'
        if data_quality le 3 then goto, skip12
    endif

    ; ========== limb data の判定をする==============
    restore, path_sav + fileorbit + '.sav'
    tgh = reform(geocube(*,12,*))/1000. ;km
    ind_limb = where_xyz(tgh le 65.535,xind=xind, yind=yind)
    if ind_limb(0) eq -1 then print, 'This is all limb data'
    if ind_limb(0) eq -1 then goto, skip12

    ; longitude, latitude
    lat1 = lati
    lon1 = longi

    ; altitude +- 10 km
    alt = altitude

    ; ls
    ls = SOLAR_LONGITUDE

    ; retrieval surface pressure [Done altitude correction]
    p1 = p1_slev
    
    ; 高度が±10 kmの場合はNanが入るように設定する
    for i = 0, ip-1 do begin
        for j = 0, io- 1 do begin

            ; bad pixelにnanを格納
            if p1(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN

            ; 高度が±10 kmの場合はNanが入るように設定する
            if alt(i,j) gt 10000d then p1(i,j) = !VALUES.F_NAN
            if alt(i,j) lt -10000d then p1(i,j) = !VALUES.F_NAN

            ; limb dataを除外する
            if tgh(i,j) ge 65.535 then p1(i,j) = !VALUES.F_NAN

        endfor
    endfor

    ; +++++ INDEXを使用してデータを選定する ++++++++++
    restore, path_index + 'INDEX_' + fileorbit + '.sav'

    ; dust index
    ; 0.4以上はダストがあるとみなして、使用しない
    dust_index = rad_2777

    ; ice cloud index
    ; 0.7以下は氷があるとみなして、pixelを除外する
    ice_index = ice_index

    ; albedo index
    ; albedo 0.1から0.5までを採用
    albedo_index = albedo

    ; surface ice index
    R_surface_CO2ice = R_surface_CO2ice

    for i = 0, ip-1 do begin
        for j = 0, io- 1 do begin

            ; dust index
            if dust_index(i,j) ge 0.4d then p1(i,j) = !VALUES.F_NAN

            ; ice cloud
            if ice_index(i,j) le 0.7d then p1(i,j) = !VALUES.F_NAN

            if R_surface_CO2ice(i,j) ge 1.0d then p1(i,j) = !VALUES.F_NAN

            ; albedo index
            if albedo_index(i,j) ge 0.4d then p1(i,j) = !VALUES.F_NAN
            if albedo_index(i,j) le 0.2d then p1(i,j) = !VALUES.F_NAN


        endfor
    endfor

    ; ここで1σに入らないデータは弾くことをする
    pm = 1 ; 何σに入るかをここで定義する
    std = stddev(p1,/Nan)
    medi = median(p1)

    min_range = medi - (pm * std)
    max_range = medi + (pm * std)

    new_p1 = dblarr(ip,io)
    ind2 = where_xyz(p1 ge min_range and p1 le max_range, xind=xind, yind=yind)
    new_p1(ind2) = p1(ind2)

    for i = 0, ip-1 do for j = 0, io-1 do if new_p1(i,j) eq 0 then new_p1(i,j) = !VALUES.F_NAN

    ; MCD surface pressure
    ;p1_mcd = mcdpressure 

    ; bad pixelにNanを格納
    ;for i = 0, ip-1 do for j = 0, io-1 do if p1_mcd(i,j) eq 0 then p1_mcd(i,j) = !VALUES.F_NAN

    ; mcd pressure altutude corecction
    ;ip = n_elements(lat1(*,0))
    ;io = n_elements(lat1(0,*))
    ;p1 = dblarr(ip,io)
    ;p5 = altitude
    ;TAmap = TAMAP

    ;for i = 0, ip-1 do for j = 0, io-1 do if p1_mcd(i,j) gt 0 then p1(i,j) = p1_mcd(i,j) * exp( (p5(i,j)) / (R*TAmap(i,j)/(-Gconst*MMars/(-1d*(RMars+p5(i,J))*(RMars+p5(i,J)) ))))
    ;for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN

    ; ====== 任意の緯度選択 ======
    ; lat: -70 to -30、-30 to 0、0 to 30、30 to 70
    minlon = min(lon1)
    maxlon = max(lon1)
    minlat = -30
    maxlat = 0

    ; 上記の任意の緯度経度内の圧力を持ってくる
    ind = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat and lati le maxlat, xind=xind, yind=yind)
    if ind(0) eq -1 then goto, skip1

    ; その中で平均を取る
    ;various_p1 = mean(new_p1(ind), /Nan)
    various_p1 = median(new_p1(ind))
 
    ; plot
    if various_p1 gt 0 and various_p1 lt 1000 and ls gt minls and ls lt maxls then begin
        loadct, 39
        plots, ls, various_p1, color=0, psym=8, symsize=2
        if input_file gt 2000 then plots, ls, various_p1, color=254, psym=8, symsize=2
        print, 'pressure: ', various_p1
        print, 'ls: ', ls
        print, 'local time: ', mean(local_time)
        print, 'ORB number: ', fileorbit
        print, 'orbloop: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1

        ;Ls_ind(loop) = ls
        ;mean_p1(loop) = various_p1

    endif

    skip1:

    ; ======= 領域ごとに色を変化させてみてみる =======
    ; ② lat: 0 to 30
    minlat2 = 0
    maxlat2 = 30

    ind2 = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat2 and lati le maxlat2, xind=xind, yind=yind)
    if ind2(0) eq -1 then goto, skip10
    
    sample2 = 

    ;various_p2 = mean(new_p1(ind2), /Nan)
    various_p2 = median(new_p1(ind2))
 
    if various_p2 gt 0 and various_p2 lt 1000 and ls gt minls and ls lt maxls then begin
        loadct, 39
        plots, ls, various_p2, color=254, psym=8, symsize=2
        ;plots, ls, various_p2, color=0, psym=8, symsize=2
        if input_file gt 2000 then plots, ls, various_p2, color=254, psym=8, symsize=2
        print, 'pressure: ', various_p2
        print, 'ls: ', ls
        print, 'local time: ', mean(local_time)
        print, 'ORB number: ', fileorbit
        print, 'orbloop: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1


    endif

    skip10:

    ; ======= 領域ごとに色を変化させてみてみる =======
    ; ③ lat: 30 to 70

    minlat3 = 30
    maxlat3 = 70

    ind3 = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat3 and lati le maxlat3, xind=xind, yind=yind)
    if ind3(0) eq -1 then goto, skip11

    ;various_p3 = mean(new_p1(ind3), /Nan)
    various_p3 = median(new_p1(ind3))
 
    if various_p3 gt 0 and various_p3 lt 1000 and ls gt minls and ls lt maxls then begin
        loadct, 39
        plots, ls, various_p3, color=60, psym=8, symsize=2
        ;plots, ls, various_p3, color=0, psym=8, symsize=2
        if input_file gt 2000 then plots, ls, various_p3, color=254, psym=8, symsize=2
        print, 'pressure: ', various_p3
        print, 'ls: ', ls
        print, 'local time: ', mean(local_time)
        print, 'ORB number: ', fileorbit
        print, 'orbloop: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1


    endif

    skip11:

    ; ======= 領域ごとに色を変化させてみてみる =======
    ; ④ lat: -70 to -30

    minlat4 = -70
    maxlat4 = -30

    ind4 = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat4 and lati le maxlat4, xind=xind, yind=yind)
    if ind4(0) eq -1 then goto, skip12

    various_p4 = mean(new_p1(ind4), /Nan)
    various_p4 = median(new_p1(ind4))
 
    if various_p4 gt 0 and various_p4 lt 1000 and ls gt minls and ls lt maxls then begin
        loadct, 39
        plots, ls, various_p4, color=130, psym=8, symsize=2
        ;plots, ls, various_p4, color=0, psym=8, symsize=2
        if input_file gt 2000 then plots, ls, various_p4, color=254, psym=8, symsize=2
        print, 'pressure: ', various_p4
        print, 'ls: ', ls
        print, 'local time: ', mean(local_time)
        print, 'ORB number: ', fileorbit
        print, 'orbloop: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1


    endif

    skip12:


endfor
stop

loadct,39
plots, curiosity_ls, cs_slev, color=205, thick = 6
plots, curiosity_ls, cs_slev, color=205, psym = 8, symsize=2

snapshot = TVRD(True=1)
Write_JPEG, path_ql + 'ret_pressure_seasonal_MY29_cs.jpg', snapshot, True=1, Quality=75


stop

end