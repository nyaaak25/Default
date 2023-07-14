;------------------------------------------------------------------------------------------------------------------------
;
; 熱潮汐波の解析をするためのプログラム
;
; create by Akira Kazama
; 
; data_save_ps_ls　 ::2023.6.13 Tue 12:00:00
; MCD pressureをsaveするプログラム

; data_save_ps_ls_1　 ::2023.7.4 Tue 13:35:00
; MCD pressureをsaveするプログラム
; MCD LT12:00の圧力をORBごとに作成するプログラム
;
; data_save_ps_ls_v1　 ::2023.6.20 Tue 10:01:00
; 平均化した圧力をsave fileに格納するプログラム
;
;------------------------------------------------------------------------------------------------------------------------


pro data_save_ps_ls_1

device, decomposed = 0, retain = 2

; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
pathwork = '/work1/LUT/SP/mcd12/'
path_ql = '/work1/LUT/SP/QL_datacover/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

; MY 27 Ls plot
;loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB0181_0.sav') ; MY27 Ls 0
;loop_b = fix(loop_b(0))
;loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB1474_5.sav') ; MY27 Ls 360
;loop_a = fix(loop_a(0))

;loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB1478_1.sav') ; MY27 Ls 0
;loop_b = fix(loop_b(0))
;loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB2595_0.sav') ; MY27 Ls 360
;loop_a = fix(loop_a(0))

; MY 28 Ls plot
;loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB2607_0.sav') ; MY28 Ls 0
;loop_b = fix(loop_b(0))
;loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB4551_4.sav') ; MY28 Ls 360
;loop_a = fix(loop_a(0))

;loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB4553_1.sav') ; MY28 Ls 0
;loop_b = fix(loop_b(0))
;loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB5055_4.sav') ; MY28 Ls 360
;loop_a = fix(loop_a(0))

; MY 29 Ls plot
;loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB5062_0.sav') ; MY29 Ls 0
;loop_b = fix(loop_b(0))
;loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB5814_3.sav') ; MY29 Ls 360
;loop_a = fix(loop_a(0))

loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB5826_1.sav') ; MY29 Ls 0
loop_b = fix(loop_b(0))
loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB7454_4.sav') ; MY29 Ls 360
loop_a = fix(loop_a(0))


; ===== altitude correction ======
g = 3.72 ;m/s2
Gconst = 6.67430d-11
MMars = 6.42d23
RMars = 3.4d6
g = -Gconst*MMars/(-RMars*RMars)
R = 192d ;J/K/kg

; ===== fileのrestore loopを始める =======
for loop = loop_b, loop_a do begin
    
    ; initial data
    ls = -0d/0d
    
    ; fileのrestore
    file = files(loop)
    restore, file
    fileorbit = strupcase(strmid(file,9,4,/REVERSE_OFFSET))
    file_orbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))

    ; n_indを格納
    ip = n_elements(LATI(*,0))
    io = n_elements(LATI(0,*))

    skip_ind = where(pressure eq 0)
    n_skip = n_elements(skip_ind)
    skip_ratio = n_skip *100d / (ip * io)
    if skip_ratio gt 97d then print, 'warning: skip_ratio'
    if skip_ratio gt 97d then goto, skip1

    ; longitude, latitude
    lat1 = lati
    lon1 = longi

    ; ls
    ls = SOLAR_LONGITUDE

    ; ====== 任意の緯度選択 ======
    minlon = min(lon1)
    maxlon = max(lon1)
    minlat = min(lat1)
    maxlat = max(lat1)
    str_min = strupcase(strmid(minlat,1,4,/REVERSE_OFFSET))
    str_max = strupcase(strmid(maxlat,1,4,/REVERSE_OFFSET))
    
    ; 上記の任意の緯度経度内の圧力を持ってくる
    ind = where_xyz(longi ge minlon and longi le maxlon and lati ge minlat and lati le maxlat, xind=xind, yind=yind)
    if ind(0) eq -1 then print, 'index error: goto skip'
    if ind(0) eq -1 then goto, skip1

    ; create index
    mcd_ps_12 = dblarr(ip,io)
    mcd_ta_12 = dblarr(ip,io)

    ; retreival surface pressure altitude correct
    p10 = p1_slev
    ret_ps_alt = p1_slev
    
    ; loopを回す場所を限定的にする
    ip1 = min(xind)
    ip2 = max(xind)
    io1 = min(yind)
    io2 = max(yind)

    ; dust indexを作成する
    if fileorbit ge 0181 and fileorbit le 2595 then d_index = 27
    if fileorbit ge 2607 and fileorbit le 5055 then d_index = 28
    if fileorbit ge 5062 and fileorbit le 7454 then d_index = 29
    
    ; indexが0だったときに間違ったloopを回さないようにするためのおまじない
    if ip2 eq 0 then ip2 = ip2 + 1
    if io2 eq 0 then io2 = io2 + 1

    ; MCD loct.12 sea leveled pressure
    for l = ip1, ip2 -1 do begin ;loop for slit scan
        for k = io1, io2 -1 do begin ;test getting surface feature
            if p10(l,k) eq 0 then goto, skip10

        ; =========== MCD version 6.1 ===========
            lat = lati(l,k)
            lon = longi(l,k)
            Ls = SOLAR_LONGITUDE
            Loct = 12d
            hrkey = 1 ; set high resolution mode on (hrkey=0 to set high resolution off)
            zkey = 3    ; specify that xz is the altitude above surface (m)
            xz = 0. ; (array of altitude: step=2000 m start: 5m)
            datekey = 1       ; <integer> type of input date (1=Mars date)
            xdate = Ls        ; <double precision> date (IF datekey = 1 : Value of Ls [deg.])
            dset ='/work1/LUT/MCDv6.1/MCD_6.1_queleclim/data/'
            ;dset ='/work1/LUT/MCD/MCD5.3/data/' ;‘MCD_DATA/’  ; <character*50> data set
            scena = d_index         ; <integer> scenario (1 = Climatology ave solar)
            perturkey = 1     ; <integer>  perturbation type (1= none)
            seedin   = 7.0    ; <real>
            gwlength = 0.0    ; <real>  for small scale (ie: gravity wave) perturbations;
            extvarkeys = LONARR(101) ; add 1 element because indexes in IDL start at 0
            for i0 = 0, 100 do extvarkeys[i0] = 1 ; <integer> array output type (extvar(i) = 0 : don’t compute)
            meanvar = FLTARR(6) ; add 1 element because indexes in IDL start at 0
            for i0 = 0, 5 do meanvar[i0] = 0.0 ; <real> mean unperturbed values (array of 5)
            extvar = FLTARR(101) ; add 1 element because indexes in IDL start at 0
            for i0 = 0, 100 do extvar[i0] = 0.0 ; <real>  extra variables (array of 100)
            pres = 0.0        ; <real> atmospheric pressure (Pa)
            dens = 0.0        ; <real> atmospheric density (kg/m^3)
            temp = 0.0        ; <real> atmospheric temperature (K)
            zonwind = 0.0     ; <real> zonal wind component (East-West)
            merwind = 0.0     ; <real> meridional wind component (North-South)
            seedout = 0.0     ; <real> current value of the seed of the random number generator
            ierr = 0          ; <integer> error flag (0 = OK)
            ;    cd, ‘/Users/Shohei/Documents/Programs/MCD_ver4.3/mcd/idl/’
            ;   mcd_idl, Ls, Loct, lat, lon, dust, zkey, hrkey, xz, meanvarz, extvarz
            cd, '/work1/LUT/MCDv6.1/MCD_6.1_queleclim/mcd/interfaces/idl/'
            ;cd, '/work1/LUT/MCD/MCD5.3/mcd/idl/'
            a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
            dset,scena,perturkey,seedin,gwlength,extvarkeys, $
            pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)

            scaleH = extvar(63)
            MOLA_altitude = extvar(4)
            Z1 = scaleH * 0.1d
            Z2 = scaleH * 4d
            dust_opacity = extvar(40) * 1.2090217E+00/4.6232791E+00
            ice_opacity = 0.d;*4.2410289E-01/***  ;!TBD!
            SP = extvar(21)
            mcd_ps_12(l,k) = pres

            xz=[Z1] ; (array of altitude with step=2000 m)
            a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
            dset,scena,perturkey,seedin,gwlength,extvarkeys, $
            pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
            TA = temp
            mcd_ta_12(l,k) = TA

            skip10:

        endfor
    endfor

    ; 選択された緯度・経度内のリトリーバル結果が全て0だった場合は、高度補正をしないでskipする
    zero_index = where(mcd_ps_12 eq 0)
    zero = n_elements(zero_index)
    zero_ratio = zero *100d / (ip *io)
    if zero_ratio eq 100d then print, 'warnig: zero skip'
    if zero_ratio eq 100d then goto, skip1
    
    
    ; retrieval surface pressure [Done altitude correction]
    ;p1 = p1_slev

    ; bad pixelにNanを格納
    ;for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN

    ; MCD surface pressure
    p1_mcd = mcdpressure 

    ; bad pixelにNanを格納
    for i = 0, ip-1 do for j = 0, io-1 do if p1_mcd(i,j) eq 0 then p1_mcd(i,j) = !VALUES.F_NAN

    ; mcd pressure altutude corecction
    mcd_ps_alt = dblarr(ip,io)
    mcd_ps12_alt = dblarr(ip,io)
    p5 = altitude
    TAmap = TAMAP

    for i = 0, ip-1 do for j = 0, io-1 do if p1_mcd(i,j) gt 0 then mcd_ps_alt(i,j) = p1_mcd(i,j) * exp( (p5(i,j)) / (R*TAmap(i,j)/(-Gconst*MMars/(-1d*(RMars+p5(i,J))*(RMars+p5(i,J)) ))))
    for i = 0, ip-1 do for j = 0, io-1 do if mcd_ps_alt(i,j) eq 0 then mcd_ps_alt(i,j) = !VALUES.F_NAN

    ; sea level mcd 12. pressure
    for i = 0, ip-1 do for j = 0, io-1 do if mcd_ps_12(i,j) gt 0 then mcd_ps12_alt(i,j) = mcd_ps_12(i,j) * exp( (p5(i,j)) / (R*mcd_ta_12(i,j)/(-Gconst*MMars/(-1d*(RMars+p5(i,J))*(RMars+p5(i,J)) ))))
    for i = 0, ip-1 do for j = 0, io-1 do if mcd_ps12_alt(i,j) eq 0 then mcd_ps12_alt(i,j) = !VALUES.F_NAN


    save, mcd_ps_12, mcd_ta_12, mcd_ps12_alt, mcd_ps_alt, ret_ps_alt, ls, longi, lati, local_time, filename = pathwork + file_orbit + '.sav'

    skip1:
    print, "loop:", (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1

endfor
stop


end