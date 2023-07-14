;------------------------------------------------------------------------------------------------------------------------
;
; dust / water ice / surface ice: index
;
; create by Akira Kazama
; 
; parameter_index　 ::2023.6.30 Fri 11:36:00
; dust/ water ice/ surface iceを判別するためのプログラム
;
; parameter_index_v1　 ::2023.6.30 Fri 16:45:00
; dust/ water ice/ surface iceを判別するためのプログラム
; retrievalされた圧力のみのQL/WORK fileを作成する
; 
;
;------------------------------------------------------------------------------------------------------------------------

; 青木さんに聞く
function BT,B,lambda
  C1 = 1.1911D-08
  C2 = 1.4387D+00
  T = 1.d + C1*(lambda^3) / B * 1e3
  T = 1.d / alog(T)
  T = T * C2 * lambda
  return,T
end

pro parameter_index_v1

; =========== path ===========
path_file = '/data2/omega/sav/'
pathfile = '/work1/LUT/SP/EWwork/'
path_work = '/work1/LUT/SP/EW_INDEX/'
path_ql = '/work1/LUT/SP/QL_EW_index/'

files=file_search(pathfile + '*.sav',count=count) ; sav fileすべてを持ってくる

loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB6439_1.sav') ; MY27 Ls 0
loop_b = fix(loop_b(0))

; MY27
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
;loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB5887_0.sav') ; MY29 Ls 0
;loop_b = fix(loop_b(0))
loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB7454_4.sav') ; MY29 Ls 360
loop_a = fix(loop_a(0))

; =========== file loop start =========== 
for loop = loop_b, loop_a do begin

    ; ++++++ restore file ++++++
    ; retrieval するfileをrestore
    file = files(loop)
    restore, file

    ; 全部indexを持ってくる
    ip = n_elements(LATI(*,0))
    io = n_elements(LATI(0,*))

    ; getting orb file name
    fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))
    skip_ind = where(pressure eq 0)
    n_skip = n_elements(skip_ind)
    skip_ratio = n_skip *100d / (ip * io)
    if skip_ratio gt 97d then print, 'warning: skip_ratio'
    if skip_ratio gt 97d then goto, skip1

    ; getting sav file
    OMEGA_file = file_search(path_file + fileorbit + '.sav',count=count) 
    restore, OMEGA_file

    ; orb information
    timei=reform(geocube(0:6,1,*))
    time=timei(*,*)
    local_time=(longi-sub_solar_longitude)*24/360+12
    indLT1=where(local_time lt 0)
    if n_elements(indLT1) gt 1 then local_time(indLT1)=24+local_time(indLT1)
    indLT2=where(local_time ge 24)
    if n_elements(indLT2) gt 1 then local_time(indLT2)=local_time(indLT2)-24

    localtime = strmid(mean(local_time),5,8)
    Ls=strmid(SOLAR_LONGITUDE,5,7)
    Latimean=strmid(mean(lati),5,6)
    Longimean=strmid(mean(longi),5,6)
    year=strmid(time(0),8,4)
    month=strmid(time(1),10,2)
    day=strmid(time(2),10,2)
    hour=strmid(time(3),10,2)
    minit=strmid(time(4),10,2)


    ; 太陽輝度fileをrestore
    restore, path_file + 'specmars.sav'

    ; 青木さんに聞く
    st1 = reform(bt(1d7*jdat(*,251,*)/(1d4/wvl(251))^2d,1d4/wvl(251)))
    st2 = reform(bt(1d7*jdat(*,252,*)/(1d4/wvl(252))^2d,1d4/wvl(252)))
    st3 = reform(bt(1d7*jdat(*,253,*)/(1d4/wvl(253))^2d,1d4/wvl(253)))
    st4 = reform(bt(1d7*jdat(*,254,*)/(1d4/wvl(254))^2d,1d4/wvl(254)))
    st5 = reform(bt(1d7*jdat(*,255,*)/(1d4/wvl(255))^2d,1d4/wvl(255)))

    ; wvl(26) 1.299 um
    albedo = reform(jdat(*,26,*)/specmars(26))

    ; wvl(140) 2.777 um
    rad_2777 = reform(jdat(*,140,*))

    ; wvl(170) 3.401 um, wvl(176) 3.525 um
    ; ice_indexは3.401/3.525 umの比率でみている
    ice_index = reform( (jdat(*,170,*)/specmars(170)) / (jdat(*,176,*)/specmars(176)) )
    
    ; indexを入れるarrayを作成している
    trans = reform(jdat(*,0,*))
    trans_norm = reform(jdat(*,0,*))
    CO = reform(jdat(*,0,*))
    CO2 = reform(jdat(*,0,*))
    O2 = reform(jdat(*,0,*))
    st = reform(jdat(*,0,*))
    LST = reform(jdat(*,0,*))
    R_surface_CO2ice = reform(jdat(*,0,*))

    ; flagの意味を青木さんに聞く
    flag = intarr(io)
    flag(*) = 0

    ; geometry index
    h = 10. ;scale height
    EM = reform(geocube(*, 9, *))*1.e-4
    IA = reform(geocube(*, 8, *))*1.e-4
    alt = reform(geocube(*,12,*)*1.e-3)


    ;wvl(115) 2.54 um, wvl(120) 2.60 um, wvl(123) 2.64 um
    X = [wvl(115), wvl(123)]
    X_CO = [wvl(94), wvl(108)]
    X_CO2 = [wvl(61), wvl(62), wvl(63), wvl(85), wvl(86), wvl(87)]
    CO2_wl = where(findgen(100) ge 61 and findgen(100) le 87 and findgen(100) ne 78)
    X_O2 = [wvl(23), wvl(25)]

    ;CO2 ice on the ground
    ;1.385 um = wvl(32)
    ;1.429 um = wvl(35)
    ;1.443 um = wvl(36)

    for i = 0, io-1 do begin
        for j = 0, ip-1 do begin 
            Y = [jdat(j,115,i), jdat(j,123,i)]
            Y_CO = [jdat(j,94,i), jdat(j,108,i)]
            Y_CO2 = [jdat(j,61,i)/specmars(61), jdat(j,62,i)/specmars(62), jdat(j,63,i)/specmars(63), $
                        jdat(j,85,i)/specmars(85), jdat(j,86,i)/specmars(86), jdat(j,87,i)/specmars(87)]
            Y_O2 = [jdat(j,23,i), jdat(j,25,i)]

            ; 吸収の深さを調べている 
            coef = linfit(X,Y)
            coef_CO = linfit(X_CO,Y_CO)
            coef_CO2 = linfit(X_CO2,Y_CO2)
            coef_O2 = linfit(X_O2, Y_O2)

            cont = coef(0) + coef(1)*wvl(120)
            cont_CO = coef_CO(0) + coef_CO(1)*wvl(100)
            cont_CO2 = coef_CO2(0) + coef_CO2(1)*wvl(CO2_wl)
            cont_O2 = coef(0) + coef(1)*wvl(24)

            CO2(j,i) = total(1.0 - (jdat(j,CO2_wl,i)/specmars(CO2_wl))/cont_CO2,/double,/nan)
            RF1429 = jdat(j,35,i)/specmars(35)*cos(geocube(j,8,i)*1e-4*!DTOR)
            RF1385 = jdat(j,32,i)/specmars(32)*cos(geocube(j,8,i)*1e-4*!DTOR)
            RF1443 = jdat(j,36,i)/specmars(36)*cos(geocube(j,8,i)*1e-4*!DTOR)
            R_surface_CO2ice(j,i) = RF1429/((RF1385^0.5d)*(RF1443^0.5d))

        ;        plot, wvl(CO2_wl), jdat(j,CO2_wl,i)/specmars(CO2_wl),psym=-1,xs=1
        ;        oplot, wvl(CO2_wl), cont_CO2, color=254
        ;        stop
        ;        
            trans(j,i) = 1.0 - jdat(j,120,i)/cont
            CO(j,i) = 1.0 - jdat(j,100,i)/cont_CO
            O2(j,i) = 1.0 - jdat(j,24,i)/cont_O2
        ;        airmass1 = 1.0 / cos(geocube(j,4,i)*1e-4*!DTOR)
        ;        airmass2 = 1.0 / cos(geocube(j,5,i)*1e-4*!DTOR)
            airmass1 = 1.0 / cos(geocube(j,8,i)*1e-4*!DTOR)
            ;airmass1 = 1.0 / (1.d - cos(geocube(j,8,i)*1e-4*!DTOR))
            airmass2 = 1.0 / cos(geocube(j,9,i)*1e-4*!DTOR)
            albedo(j,i) = albedo(j,i); /cos(geocube(j,8,i)*1e-4*!DTOR)
            trans(j,i) = trans(j,i); / (airmass1 + airmass2)
            CO(j,i) = CO(j,i); / (airmass1 + airmass2)
            O2(j,i) = O2(j,i)
            CO2(j,i) = CO2(j,i) ;/ (airmass1 + airmass2)
            trans_norm(j,i) = trans(j,i) * exp(alt(j,i)/H)
            st_sets = [st1(j,i),st2(j,i),st3(j,i),st4(j,i),st5(j,i)]
            st(j,i) = median(st_sets)
            LST(j,i) = (-SUB_SOLAR_LONGITUDE + longi(j,i)) / 15. + 12.
            if EM(j,i) ge 60. then flag(i) = 1
            if IA(j,i) ge 90. then flag(i) = 1
            if JDAT(j,115,i) ge 2. then flag(i) = 1
            if JDAT(j,115,i) le 0.02 then flag(i) = 1
            if alt(j,i) ge 65.535 then flag(i) = 1
            stop
        endfor
    endfor

    ; imageの出力
    loadct, 39
    Set_Plot, 'x'

    ; プロットの枠組み
    window, 2, xs=1250,ys=800
    !P.Multi = [0, 3, 2]
    device, retain=1, decomposed=0

    ; latitude / longitude
    lat1 = lati
    minlati = min(lati)
    maxlati = max(lati)
    lon1 = longi
    minlon = min(longi)
    maxlon = max(longi)

    ; ===== CO plot =======
    min_CO = min(CO(*,*),/nan);0.
    max_CO = max(CO(*,*),/nan);0.05; max(CO(*,where(flag(*) eq 0)),/nan)*1.1
    min_CO_cg = min(CO(*,*),/nan);0.
    max_CO_cg = max(CO(*,*),/nan);0.05; max(CO(*,where(flag(*) eq 0)),/nan)*1.1
    rn_CO = max_CO - min_CO
    if rn_CO eq 0 then max_CO_cg = 1d
    if max_CO gt 10000000000 then max_CO_cg = 1d
    if min_CO lt -10000000000 then min_CO_cg = 0d

    loadct,39
    plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.55,0.31,0.92], /nodata, charsize=3, title='CO'+' Max: '+string(max_CO,format='(e0.1)' ), xtitle='Lon', ytitle='Lat'
    cgLOADCT, 39
    cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[min_CO_cg,max_CO_cg],position=[0.04,0.58,0.046,0.86], TCHARSIZE=10, /vertical 
    loadct, 39

    for i = 0, io-1 do begin
        for j = 0, ip-1 do begin
            color = (CO(j,i)-min_CO)/rn_CO*254.
            if color gt 254 then color = 254
            if color lt 0 then color = 0
            plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
        endfor
    endfor
    
    ; ======= O2 plot ========
    min_O2 = min(O2(*,*),/nan);0.
    max_O2 = max(O2(*,*),/nan);0.05; max(CO(*,where(flag(*) eq 0)),/nan)*1.1
    min_O2_cg = min(O2(*,*),/nan);0.
    max_O2_cg = max(O2(*,*),/nan);0.05; max(CO(*,where(flag(*) eq 0)),/nan)*1.1
    rn_O2 = max_O2 - min_O2
    if rn_O2 eq 0 then max_O2_cg =  min_O2_cg + 1d
    if max_O2 gt 10000000000 then max_O2_cg = 1d
    if min_O2 lt -10000000000 then min_O2_cg = 0d
    
    loadct,39
    plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.55,0.65,0.92], /nodata, charsize=3, title='O2'+' Max: '+string(max_O2,format='(e0.1)' ), xtitle='Lon', ytitle='Lat'
    cgLOADCT, 39
    cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[min_O2_cg,max_O2_cg],position=[0.38,0.58,0.386,0.86], TCHARSIZE=10, /vertical 
    loadct, 39

    for i = 0, io-1 do begin
        for j = 0, ip-1 do begin
            color = (O2(j,i)-min_O2)/rn_O2*254.
            if color gt 254 then color = 254
            if color lt 0 then color = 0
            plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
        endfor
    endfor

    ; ===== albedo plot =====
    max_albedo = max(albedo(*,*),/nan)
    min_albedo = min(albedo(*,*),/nan)
    max_albedo_cg = max(albedo(*,*),/nan)
    min_albedo_cg = min(albedo(*,*),/nan)
    rn_albedo = max_albedo - min_albedo
    if rn_albedo eq 0 then max_albedo_cg =  1d
    if max_albedo gt 10000000000 then max_albedo_cg = 1d
    if min_albedo lt -10000000000 then min_albedo_cg = 0d
    
    loadct,39
    plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.55,0.99,0.92], /nodata, charsize=2, title='Albedo (I/F@1.3um): '+string(min_albedo,format='(f5.3)')+' - '+string(max_albedo,format='(f5.3)'), xtitle='Lon', ytitle='Lat'
    cgLOADCT, 39
    cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[min_albedo_cg,max_albedo_cg],position=[0.72,0.58,0.726,0.86], TCHARSIZE=10, /vertical 

    for i = 0, io-1 do begin
        for j = 0, ip-1 do begin
            color = (albedo(j,i)-min_albedo)/rn_albedo*254.
            if color gt 254 then color = 254
            if color lt 0 then color = 0
            plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
        endfor
    endfor

    ; ==== ice plot =====
    max_ice = max(ice_index(*,*));1.0
    min_ice = min(ice_index(*,*));0.5
    max_ice_cg = max(ice_index(*,*));1.0
    min_ice_cg = min(ice_index(*,*));0.5
    rn_ice = max_ice - min_ice
    if rn_ice eq 0 then max_ice_cg =  1d
    if max_ice gt 10000000000 then max_ice_cg = 1d
    if min_ice lt -10000000000 then min_ice_cg = 0d

    loadct,39
    plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.08,0.31,0.45], /nodata, charsize=3, title='Ice index (3.40/3.52): '+string(min_ice,format='(F3.1)')+' - '+string(max_ice,format='(F3.1)'), xtitle='Lon', ytitle='Lat'
    cgLOADCT, 39
    cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[min_ice_cg,max_ice_cg],position=[0.04,0.12,0.046,0.4], TCHARSIZE=10, /vertical 
    loadct, 39

    for i = 0, io-1 do begin
        for j = 0, ip-1 do begin
            color = (ice_index(j,i)-min_ice)/rn_ice*254.
            if color gt 254 then color = 254
            if color lt 0 then color = 0
            plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
        endfor
    endfor
    

    ; ===== CO2 ice =======
    max_CO2ice = max(R_surface_CO2ice(*,*));1.0
    min_CO2ice = min(R_surface_CO2ice(*,*));0.5
    max_CO2ice_cg = max(R_surface_CO2ice(*,*));1.0
    min_CO2ice_cg = min(R_surface_CO2ice(*,*));0.5
    rn_CO2ice = max_CO2ice - min_CO2ice
    if rn_CO2ice eq 0 then max_CO2ice_cg =  1d
    if max_CO2ice gt 10000000000 then max_CO2ice_cg = 1d
    if min_CO2ice lt -10000000000 then min_CO2ice_cg = 0d

    loadct,39
    plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.08,0.65,0.45], /nodata, charsize=3, title='R_surface_CO2ice: '+string(min_CO2ice,format='(F3.1)')+' - '+string(max_CO2ice,format='(F3.1)'), xtitle='Lon', ytitle='Lat'
    cgLOADCT, 39
    cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[min_CO2ice_cg,max_CO2ice_cg],position=[0.38,0.12,0.386,0.4], TCHARSIZE=10, /vertical 
    loadct, 39

    for i = 0, io-1 do begin
        for j = 0, ip-1 do begin
            color = (R_surface_CO2ice(j,i)-min_CO2ice)/rn_CO2ice*254.
            if color gt 254 then color = 254
            if color lt 0 then color = 0
            plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
        endfor
    endfor

    ; ==== dust plot =====
    max_rad_2777 = max(rad_2777(*,*),/nan)
    min_rad_2777 = min(rad_2777(*,*),/nan)
    max_rad_2777_cg = max(rad_2777(*,*),/nan)
    min_rad_2777_cg = min(rad_2777(*,*),/nan)
    rn_rad_2777 = max_rad_2777 - min_rad_2777
    if rn_rad_2777 eq 0 then max_rad_2777_cg = 1d
    if max_rad_2777 gt 10000000000 then max_rad_2777_cg = 1d
    if min_rad_2777 lt -10000000000 then min_rad_2777_cg = 0d

    loadct,39
    plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.45], /nodata, charsize=2, title='rad@2777 (dust): '+string(min_rad_2777,format='(e0.1)')+' - '+string(max_rad_2777,format='(e0.1)'), xtitle='Lon', ytitle='Lat'
    cgLOADCT, 39
    cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[min_rad_2777_cg,max_rad_2777_cg],position=[0.72,0.12,0.726,0.4], TCHARSIZE=10, /vertical 
    loadct, 39

    for i = 0, io-1 do begin
        for j = 0, ip-1 do begin
            color = (rad_2777(j,i)-min_rad_2777)/rn_rad_2777*254.
            if color gt 254 then color = 254
            if color lt 0 then color = 0
            plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
        endfor
    endfor

    loadct, 39
    xyouts,0.06,0.975,fileorbit,charsize=2.5,color=0,/normal
    xyouts,0.2,0.975,'Time : '+year+'/'+month+'/'+day+' '+hour+':'+minit,color=0,/normal,charsize=2
    xyouts,0.45,0.975,'Ls : '+Ls,color=0,/normal,charsize=2
    xyouts,0.57,0.975,'Local Time : '+localtime,color=0,/normal,charsize=2

  ; =========== save work file and QL ===========
    ; create work file
    save, lati, longi, rad_2777, R_surface_CO2ice, ice_index, albedo, O2, CO, CO2, trans, Ls, year, month, day, hour, minit, localtime, filename = path_work + 'INDEX_' + fileorbit + '.sav'

    snapshot = TVRD(True=1)
    Write_JPEG, path_ql  + fileorbit +'.jpg', snapshot, True=1, Quality=75

    skip1:

    print, 'ORB LOOP: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1
    print, 'ORB number: ', fileorbit

endfor
stop

end