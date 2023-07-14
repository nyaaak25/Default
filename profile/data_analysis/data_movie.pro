
;------------------------------------------------------------------------------------------------------------------------
; +++
; create by Akira Kazama
; 
; data_movie　 ::2023.4.13 Thu 13:00:00
; リトリーバルされた圧力をパラパラ漫画風にplotするためのプログラム
; +++
;------------------------------------------------------------------------------------------------------------------------


;------------------------------------------------------------------------------------------------------------------------
function Ls_search, Ls
;------------------------------------------------------------------------------------------------------------------------
; Lsを探すためのfunction

path_work = '/work1/LUT/SP/EWwork/'
pathfile = '/data2/omega/sav/'

files=file_search(path_work + '*.sav',count=count)

if Ls le 90 then loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB2607_0.sav') ; MY28 Ls0
if Ls gt 90 and Ls le 180 then loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB3316_1.sav') ; MY28 Ls90
if Ls gt 180 and Ls le 270 then loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB3843_3.sav') ; MY28 Ls180
if Ls gt 270 then loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB4495_3.sav') ; MY28 Ls270

loop_b = fix(loop_b(0))
loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB5055_5.sav') ; MY28おわりのORB
loop_a = fix(loop_a(0))

for loop = loop_b, loop_a do begin
    file = files(loop)
    restore, file

    fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))

    LS_file = pathfile + fileorbit + '.sav'
    restore, Ls_file
    LS_ind = SOLAR_LONGITUDE
    
    dev = LS_ind - Ls

    if (dev gt 0.0d0) then ind = loop
    if (dev gt 0.0d0) then goto, skip1

endfor
skip1:

file_name = files(ind)
fileorbit = strupcase(strmid(file_name,12,9,/REVERSE_OFFSET))

return, fileorbit
end


pro data_movie

device, decomposed = 0, retain = 2
loadct,39

; restore work file
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_test/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

; +++ ここでplot設定 ++++   !!! 基本的にパラメーターをいじる場所はここだけ !!!
; 始まりと終わりのLsを指定する
start_ls = 0
end_ls = 360
; 何度刻みのパラパラ漫画を作りたいかを決める
grid = 5
; ++++++++++++++++++++


endloop = (end_ls - (start_ls + grid)) / grid

; file loop start
for n = 0, endloop do begin

    ; 任意のLsのファイルを持ってくる
    Lsb = start_ls + grid*n ; 何度刻みのパラパラ漫画を作成したいかをnの前の数字で指定
    file_ind_before = Ls_search(Lsb)

    Lsa = (start_ls + grid) + grid*n ; 何度刻みのパラパラ漫画を作成したいかをnの前の数字で指定
    file_ind_after = Ls_search(Lsa)

    loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_' + file_ind_before + '.sav') 
    loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_' + file_ind_after + '.sav')
    loop_b = fix(loop_b(0))
    loop_a = fix(loop_a(0))

    ; プロットの枠組み
    A = FINDGEN(17) * (!PI*2/16.)
    USERSYM, COS(A), SIN(A), /FILL

    ;; plot
    ;window, 1, xs=2600, ys=600
    window, 25, xs=4000, ys=800 ; connect display 
    !P.Multi = [0, 2, 2]
    loadct,39

    ; === ret SP plot === 
    minp1 = 400
    maxp1 = 1000

    plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360],  position=[0.1,0.6,0.49,0.9], back=255, color=0, /nodata, charsize=2, title='Retrieved pressure', xtitle='Longitude [deg]', ytitle='Latitude [deg]'
    xyouts,0.08,0.97,'Ls = ' +strmid(Lsb,5,6) + '  to  ' +strmid(Lsa,5,6),charsize=2.5,color=0,/normal
    cgLOADCT, 39
    cgColorbar, Divisions=4, Minor=5, charsize=1, Range=[minp1,maxp1],position=[0.04,0.62,0.044,0.86], title='Pa', TCHARSIZE=1, /vertical

    for loop = loop_b, loop_a do begin
        ; fileのrestore
        file = files(loop)
        restore, file

        ; n_indを格納
        ip = n_elements(LATI(*,0))
        io = n_elements(LATI(0,*))

        lat1 = lati
        lon1 = longi

        ; === Retrieval pressure map ===
        p1 = pressure
        p1 = exp(p1)

        for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 1 then p1(i,j) = !VALUES.F_NAN

            loadct, 39
        for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1
        print, 'ORB LOOP: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1
    
    endfor

    ; === dev ret - mcd ===
    plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360],  position=[0.1,0.15,0.49,0.45], back=255, color=0, /nodata, charsize=2, title='retSP - MCD', xtitle='Longitude [deg]', ytitle='Latitude [deg]'

    for loop = loop_b, loop_a do begin
        ; fileのrestore
        file = files(loop)
        restore, file

        ; n_indを格納
        ip = n_elements(LATI(*,0))
        io = n_elements(LATI(0,*))

        lat1 = lati
        lon1 = longi

        p2 = pressure
        p2 = exp(p2)
        medip2 = median(p2)
        p3 = MCDpressure

        p1 = p2 - p3
        medip1 = median(p1)

        ; 何%変動までみたいかによってcolor rangeを変更
        dev_range = 0.01*medip2
        minp1 = medip1 - dev_range
        maxp1 = medip1 + dev_range
        
        for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 1 then p1(i,j) = !VALUES.F_NAN
    
        loadct, 39
        for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1
        print, 'ORB LOOP: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1

    endfor

    ; === MCD SP plot ===
    minp2 = 400
    maxp2 = 1000

    plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360],  position=[0.61,0.6,0.99,0.9], back=255, color=0, /nodata, charsize=2, title='MCD v5.3 pressure', xtitle='Longitude [deg]', ytitle='Latitude [deg]'
    cgLOADCT, 39
    cgColorbar, Divisions=4, Minor=5, charsize=1, Range=[minp2,maxp2],position=[0.55,0.62,0.554,0.86], title='Pa', TCHARSIZE=1, /vertical

    for loop = loop_b, loop_a do begin
        ; fileのrestore
        file = files(loop)
        restore, file

        ; n_indを格納
        ip = n_elements(LATI(*,0))
        io = n_elements(LATI(0,*))

        lat1 = lati
        lon1 = longi
    
        ; === MCD prrediction pressure map ===
        p1 = MCDpressure

        p2 = pressure
        for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN

        loadct, 39
        for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp2)/(maxp2-minp2)*254., psym=8, symsize=1
        print, 'ORB LOOP: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1

    
    endfor

    ; === dust opacity plot ===
    minp3 = 0d
    maxp3 = 0.25d

    plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360],  position=[0.61,0.15,0.99,0.45], back=255, color=0, /nodata, charsize=2, title='Dust opacity', xtitle='Longitude [deg]', ytitle='Latitude [deg]'
    cgLOADCT, 39
    cgColorbar, Divisions=4, Minor=5, charsize=1, Range=[minp3,maxp3],position=[0.55,0.17,0.554,0.41], title='Dust opacity', TCHARSIZE=1, /vertical

    for loop = loop_b, loop_a do begin
        ; fileのrestore
        file = files(loop)
        restore, file

        ; n_indを格納
        ip = n_elements(LATI(*,0))
        io = n_elements(LATI(0,*))

        lat1 = lati
        lon1 = longi

        ; === Retrieval pressure map ===
        p2 = pressure

        ; === dust opacity ===
        p1 = dustmap


        for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN

        loadct, 39
        for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp3)/(maxp3-minp3)*254., psym=8, symsize=1
        print, 'ORB LOOP: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1

    endfor

    snapshot = TVRD(True=1)
    Write_JPEG, path_ql + '52_EW_MY28_Ls'+strmid(Lsb,5,6)+ '-' +strmid(Lsa,5,6)+ '.jpg', snapshot, True=1, Quality=75

endfor

stop
end