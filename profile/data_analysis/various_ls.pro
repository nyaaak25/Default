pro various_ls

device, decomposed = 0, retain = 2
loadct,39

; === restore work file ===
pathfile = '/data2/omega/sav/'
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

; === MY28 Ls270-300 plot, Look LS variations ===
loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB4495_3.sav') ; MY28 Ls 270
loop_b = fix(loop_b(0))
loop_a = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB4667_5.sav') ; MY28 Ls 300
loop_a = fix(loop_a(0))

; === mcd SPmap ===
;minp1 = 400
;maxp1 = 1000

; === ret SPmap ===
;minp1 = 800
;maxp1 = 1400

; === data plot range ===
;minlon = 290d
;maxlon = 320d

;minlon = 190d
;maxlon = 210d

minlon = 230d
maxlon = 250d

; === プロットの枠組み ===
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

; === plot ===
window, 0, xs=1800, ys=1300
loadct,39
plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[270, 300],  position=[0.16,0.15,0.99,0.9], back=255, color=0, /nodata, charsize=3, title='Retrieved pressure Ls 270-300', xtitle='Ls [deg]', ytitle='Latitude [deg]'
;plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360],  position=[0.16,0.15,0.99,0.9], back=255, color=0, /nodata, charsize=3, title='MCD v5.3 pressure Ls 300-360', xtitle='Longitude [deg]', ytitle='Latitude [deg]'
;plot, findgen(100),xs=1, ys=1, yr=[-90, 90], xr=[0, 360], back=255, color=0, /nodata, charsize=3, title='Ret - mcd pressure Ls 300-360', xtitle='Longitude [deg]', ytitle='Latitude [deg]'

;cgLOADCT, 39
;cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.051,0.12,0.055,0.86], title='Pa', TCHARSIZE=10, /vertical 

; === LOOP START ===
for loop = loop_b, loop_a do begin
    ; fileのrestore
    file = files(loop)
    restore, file

    ; ORB名を取ってくる
    fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))

    ; n_indを格納
    ip = n_elements(LATI(*,0))
    io = n_elements(LATI(0,*))

    lat1 = lati
    lon1 = longi

    ; === Retrieval pressure map ===
    ;p1 = pressure
    ;p1 = exp(p1)

    ; === MCD prrediction pressure map ===
    ;p1 = MCDpressure

    ; === dev ret - mcd ===
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

    LS_file = pathfile + fileorbit + '.sav'
    restore, LS_file
    LS = SOLAR_LONGITUDE

    ls_arr = dblarr(ip,io)
    ls_arr = replicate(ls,ip,io)

    ;for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 1 then p1(i,j) = !VALUES.F_NAN
    ;for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 0 then p1(i,j) = !VALUES.F_NAN
    for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 1 then p1(i,j) = !VALUES.F_NAN

    ; === ORB名を記入したい！ ===
    ; x plot range
    Long_Ls = (300d - 270d) ; Lsデータのカバレージ。Ls270からLs300までの間の分割数
    Long_pix = (0.96d - 0.12d) ; ピクセルデータのカバレージ。 0.96はLs300の軸中心、0.12はLs270の中心
    xran = ((Long_pix/Long_Ls)*(Ls - 270d)) + 0.126d

    ; y plot range
    lat_plot = lat1 + 90d
    for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 1 then lat_plot(i,j) = !VALUES.F_NAN

    maxlat = max(lat_plot)
    Long_lat = (90d + 90d) ; Latデータのカバレージ。Lat -90から90°までの間
    Long_pixel = (0.9d - 0.15d) ; ピクセルデータのカバレージ。0.9はLat90°の中心、0.15はLat-90の中心
    yran = ((Long_pixel/Long_lat)*(maxlat - 0d)) + 0.15d
    
    loadct, 39
    for i = 0, ip-1 do begin 
        for j = 0, io-1 do begin 
            if p1(i,j) gt 0 and lon1(i,j) gt minlon and lon1(i,j) lt maxlon then begin
                plots, ls_arr(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1
                xyouts,xran, yran, fileorbit, charsize=1.5, color=0, /normal
            endif
        endfor
    endfor
    print, 'ORB LOOP: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1

endfor

snapshot = TVRD(True=1)
;Write_JPEG, path_ql + 'EW_MY28_Ls270-300_Lon230-250_ret-mcd.jpg', snapshot, True=1, Quality=75
Write_JPEG, path_ql + 'EW_MY28_Ls270-300_Lon230-250_ret-mcd_inORB.jpg', snapshot, True=1, Quality=75


stop
end