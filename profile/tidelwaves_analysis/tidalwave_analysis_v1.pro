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
;
;------------------------------------------------------------------------------------------------------------------------


pro tidalwave_analysis_v1

device, decomposed = 0, retain = 2

; =============== restore working file ===============
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL_datacover/'

files=file_search(path_work + '*.sav',count=count) ; sav fileすべてを持ってくる

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
loop_b = where(files eq '/work1/LUT/SP/EWwork/EW_work_ORB5062_0.sav') ; MY29 Ls 0
loop_b = fix(loop_b(0))
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
window, 0, xs=1500, ys=1000
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
Write_JPEG, path_ql + 'MY29_ret-mcd_70-30.jpg', snapshot, True=1, Quality=75


stop

end