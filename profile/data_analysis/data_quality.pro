
;------------------------------------------------------------------------------------------------------------------------
; +++
; create by Akira Kazama
; 
; data_movie　 ::2023.4.13 Thu 13:43:00
; 死んでいる素子を判断するためにplotをする
; +++
;------------------------------------------------------------------------------------------------------------------------


pro data_quality

device, decomposed = 0, retain = 2
loadct,39

; .sav fileをdirectoryに入っているfileをすべて持ってくる
pathfile = '/data2/omega/sav/'
files=file_search(pathfile + '*.sav',count=count)

; MY27のdataをみる
;loop_b = where(files eq '/data2/omega/sav/ORB0920_3.sav') ;MY27はじまりのORB
;loop_a = where(files eq '/data2/omega/sav/ORB2595_0.sav') ; MY27おわりのORB

; MY28のdataをみる
;loop_b = where(files eq '/data2/omega/sav/ORB2607_0.sav') ;MY28はじまりのORB
;loop_a = where(files eq '/data2/omega/sav/ORB5055_5.sav') ; MY28おわりのORB

; MY29
;loop_b = where(files eq '/data2/omega/sav/ORB5062_0.sav') ; MY29はじまりのORB
;loop_b = where(files eq '/data2/omega/sav/ORB7410_3.sav')
loop_b = where(files eq '/data2/omega/sav/ORB0920_3.sav')
loop_a = where(files eq '/data2/omega/sav/ORB7454_4.sav') ; MY29おわりのORB


loop_b = fix(loop_b(0))
loop_a = fix(loop_a(0))

; *** memo ***
; 最初から死んでいる素子は2つ
; 死んでいる素子は、 wvl(17), wvl(27)
; (17) 2.0392201, (27) 2.1776299 μm
;
; ORB2432_1から死んでいる素子が3つに増加
; 死んでいる素子は、wvl(11), wvl(17), wvl(27)
; (11) 1.9553300, (17) 2.0392201, (27) 2.1776299 μm
;
; ORB2618_0から死んでいる素子が4つに増加！
; 死んでいる素子は、wvl(11), wvl(17), wvl(19), wvl(27)
; (11) 1.9553300, (17) 2.0392201, (19) 2.0670500, (27) 2.1776299 μm

; ORB5306_0から死んでいる素子が5つに増加！
; 死んでいる素子は、wvl(1), wvl(11), wvl(17), wvl(19), wvl(27)
; (1) 1.8143300, (11) 1.9553300, (17) 2.0392201, (19) 2.0670500, (27) 2.1776299 μm
;
; ORB5408_0から死んでいる素子が6つに増加！
; 死んでいる素子は、wvl(1), wvl(8), wvl(11), wvl(17), wvl(19), wvl(27)
; (1) 1.8143300, (8) 1.9131700, (11) 1.9553300, (17) 2.0392201, (19) 2.0670500, (27) 2.1776299 μm

; ORB5692_0から死んでいる素子が7つに増加！
; 死んでいる素子は、wvl(1), wvl(8), wvl(7), wvl(11), wvl(17), wvl(19), wvl(27)
; (1) 1.8143300, (7) 1.8990901,(8) 1.9131700, (11) 1.9553300, (17) 2.0392201, (19) 2.0670500, (27) 2.1776299 μm

; ORB6008_0から死んでいる素子が9つに増加！
; 死んでいる素子は、wvl(1), wvl(2), wvl(6), wvl(7), wvl(8), wvl(11), wvl(17), wvl(19), wvl(27)
; (1) 1.8143300, (2) 1.8284900, (6) 1.8850000, (7) 1.8990901,(8) 1.9131700, (11) 1.9553300, (17) 2.0392201, (19) 2.0670500, (27) 2.1776299 μm



; plotしてみる
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

;for loop = 0, count-1 do begin
for loop = loop_b, loop_a do begin
    file = files(loop)

    restore, file
    ;restore, '/data2/omega/sav/ORB5055_5.sav'

    fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))

    ip = n_elements(LATI(*,0))
    io = n_elements(LATI(0,*))

    ; CO2 absorption
    CO2=where(wvl gt 1.8 and wvl lt 2.2) ; 波長範囲はここで変更可能
    wvl=wvl[CO2]
    jdat = jdat(*,CO2,*)

    index = ip * io

    notgood_ind = where(jdat le 0.00000001)
    jdat(notgood_ind) = -0d/0d
    n_ind = n_elements(notgood_ind)

    band=where(wvl gt 1.94 and wvl lt 2.09)

    bad_ind = n_ind/index
    cont = [0,1,2, 23, 25, 26]
    ;LUT7 = [1,2,5,6,7,8,11,17,18,19,24,27]
    ;LUT7 = [1,5,7,8,11,17,18,19,24,27]
    LUT7 = [5,8,17,18,19,24,27]
    jdat_7 = jdat(0,*,0)
    jdat_7(LUT7) = -0d/0d

    jdat_cont = jdat(0,cont,0)
    jdat_band = jdat(0,band,0)
    
    window, 2
    plot, wvl, jdat(0,*,0), back=255, color=0;, psym = 2, symsize = 3
    oplot, wvl(cont), jdat_cont(*),color=60, psym = 1, symsize = 2
    oplot, wvl(band), jdat_band(*),color=150, psym = 1, symsize = 2
    oplot, wvl, jdat_7(*),color=254, psym = 1

    stop


    ; search for broken spectrel [下参照]
    ; ORBの各indexに対してbad pixel(bad data)を探す
    ;window, 0, xs=1000, ys=800
    ;loadct,39
    ;plot, findgen(100),xr=[0,ip-1],yr=[0,30],ys=1, xs=1, back=255, color=0, /nodata, charsize=3, title='n_ind plot', xtitle = 'cube index'

    ;nan = where(jdat le 0.000001)
    ;for i = 0, ip-1 do begin
    ;    notgood_ind = where(jdat(i,*,0) le 0.000001)
        ;sza = reform(geocube(i,10,0))*1.e-4
    ;    n_ind = n_elements(notgood_ind)

    ;    loadct, 39
    ;    plots, i, n_ind, color=60, psym=2,symsize=2
    ;    print, i
    ;    print, n_ind
    ;endfor
    ;stop
    

    if bad_ind gt 9 then print, 'file orbit: ' + fileorbit + ',   Bad pixel: ', bad_ind
    ;if bad_ind lt 4 then stop 
    if bad_ind eq 9 then print, 'file orbit: ' + fileorbit + ',   Hello :) ' 
    ;if bad_ind eq 4 then stop 
    ;if bad_ind gt 4 then print, 'file orbit: ' + fileorbit + ',   Bad pixel: ', bad_ind 


endfor



stop
end
