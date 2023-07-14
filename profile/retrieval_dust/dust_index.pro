;------------------------------------------------------------------------------------------------------------------------
;
; ダストのインデックスを作成するためのプログラム
; create by Akira Kazama
; 
; dust_index　 ::2023.7.13 Thu 15:55:00
; 緯度平均を取って、ダスト量を導出するためのプログラム
;
;
;------------------------------------------------------------------------------------------------------------------------


pro dust_index

file = '/data2/omega/sav/ORB0923_3.sav'
restore, file

fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))

A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

pathfile = '/data2/omega/sav/'
restore, pathfile + 'specmars.sav'
specmars = specmars

; cal dust indexを使用して作ったdust indexをここでrestoreをする
path = '/work1/LUT/dust/table/output/'
file_index = path + 'dust_index_cal.sav'
restore, file_index
dust_index_file = dust_index_file

;0:127がSWIR
;128:255がLWIR
;256:352がVNIR

; 緯度1度ごともってくる
min_lati = fix(min(lati))
max_lati = fix(max(lati)) + 1
dev = 1
def = max_lati - min_lati

dust_IF = dblarr(def)
dust_CD = dblarr(def)

dust_grid = dblarr(31)
for j = 0, 30 do dust_grid(j) = j * 0.01d

for loop = 0, def -1 do begin
    ; 緯度1度ごとに平均をもってくる
    ind = where_xyz(lati gt min_lati + loop and lati lt min_lati + loop + dev, xind=xind, yind=yind)

    ; wvl(140): 2.777 μm
    jdat_ind = jdat(xind, 140, yind)
    xp = n_elements(xind)
    yp = n_elements(yind)
    jdat_ind1 = dblarr(xp,yp)

    for i = 0, xp-1 do for j = 0, yp-1 do jdat_ind1(i,j) = jdat_ind(i,0,j) / specmars(140)

    pm = 1 ; 何σに入るかをここで定義する
    std = stddev(jdat_ind1,/nan)
    medi = median(jdat_ind1)

    min_range = medi - (pm * std)
    max_range = medi + (pm * std)

    ; 1σに入るものをいれよう
    for i = 0, xp-1 do begin
        for j = 0, yp-1 do begin 
            if jdat_ind1(i,j) le min_range then jdat_ind1(i,j) = !VALUES.F_NAN
            if jdat_ind1(i,j) ge max_range then jdat_ind1(i,j) = !VALUES.F_NAN
        endfor
    endfor

    dust_mean = mean(jdat_ind1,/nan)
    dust_IF(loop) = dust_mean

    dust_CD(loop) = interpol(dust_grid, dust_index_file, dust_mean)

endfor

path_sav = '/work1/LUT/dust/table/index/'
save, dust_CD, filename = path_sav + fileorbit + '.sav'

stop

end



