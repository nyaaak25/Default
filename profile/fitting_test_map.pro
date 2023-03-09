
; 各地点でのfitting具合を表示させるプロット

pro fitting_test_map

Set_Plot, 'x'
device, retain=1, decomposed=0
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

; restore file
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_albedo_update.sav'
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_B2_update.sav'
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0931_3_B1.sav'
fit = bestfit
cont = continuum

dust = dustmap
albedo = albedomap
pressure = pressure

;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0931_3.sav'
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0920_3.sav'

;Block 1
;ind = where_xyz(longi ge 274 and longi le 277 and lati ge 56 and lati le 58, xind=xind, yind=yind)
;
; Block 2
ind = where_xyz(longi ge 274 and longi le 277 and lati ge 54 and lati le 56, xind=xind, yind=yind)

; Block 3
;ind = where_xyz(longi ge 274 and longi le 277 and lati ge 52 and lati le 54, xind=xind, yind=yind)


; count number of array
ip = n_elements(LATI(*,0))
io = n_elements(LATI(0,*))
ip_b = min(xind)
ip_a = max(xind)
io_b = min(yind)
io_a = max(yind)

; CO2 absorption line selection
CO2=where(wvl ge 1.8 and wvl le 2.2)
wvl_CO2=wvl[CO2]
jdat_CO2=jdat(*,CO2,*)

; I/F inversion
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/specmars.sav'
specmars = specmars
specmars_CO2 = specmars[CO2]

; not use spectrum selection
jdat_CO2(*,0:3,*)= !VALUES.F_NAN
jdat_CO2(*,8,*)= !VALUES.F_NAN
jdat_CO2(*,17,*)= !VALUES.F_NAN
jdat_CO2(*,24,*)= !VALUES.F_NAN
jdat_CO2(*,27,*)= !VALUES.F_NAN

; create array
nx = n_elements(CO2)
sumchi = dblarr(ip,io)
m_o = dblarr(ip,io,21)
y_all = dblarr(ip,io,21)
con_all = dblarr(ip,io,21)
err = dblarr(nx)


for l = ip_b, ip_a -1 do begin ;loop for slit scan
  for k = io_b, io_a -1 do begin ;test getting surface feature
    x = wvl_CO2
    y_obs = reform(jdat_CO2(l, *, k))
    y_obs = y_obs/specmars_CO2
   
    f = reform(fit(l,k,*))
    f(0:3) = -0d/0d
    f(8) = -0d/0d
    f(17) = -0d/0d
    f(24) = -0d/0d
    f(27) = -0d/0d
    
    con = reform(cont(l,k,*))
    con(0:3) = -0d/0d
    con(8) = -0d/0d
    con(17) = -0d/0d
    con(24) = -0d/0d
    con(27) = -0d/0d
    
    good = where(FINITE(y_obs) eq 1)
    y_all(l,k,*) = y_obs(good)
    m_o(l,k,*) =  y_obs(good) - f(good)
    con_all(l,k,*) = con(good)
       
    ;chi-square
    ; 残差2乗和
    N = n_elements(f(good))
    rss = (y_obs(good) - f(good))^2
    err(*) = median(y_obs)*1d-3
    chi_in = total((rss/err)^(0.5d))
    chi = chi_in*(N^(0.5d))
    
    sumchi(l,k) = chi

    
    ; R^2
;    rss = (y_obs(good) - f(good))^2
;    sum_rss = total(rss)    
;    ave_y = mean(y_obs(good))
;    tss = (y_obs(good) - ave_y)^2
;    sum_tss = total(tss)
;    
;    sumchi(l,k) = 1 - (sum_rss/sum_tss)
    

;    1地点1地点のfitting具合を表示する         
;    window,5, xs=500,ys=500
;    plot, x(good), f(good) - y_obs(good), back=255, color=0, thick=3, psym=-1, xr=[1.85, 2.2], xs=1
;    ;plot, x(good), y_obs(good), yr=[-0.1,0.4], back=255, color=0, thick=3, psym=-1, xr=[1.85, 2.2], xs=1
;    plot, x(good), y_obs(good),  yr=[-0.05,0.3], back=255, color=0, thick=3, psym=-1, xr=[1.85, 2.2], xs=1
;    oplot, x(good), f(good), color=254, thick=2, psym=-1, linestyle=2
;    oplot, x(good), f(good) - y_obs(good), color=0, thick=2, psym=-1, linestyle=2
;    xyouts, 2.05, -0.08, 'SP='+strcompress(pressure(l,k)), charsize=1.5, color=0
;    xyouts, 2.05, -0.06, 'Albedo='+strcompress(albedo(l,k)), charsize=1.5, color=0
;    xyouts, 2.05, -0.04, 'Dust='+strcompress(dust(l,k)), charsize=1.5, color=0
       
  endfor
endfor

; color bar表示に使用するためのパラメータ
min_chi = min(sumchi(ip_b:ip_a-1,io_b:io_a-1))
max_chi = max(sumchi(ip_b:ip_a-1,io_b:io_a-1))
mean_chi = max_chi - min_chi

loadct, 39

; fitting具合の図をプロット
;window,7, xs=500,ys=500
;plot, findgen(10), yr=[53.5, 56.5], xr=[273, 277.5], back=255, color=0, /nodata, charsize=2
;for i = 0, 127 do for j = 0, 595 do if sumchi(i,j) gt 0 then plots, longi(i,j), lati(i,j), color=(sumchi(i,j)-min_chi)/mean_chi * 254., psym=2

; スペクトルの残差を重ねてプロット
window,5, xs=500,ys=500
plot, findgen(10), yr=[0.9, 1.1], xr=[1.85, 2.2], back=255, color=0, /nodata, charsize=2, xs=1
for i = 0, 127 do for j = 0, 595 do if y_all(i,j,0) gt 0 then oplot, x(good), m_o(i,j,*), color=254, thick=1, linestyle=0
for i = 0, 127 do for j = 0, 595 do if y_all(i,j,0) gt 0 then oplot, x(good), y_all(i,j,*), color=0, thick=1,linestyle=0
for i = 0, 127 do for j = 0, 595 do if y_all(i,j,0) gt 0 then oplot, x(good), con_all(i,j,*), color=60, thick=1,linestyle=0
;for i = 0, 127 do for j = 0, 595 do if m_o(i,j,0) gt 0 then oplot, x(good), m_o(i,j,*) / median(y_all(i,j,*)), color=254, thick=1,linestyle=0


stop
end