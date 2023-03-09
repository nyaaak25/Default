pro sys_err_calc

; 系統誤差を引くFactor計算するためのプログラム

Set_Plot, 'x'
device, retain=1, decomposed=0
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL

; ==== restore file ====
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_albedo_update.sav'
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0931_3_albedo_update.sav'
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0920_3_B1.sav'
;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/SPmap_ORB0931_3_B1.sav'
fit = bestfit
dust = dustmap
albedo = albedomap
pressure = pressure

;restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0920_3.sav'
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0931_3.sav'

;Block 1
;ind = where_xyz(longi ge 274 and longi le 277 and lati ge 56 and lati le 58, xind=xind, yind=yind)
;
; Block 2
ind = where_xyz(longi ge 274 and longi le 277 and lati ge 54 and lati le 56, xind=xind, yind=yind)

; Block 3
;ind = where_xyz(longi ge 274 and longi le 277 and lati ge 52 and lati le 54, xind=xind, yind=yind)

; === create sys err factor ===
; syserr_*** : Model - OBSの残差をmedian(OBS)で割ったもの
; syserr_***_dif : Model - OBSの残差をOBSで割ったもの
; syserr_***_dif_2 : Model - OBS
; syserr_***_dif_3 : Model / OBS
 
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/syserr_ORB0920_3_B2_dif_2.sav'
B2_920_sys = sys_err
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/syserr_ORB0931_3_B2_dif_2.sav'
B2_931_sys = sys_err
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/syserr_ORB0920_3_B1_dif_2.sav'
B1_920_sys = sys_err
restore, '/Users/nyonn/IDLWorkspace/Default/savfile/syserr_ORB0931_3_B1_dif_2.sav'
B1_931_sys = sys_err

; average
sys_err_factor = ( B1_931_sys+ B1_920_sys + B2_931_sys + B2_920_sys ) / 4d
print, sys_err_factor
stop
save, sys_err_factor, filename = '/Users/nyonn/IDLWorkspace/Default/savfile/sys_err_factor.sav'
stop

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
m_o_after = dblarr(ip,io,21)
y_all = dblarr(ip,io,21)
err = dblarr(nx)

; jdatとかを計算する
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

    good = where(FINITE(y_obs) eq 1)
    y_all(l,k,*) = y_obs(good)

    m_o_after(l,k,*) =  (y_obs(good) - sys_err_factor) - f(good)
    m_o(l,k,*) =  y_obs(good) - f(good)


  endfor
endfor

; 系統誤差を波長方向に計算させている
;n = 0
;for i = 0, 127 do for j = 0, 595 do if m_o(i,j,0) gt 0 then n = n+1
;
;sys_err =dblarr(21, n)
;n = 0
;
;for i = 0, 127 do begin
;  for j = 0, 595 do begin
;    if m_o(i,j,0) gt 0 then begin
;    sys_err(*,n) = m_o(i, j, *) ;/ y_all(i, j, *) ;/ median(y_all(i, j, *))
;    n = n + 1
;    endif
;  endfor
;endfor

;sys_err = mean(sys_err,dimension=2, /double)
;
;loadct, 39
;window,6, xs=500,ys=500
;plot, findgen(10), yr=[0.8, 1.1], xr=[1.85, 2.2], back=255, color=0, /nodata, charsize=2, xs=1
;for i = 0, 127 do for j = 0, 595 do if y_all(i,j,0) gt 0 then oplot, x(good), y_all(i,j,*), color=0, thick=1,linestyle=0
;for i = 0, 127 do for j = 0, 595 do if y_all(i,j,0) gt 0  then oplot, x(good), m_o(i,j,*), color=254, thick=1,linestyle=0
;oplot, x(good),sys_err, color=60, thick=2,linestyle=0
;stop

;save, sys_err, filename= '/Users/nyonn/IDLWorkspace/Default/savfile/syserr_ORB0931_3_B2_dif_3.sav'
;stop

; === sys err plot ===
;loadct, 39
;window,1, xs=500,ys=500
;plot, findgen(10), yr=[0.8, 1.1], xr=[1.85, 2.2], back=255, color=0, /nodata, charsize=2, xs=1
;oplot, x(good), B1_931_sys, color=0, thick=2,linestyle=0
;oplot, x(good), B1_920_sys, color=254, thick=2,linestyle=0
;oplot, x(good), B2_931_sys, color=60, thick=2,linestyle=0
;oplot, x(good), B2_920_sys, color=100, thick=2,linestyle=0



;;
loadct, 39
window,5, xs=500,ys=500
plot, findgen(10), yr=[-0.03, 0.03], xr=[1.85, 2.2], back=255, color=0, /nodata, charsize=2, xs=1
;for i = 0, 127 do for j = 0, 595 do if y_all(i,j,0) gt 0 then oplot, x(good), y_all(i,j,*), color=0, thick=1,linestyle=0
for i = 0, 127 do for j = 0, 595 do if y_all(i,j,0) gt 0  then oplot, x(good), m_o(i,j,*), color=60, thick=1,linestyle=0
for i = 0, 127 do for j = 0, 595 do if y_all(i,j,0) gt 0  then oplot, x(good), m_o_after(i,j,*), color=254, thick=1,linestyle=0

;oplot, x(good), sys_err_factor, color=60, thick=1,linestyle=0



stop
end