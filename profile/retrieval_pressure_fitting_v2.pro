;------------------------------------------------------------------------------------------------------------------------
;
; ある1地点でのLUT + fitting をするプログラム
;
; create by Akira Kazama
; 
; retrieval_pressure_fitting_v1　 ::2022.11.16 Mon 10:52:00
; ver1: LUT update前に対応するプログラム
;
; retrieval_pressure_fitting_v2　 ::2022.11.16 Mon 11:52:00
; ver2: LUT update後に対応するプログラム
; 現状、きちんとfittingされるのはsimple continuumの場合のみ
; 
; +++
;------------------------------------------------------------------------------------------------------------------------


;---------------------------
function forward, x, p
;---------------------------
;x: wavenumber
;p: 9 paprameters, surface pressure, T1, T2, SZA, EA, PhaseA, Dust, Water ice clouds, Albedo
;return: radiance
;---------------------------

nx = n_elements(x)
y = dblarr(nx)
y_tmp = dblarr(15, nx)

;define index
T1 = P(1)
T2 = P(2)
SZA = P(3)
EA = P(4)
PA = P(5)
Dust = P(6)
Waterice = P(7)
Albedo = P(8)

if T1 gt 285. then begin
  print,'Warning: T1 > limit'
  stop
endif
if T1 lt 135. then begin
  print,'Warning: T1 < limit'
  stop
endif
if T2 gt 200. then begin
  print,'Warning: T2 > limit'
  stop
endif
if T2 lt 80. then begin
  print,'Warning: T2 < limit'
  stop
endif
if SZA gt 75. then begin
  print,'Warning: SZA > limit'
  stop
endif
if EA gt 50. then begin
  print,'Warning: EA > limit'
  stop
endif
if PA gt 180. then begin
  print,'Warning: PA > limit'
  stop
endif
if Dust gt 1.5 then begin
  print,'Warning: DUST > limit'
  stop
endif
if Waterice gt 1.0 then begin
  print,'Warning: ICE > limit'
  stop
endif
if Albedo gt 0.6 then begin
  print,'Warning: Albedo > limit'
  stop
endif
if Albedo lt 0.05 then begin
  print,'Warning: Albedo < limit'
  stop
endif

x1a = dblarr(5) ; T1
x1a(0) = (135.d)^1.5d
x1a(1) = (160.d)^1.5d
x1a(2) = (213.d)^1.5d
x1a(3) = (260.d)^1.5d
x1a(4) = (285.d)^1.5d

x2a = dblarr(3) ; T2
x2a(0) = (80.d)^1.5d
x2a(1) = (146.d)^1.5d
x2a(2) = (200.d)^1.5d

x3a = dblarr(6) ; SZA
x3a(0) = exp(-cos(0d/180d*!dpi))
x3a(1) = exp(-cos(15d/180d*!dpi))
x3a(2) = exp(-cos(30d/180d*!dpi))
x3a(3) = exp(-cos(45d/180d*!dpi))
x3a(4) = exp(-cos(60d/180d*!dpi))
x3a(5) = exp(-cos(75d/180d*!dpi))

x4a = dblarr(5) ; EA
x4a(0) = exp(-cos(0d/180d*!dpi))
x4a(1) = exp(-cos(5d/180d*!dpi))
x4a(2) = exp(-cos(10d/180d*!dpi))
x4a(3) = exp(-cos(30d/180d*!dpi))
x4a(4) = exp(-cos(50d/180d*!dpi))

x5a = dblarr(5) ; PA
x5a(0) = -cos(0d/180d*!dpi) + 1.d
x5a(1) = -cos(45d/180d*!dpi) + 1.d
x5a(2) = -cos(90d/180d*!dpi) + 1.d
x5a(3) = -cos(135d/180d*!dpi) + 1.d
x5a(4) = -cos(180d/180d*!dpi) + 1.d

x6a = dblarr(6) ; Dust
x6a(0) = 0.0d
x6a(1) = 0.3d
x6a(2) = 0.6d
x6a(3) = 0.9d
x6a(4) = 1.2d
x6a(5) = 1.5d

x7a = dblarr(3) ; Water-ice
x7a(0) = 0.0d
x7a(1) = 0.5d
x7a(2) = 1.0d

x8a = dblarr(7) ; Surface Albedo
x8a(0) = 0.05d
x8a(1) = 0.1d
x8a(2) = 0.2d
x8a(3) = 0.3d
x8a(4) = 0.4d
x8a(5) = 0.5d
x8a(6) = 0.6d

X1 = (T1)^1.5d ;T1 as test
X2 = (T2)^1.5d ;T2 as test
x3 = exp(-cos(SZA/180d*!dpi))
x4 = exp(-cos(EA/180d*!dpi))
x5 = -cos(PA/180d*!dpi) + 1.d
x6 = Dust ;reduce_dust(loop)*1.86
x7 = Waterice; reduce_waterice(loop)*3.87
x8 = Albedo

for I = 0, 5-1 do begin
  F = x1a(I) - X1
  if (F gt 0.0d0) then j1 = I
  if (F gt 0.0d0) then goto, skip1
endfor
skip1:

for I = 0, 3-1 do begin
  F = x2a(I) - X2
  if (F gt 0.0d0) then j2 = I
  if (F gt 0.0d0) then goto, skip2
endfor
skip2:

for I = 0, 6-1 do begin
  F = x3a(I) - X3
  if (F gt 0.0d0) then j3 = I
  if (F gt 0.0d0) then goto, skip3
endfor
skip3:

for I = 0, 5-1 do begin  ; !!TBD
  F = x4a(I) - X4
  if (F gt 0.0d0) then j4 = I
  if (F gt 0.0d0) then goto, skip4
endfor
skip4:

for I = 0, 5-1 do begin
  F = x5a(I) - X5
  if (F gt 0.0d0) then j5 = I
  if (F gt 0.0d0) then goto, skip5
endfor
skip5:

for I = 0, 6-1 do begin
  F = x6a(I) - X6
  if (F gt 0.0d0) then j6 = I
  if (F gt 0.0d0) then goto, skip6
endfor
skip6:

for I = 0, 3-1 do begin
  F = x7a(I) - X7
  if (F gt 0.0d0) then j7 = I
  if (F gt 0.0d0) then goto, skip7
endfor
skip7:

for I = 0, 7-1 do begin
  F = x8a(I) - X8
  if (F gt 0.0d0) then j8 = I
  if (F gt 0.0d0) then goto, skip8
endfor
skip8:

j1 = j1 -1
j2 = j2 -1
j3 = j3 -1
j4 = j4 -1
j5 = j5 -1
j6 = j6 -1
j7 = j7 -1
j8 = j8 -1

if j1 lt 0 then stop
if j2 lt 0 then stop
if j3 lt 0 then stop
if j4 lt 0 then stop
if j5 lt 0 then stop
if j6 lt 0 then stop
if j7 lt 0 then stop
if j8 lt 0 then stop
if j1+1 ge 5 then stop
if j2+1 ge 3 then stop
if j3+1 ge 6 then stop
if j4+1 ge 5 then stop
if j5+1 ge 5 then stop
if j6+1 ge 6 then stop
if j7+1 ge 3 then stop
if j8+1 ge 7 then stop

;mutli-dimensional interpolations
for i = 0, 15-1 do begin
  for j = 0, nx-1 do begin
    
    if j eq 0 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad0.sav'
    if j eq 1 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad1.sav'
    if j eq 2 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad2.sav'
    if j eq 3 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad3.sav'
    if j eq 4 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad4.sav'
    if j eq 5 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad5.sav'
    if j eq 6 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad6.sav'
    if j eq 7 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad7.sav'
    if j eq 8 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad8.sav'
    if j eq 9 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad9.sav'
    if j eq 10 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad10.sav'
    if j eq 11 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad11.sav'
    if j eq 12 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad12.sav'
    if j eq 13 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad13.sav'
    if j eq 14 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad14.sav'
    if j eq 15 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad15.sav'
    if j eq 16 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad16.sav'
    if j eq 17 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad17.sav'
    if j eq 18 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad18.sav'
    if j eq 19 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad19.sav'
    if j eq 20 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad20.sav'
    if j eq 21 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad21.sav'
    if j eq 22 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad22.sav'
    if j eq 23 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad23.sav'
    if j eq 24 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad24.sav'
    if j eq 25 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad25.sav'
    if j eq 26 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad26.sav'
    if j eq 27 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad27.sav'
    if j eq 28 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad28.sav'

    ;8-dimensional_interpolation
    if i eq 0 then Table_Equivalent_width = Table_Equivalent_Pressure1
    if i eq 1 then Table_Equivalent_width = Table_Equivalent_Pressure2
    if i eq 2 then Table_Equivalent_width = Table_Equivalent_Pressure3
    if i eq 3 then Table_Equivalent_width = Table_Equivalent_Pressure4
    if i eq 4 then Table_Equivalent_width = Table_Equivalent_Pressure5
    if i eq 5 then Table_Equivalent_width = Table_Equivalent_Pressure6
    if i eq 6 then Table_Equivalent_width = Table_Equivalent_Pressure7
    if i eq 7 then Table_Equivalent_width = Table_Equivalent_Pressure8
    if i eq 8 then Table_Equivalent_width = Table_Equivalent_Pressure9
    if i eq 9 then Table_Equivalent_width = Table_Equivalent_Pressure10
    if i eq 10 then Table_Equivalent_width = Table_Equivalent_Pressure11
    if i eq 11 then Table_Equivalent_width = Table_Equivalent_Pressure12
    if i eq 12 then Table_Equivalent_width = Table_Equivalent_Pressure13
    if i eq 13 then Table_Equivalent_width = Table_Equivalent_Pressure14
    if i eq 14 then Table_Equivalent_width = Table_Equivalent_Pressure15
  
    y11111111=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7,j8)
    y21111111=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7,j8)
    y12111111=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7,j8)
    y11211111=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7,j8)
    y11121111=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7,j8)
    y11112111=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7,j8)
    y11111211=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7,j8)
    y22111111=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7,j8)
    y21211111=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7,j8)
    y21121111=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7,j8)
    y21112111=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7,j8)
    y21111211=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7,j8)
    y12211111=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7,j8)
    y12121111=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7,j8)
    y12112111=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7,j8)
    y12111211=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7,j8)
    y11221111=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7,j8)
    y11212111=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7,j8)
    y11211211=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7,j8)
    y11122111=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7,j8)
    y11121211=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7,j8)
    y11112211=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7,j8)
    y22211111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7,j8)
    y22121111=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7,j8)
    y22112111=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7,j8)
    y22111211=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7,j8)
    y21221111=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7,j8)
    y21212111=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7,j8)
    y21211211=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7,j8)
    y21122111=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7,j8)
    y21121211=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7,j8)
    y21112211=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7,j8)
    y12221111=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7,j8)
    y12212111=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7,j8)
    y12211211=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7,j8)
    y12122111=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7,j8)
    y12121211=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7,j8)
    y12112211=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7,j8)
    y11222111=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7,j8)
    y11221211=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7,j8)
    y11212211=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7,j8)
    y11122211=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7,j8)
    y11222211=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7,j8)
    y12122211=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7,j8)
    y12212211=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7,j8)
    y12221211=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7,j8)
    y12222111=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7,j8)
    y21122211=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7,j8)
    y21212211=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7,j8)
    y21221211=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7,j8)
    y21222111=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7,j8)
    y22112211=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7,j8)
    y22121211=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7,j8)
    y22122111=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7,j8)
    y22211211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7,j8)
    y22212111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7,j8)
    y22221111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7,j8)
    y12222211=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,j8)
    y21222211=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7,j8)
    y22122211=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7,j8)
    y22212211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7,j8)
    y22221211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7,j8)
    y22222111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7,j8)
    y22222211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,j8)
    y11111121=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7+1,j8)
    y21111121=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7+1,j8)
    y12111121=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7+1,j8)
    y11211121=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7+1,j8)
    y11121121=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7+1,j8)
    y11112121=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7+1,j8)
    y11111221=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7+1,j8)
    y22111121=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7+1,j8)
    y21211121=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7+1,j8)
    y21121121=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7+1,j8)
    y21112121=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7+1,j8)
    y21111221=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7+1,j8)
    y12211121=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7+1,j8)
    y12121121=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7+1,j8)
    y12112121=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7+1,j8)
    y12111221=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7+1,j8)
    y11221121=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7+1,j8)
    y11212121=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7+1,j8)
    y11211221=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7+1,j8)
    y11122121=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7+1,j8)
    y11121221=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7+1,j8)
    y11112221=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7+1,j8)
    y22211121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7+1,j8)
    y22121121=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7+1,j8)
    y22112121=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7+1,j8)
    y22111221=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7+1,j8)
    y21221121=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7+1,j8)
    y21212121=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7+1,j8)
    y21211221=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7+1,j8)
    y21122121=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7+1,j8)
    y21121221=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7+1,j8)
    y21112221=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7+1,j8)
    y12221121=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7+1,j8)
    y12212121=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7+1,j8)
    y12211221=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7+1,j8)
    y12122121=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7+1,j8)
    y12121221=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7+1,j8)
    y12112221=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7+1,j8)
    y11222121=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7+1,j8)
    y11221221=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7+1,j8)
    y11212221=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7+1,j8)
    y11122221=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7+1,j8)
    y11222221=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,j8)
    y12122221=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,j8)
    y12212221=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,j8)
    y12221221=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,j8)
    y12222121=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,j8)
    y21122221=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7+1,j8)
    y21212221=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7+1,j8)
    y21221221=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7+1,j8)
    y21222121=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7+1,j8)
    y22112221=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7+1,j8)
    y22121221=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7+1,j8)
    y22122121=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7+1,j8)
    y22211221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7+1,j8)
    y22212121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7+1,j8)
    y22221121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7+1,j8)
    y12222221=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,j8)
    y21222221=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,j8)
    y22122221=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,j8)
    y22212221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,j8)
    y22221221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,j8)
    y22222121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,j8)
    y22222221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,j8)
  
    y11111112=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7,j8+1)
    y21111112=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7,j8+1)
    y12111112=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7,j8+1)
    y11211112=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7,j8+1)
    y11121112=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7,j8+1)
    y11112112=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7,j8+1)
    y11111212=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7,j8+1)
    y22111112=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7,j8+1)
    y21211112=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7,j8+1)
    y21121112=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7,j8+1)
    y21112112=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7,j8+1)
    y21111212=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7,j8+1)
    y12211112=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7,j8+1)
    y12121112=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7,j8+1)
    y12112112=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7,j8+1)
    y12111212=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7,j8+1)
    y11221112=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7,j8+1)
    y11212112=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7,j8+1)
    y11211212=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7,j8+1)
    y11122112=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7,j8+1)
    y11121212=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7,j8+1)
    y11112212=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7,j8+1)
    y22211112=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7,j8+1)
    y22121112=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7,j8+1)
    y22112112=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7,j8+1)
    y22111212=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7,j8+1)
    y21221112=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7,j8+1)
    y21212112=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7,j8+1)
    y21211212=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7,j8+1)
    y21122112=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7,j8+1)
    y21121212=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7,j8+1)
    y21112212=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7,j8+1)
    y12221112=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7,j8+1)
    y12212112=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7,j8+1)
    y12211212=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7,j8+1)
    y12122112=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7,j8+1)
    y12121212=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7,j8+1)
    y12112212=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7,j8+1)
    y11222112=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7,j8+1)
    y11221212=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7,j8+1)
    y11212212=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7,j8+1)
    y11122212=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7,j8+1)
    y11222212=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7,j8+1)
    y12122212=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7,j8+1)
    y12212212=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7,j8+1)
    y12221212=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7,j8+1)
    y12222112=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7,j8+1)
    y21122212=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7,j8+1)
    y21212212=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7,j8+1)
    y21221212=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7,j8+1)
    y21222112=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7,j8+1)
    y22112212=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7,j8+1)
    y22121212=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7,j8+1)
    y22122112=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7,j8+1)
    y22211212=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7,j8+1)
    y22212112=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7,j8+1)
    y22221112=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7,j8+1)
    y12222212=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,j8+1)
    y21222212=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7,j8+1)
    y22122212=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7,j8+1)
    y22212212=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7,j8+1)
    y22221212=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7,j8+1)
    y22222112=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7,j8+1)
    y22222212=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,j8+1)
    y11111122=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7+1,j8+1)
    y21111122=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7+1,j8+1)
    y12111122=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7+1,j8+1)
    y11211122=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7+1,j8+1)
    y11121122=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7+1,j8+1)
    y11112122=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7+1,j8+1)
    y11111222=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7+1,j8+1)
    y22111122=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7+1,j8+1)
    y21211122=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7+1,j8+1)
    y21121122=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7+1,j8+1)
    y21112122=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7+1,j8+1)
    y21111222=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7+1,j8+1)
    y12211122=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7+1,j8+1)
    y12121122=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7+1,j8+1)
    y12112122=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7+1,j8+1)
    y12111222=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7+1,j8+1)
    y11221122=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7+1,j8+1)
    y11212122=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7+1,j8+1)
    y11211222=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7+1,j8+1)
    y11122122=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7+1,j8+1)
    y11121222=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7+1,j8+1)
    y11112222=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7+1,j8+1)
    y22211122=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7+1,j8+1)
    y22121122=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7+1,j8+1)
    y22112122=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7+1,j8+1)
    y22111222=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7+1,j8+1)
    y21221122=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7+1,j8+1)
    y21212122=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7+1,j8+1)
    y21211222=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7+1,j8+1)
    y21122122=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7+1,j8+1)
    y21121222=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7+1,j8+1)
    y21112222=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7+1,j8+1)
    y12221122=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7+1,j8+1)
    y12212122=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7+1,j8+1)
    y12211222=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7+1,j8+1)
    y12122122=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7+1,j8+1)
    y12121222=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7+1,j8+1)
    y12112222=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7+1,j8+1)
    y11222122=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7+1,j8+1)
    y11221222=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7+1,j8+1)
    y11212222=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7+1,j8+1)
    y11122222=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7+1,j8+1)
    y11222222=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,j8+1)
    y12122222=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,j8+1)
    y12212222=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,j8+1)
    y12221222=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,j8+1)
    y12222122=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,j8+1)
    y21122222=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7+1,j8+1)
    y21212222=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7+1,j8+1)
    y21221222=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7+1,j8+1)
    y21222122=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7+1,j8+1)
    y22112222=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7+1,j8+1)
    y22121222=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7+1,j8+1)
    y22122122=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7+1,j8+1)
    y22211222=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7+1,j8+1)
    y22212122=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7+1,j8+1)
    y22221122=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7+1,j8+1)
    y12222222=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,j8+1)
    y21222222=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,j8+1)
    y22122222=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,j8+1)
    y22212222=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,j8+1)
    y22221222=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,j8+1)
    y22222122=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,j8+1)
    y22222222=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,j8+1)
  
    t1 = (x1-x1a(j1)) / (x1a(j1+1)-x1a(j1))
    t2 = (x2-x2a(j2)) / (x2a(j2+1)-x2a(j2))
    t3 = (x3-x3a(j3)) / (x3a(j3+1)-x3a(j3))
    t4 = (x4-x4a(j4)) / (x4a(j4+1)-x4a(j4))
    t5 = (x5-x5a(j5)) / (x5a(j5+1)-x5a(j5))
    t6 = (x6-x6a(j6)) / (x6a(j6+1)-x6a(j6))
    t7 = (x7-x7a(j7)) / (x7a(j7+1)-x7a(j7))
    t8 = (x8-x8a(j8)) / (x8a(j8+1)-x8a(j8))
  
    y_tmp(i,j) = y11111111*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y21111111*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y12111111*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y11211111*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y11121111*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y11112111*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y11111211*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y22111111*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y21211111*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y21121111*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y21112111*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y21111211*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y12211111*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y12121111*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y12112111*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y12111211*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y11221111*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y11212111*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y11211211*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y11122111*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y11121211*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y11112211*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y22211111*t1*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y22121111*t1*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y22112111*t1*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y22111211*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y21221111*t1*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y21212111*t1*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y21211211*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y21122111*t1*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y21121211*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y21112211*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y12221111*(1.d - t1)*t2*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y12212111*(1.d - t1)*t2*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y12211211*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y12122111*(1.d - t1)*t2*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y12121211*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y12112211*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y11222111*(1.d - t1)*(1.d - t2)*t3*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y11221211*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y11212211*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y11122211*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y11222211*(1.d - t1)*(1.d - t2)*t3*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y12122211*(1.d - t1)*t2*(1.d - t3)*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y12212211*(1.d - t1)*t2*t3*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y12221211*(1.d - t1)*t2*t3*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y12222111*(1.d - t1)*t2*t3*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y21122211*t1*(1.d - t2)*(1.d - t3)*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y21212211*t1*(1.d - t2)*t3*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y21221211*t1*(1.d - t2)*t3*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y21222111*t1*(1.d - t2)*t3*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y22112211*t1*t2*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y22121211*t1*t2*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y22122111*t1*t2*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y22211211*t1*t2*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y22212111*t1*t2*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y22221111*t1*t2*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y12222211*(1.d - t1)*t2*t3*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y21222211*t1*(1.d - t2)*t3*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y22122211*t1*t2*(1.d - t3)*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y22212211*t1*t2*t3*(1.d - t4)*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y22221211*t1*t2*t3*t4*(1.d - t5)*t6*(1.d - t7)*(1.d - t8) +   $
      y22222111*t1*t2*t3*t4*t5*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
      y22222211*t1*t2*t3*t4*t5*t6*(1.d - t7)*(1.d - t8) +   $
      y11111121*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y21111121*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y12111121*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y11211121*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y11121121*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y11112121*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y11111221*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y22111121*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y21211121*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y21121121*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y21112121*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y21111221*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y12211121*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y12121121*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y12112121*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y12111221*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y11221121*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y11212121*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y11211221*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y11122121*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y11121221*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y11112221*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
      y22211121*t1*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y22121121*t1*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y22112121*t1*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y22111221*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y21221121*t1*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y21212121*t1*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y21211221*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y21122121*t1*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y21121221*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y21112221*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
      y12221121*(1.d - t1)*t2*t3*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y12212121*(1.d - t1)*t2*t3*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y12211221*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y12122121*(1.d - t1)*t2*(1.d - t3)*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y12121221*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y12112221*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
      y11222121*(1.d - t1)*(1.d - t2)*t3*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y11221221*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y11212221*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
      y11122221*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*t6*t7*(1.d - t8) +   $
      y11222221*(1.d - t1)*(1.d - t2)*t3*t4*t5*t6*t7*(1.d - t8) +   $
      y12122221*(1.d - t1)*t2*(1.d - t3)*t4*t5*t6*t7*(1.d - t8) +   $
      y12212221*(1.d - t1)*t2*t3*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
      y12221221*(1.d - t1)*t2*t3*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y12222121*(1.d - t1)*t2*t3*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y21122221*t1*(1.d - t2)*(1.d - t3)*t4*t5*t6*t7*(1.d - t8) +   $
      y21212221*t1*(1.d - t2)*t3*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
      y21221221*t1*(1.d - t2)*t3*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y21222121*t1*(1.d - t2)*t3*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y22112221*t1*t2*(1.d - t3)*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
      y22121221*t1*t2*(1.d - t3)*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y22122121*t1*t2*(1.d - t3)*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y22211221*t1*t2*t3*(1.d - t4)*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y22212121*t1*t2*t3*(1.d - t4)*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y22221121*t1*t2*t3*t4*(1.d - t5)*(1.d - t6)*t7*(1.d - t8) +   $
      y12222221*(1.d - t1)*t2*t3*t4*t5*t6*t7*(1.d - t8) +   $
      y21222221*t1*(1.d - t2)*t3*t4*t5*t6*t7*(1.d - t8) +   $
      y22122221*t1*t2*(1.d - t3)*t4*t5*t6*t7*(1.d - t8) +   $
      y22212221*t1*t2*t3*(1.d - t4)*t5*t6*t7*(1.d - t8) +   $
      y22221221*t1*t2*t3*t4*(1.d - t5)*t6*t7*(1.d - t8) +   $
      y22222121*t1*t2*t3*t4*t5*(1.d - t6)*t7*(1.d - t8) +   $
      y22222221*t1*t2*t3*t4*t5*t6*t7*(1.d - t8) +   $
      y11111112*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y21111112*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y12111112*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y11211112*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y11121112*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y11112112*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y11111212*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y22111112*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y21211112*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y21121112*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y21112112*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y21111212*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y12211112*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y12121112*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y12112112*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y12111212*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y11221112*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y11212112*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y11211212*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y11122112*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y11121212*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y11112212*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
      y22211112*t1*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y22121112*t1*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y22112112*t1*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y22111212*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y21221112*t1*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y21212112*t1*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y21211212*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y21122112*t1*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y21121212*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y21112212*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
      y12221112*(1.d - t1)*t2*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y12212112*(1.d - t1)*t2*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y12211212*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y12122112*(1.d - t1)*t2*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y12121212*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y12112212*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
      y11222112*(1.d - t1)*(1.d - t2)*t3*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y11221212*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y11212212*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
      y11122212*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*t6*(1.d - t7)*t8 +   $
      y11222212*(1.d - t1)*(1.d - t2)*t3*t4*t5*t6*(1.d - t7)*t8 +   $
      y12122212*(1.d - t1)*t2*(1.d - t3)*t4*t5*t6*(1.d - t7)*t8 +   $
      y12212212*(1.d - t1)*t2*t3*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
      y12221212*(1.d - t1)*t2*t3*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y12222112*(1.d - t1)*t2*t3*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y21122212*t1*(1.d - t2)*(1.d - t3)*t4*t5*t6*(1.d - t7)*t8 +   $
      y21212212*t1*(1.d - t2)*t3*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
      y21221212*t1*(1.d - t2)*t3*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y21222112*t1*(1.d - t2)*t3*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y22112212*t1*t2*(1.d - t3)*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
      y22121212*t1*t2*(1.d - t3)*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y22122112*t1*t2*(1.d - t3)*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y22211212*t1*t2*t3*(1.d - t4)*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y22212112*t1*t2*t3*(1.d - t4)*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y22221112*t1*t2*t3*t4*(1.d - t5)*(1.d - t6)*(1.d - t7)*t8 +   $
      y12222212*(1.d - t1)*t2*t3*t4*t5*t6*(1.d - t7)*t8 +   $
      y21222212*t1*(1.d - t2)*t3*t4*t5*t6*(1.d - t7)*t8 +   $
      y22122212*t1*t2*(1.d - t3)*t4*t5*t6*(1.d - t7)*t8 +   $
      y22212212*t1*t2*t3*(1.d - t4)*t5*t6*(1.d - t7)*t8 +   $
      y22221212*t1*t2*t3*t4*(1.d - t5)*t6*(1.d - t7)*t8 +   $
      y22222112*t1*t2*t3*t4*t5*(1.d - t6)*(1.d - t7)*t8 +   $
      y22222212*t1*t2*t3*t4*t5*t6*(1.d - t7)*t8 +   $
      y11111122*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y21111122*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y12111122*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y11211122*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y11121122*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y11112122*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
      y11111222*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
      y22111122*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y21211122*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y21121122*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y21112122*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
      y21111222*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
      y12211122*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y12121122*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y12112122*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
      y12111222*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
      y11221122*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y11212122*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
      y11211222*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
      y11122122*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*t7*t8 +   $
      y11121222*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*t7*t8 +   $
      y11112222*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*t7*t8 +   $
      y22211122*t1*t2*t3*(1.d - t4)*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y22121122*t1*t2*(1.d - t3)*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y22112122*t1*t2*(1.d - t3)*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
      y22111222*t1*t2*(1.d - t3)*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
      y21221122*t1*(1.d - t2)*t3*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y21212122*t1*(1.d - t2)*t3*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
      y21211222*t1*(1.d - t2)*t3*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
      y21122122*t1*(1.d - t2)*(1.d - t3)*t4*t5*(1.d - t6)*t7*t8 +   $
      y21121222*t1*(1.d - t2)*(1.d - t3)*t4*(1.d - t5)*t6*t7*t8 +   $
      y21112222*t1*(1.d - t2)*(1.d - t3)*(1.d - t4)*t5*t6*t7*t8 +   $
      y12221122*(1.d - t1)*t2*t3*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y12212122*(1.d - t1)*t2*t3*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
      y12211222*(1.d - t1)*t2*t3*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
      y12122122*(1.d - t1)*t2*(1.d - t3)*t4*t5*(1.d - t6)*t7*t8 +   $
      y12121222*(1.d - t1)*t2*(1.d - t3)*t4*(1.d - t5)*t6*t7*t8 +   $
      y12112222*(1.d - t1)*t2*(1.d - t3)*(1.d - t4)*t5*t6*t7*t8 +   $
      y11222122*(1.d - t1)*(1.d - t2)*t3*t4*t5*(1.d - t6)*t7*t8 +   $
      y11221222*(1.d - t1)*(1.d - t2)*t3*t4*(1.d - t5)*t6*t7*t8 +   $
      y11212222*(1.d - t1)*(1.d - t2)*t3*(1.d - t4)*t5*t6*t7*t8 +   $
      y11122222*(1.d - t1)*(1.d - t2)*(1.d - t3)*t4*t5*t6*t7*t8 +   $
      y11222222*(1.d - t1)*(1.d - t2)*t3*t4*t5*t6*t7*t8 +   $
      y12122222*(1.d - t1)*t2*(1.d - t3)*t4*t5*t6*t7*t8 +   $
      y12212222*(1.d - t1)*t2*t3*(1.d - t4)*t5*t6*t7*t8 +   $
      y12221222*(1.d - t1)*t2*t3*t4*(1.d - t5)*t6*t7*t8 +   $
      y12222122*(1.d - t1)*t2*t3*t4*t5*(1.d - t6)*t7*t8 +   $
      y21122222*t1*(1.d - t2)*(1.d - t3)*t4*t5*t6*t7*t8 +   $
      y21212222*t1*(1.d - t2)*t3*(1.d - t4)*t5*t6*t7*t8 +   $
      y21221222*t1*(1.d - t2)*t3*t4*(1.d - t5)*t6*t7*t8 +   $
      y21222122*t1*(1.d - t2)*t3*t4*t5*(1.d - t6)*t7*t8 +   $
      y22112222*t1*t2*(1.d - t3)*(1.d - t4)*t5*t6*t7*t8 +   $
      y22121222*t1*t2*(1.d - t3)*t4*(1.d - t5)*t6*t7*t8 +   $
      y22122122*t1*t2*(1.d - t3)*t4*t5*(1.d - t6)*t7*t8 +   $
      y22211222*t1*t2*t3*(1.d - t4)*(1.d - t5)*t6*t7*t8 +   $
      y22212122*t1*t2*t3*(1.d - t4)*t5*(1.d - t6)*t7*t8 +   $
      y22221122*t1*t2*t3*t4*(1.d - t5)*(1.d - t6)*t7*t8 +   $
      y12222222*(1.d - t1)*t2*t3*t4*t5*t6*t7*t8 +   $
      y21222222*t1*(1.d - t2)*t3*t4*t5*t6*t7*t8 +   $
      y22122222*t1*t2*(1.d - t3)*t4*t5*t6*t7*t8 +   $
      y22212222*t1*t2*t3*(1.d - t4)*t5*t6*t7*t8 +   $
      y22221222*t1*t2*t3*t4*(1.d - t5)*t6*t7*t8 +   $
      y22222122*t1*t2*t3*t4*t5*(1.d - t6)*t7*t8 +   $
      y22222222*t1*t2*t3*t4*t5*t6*t7*t8

  endfor
endfor

;Surface pressure grid
Pressure_grid = dblarr(15)
Pressure_grid(0) = alog(50d)
Pressure_grid(1) = alog(150d)
Pressure_grid(2) = alog(180d)
Pressure_grid(3) = alog(215d)
Pressure_grid(4) = alog(257d)
Pressure_grid(5) = alog(308d)
Pressure_grid(6) = alog(369d)
Pressure_grid(7) = alog(442d)
Pressure_grid(8) = alog(529d)
Pressure_grid(9) = alog(633d)
Pressure_grid(10) = alog(758d)
Pressure_grid(11) = alog(907d)
Pressure_grid(12) = alog(1096d)
Pressure_grid(13) = alog(1300d)
Pressure_grid(14) = alog(1500d)

for i = 0, nx-1 do Y(i) = interpol(Y_tmp(*,i), Pressure_grid, alog(P(0)))

restore, '/work1/LUT/Common/specmars_CO2_update.sav'
Y = Y/specmars

; 2 gaussian continuum
;Y = Y * pyroxenes(x, P(9:11))

; simple continuum
Y = Y * poly(findgen(nx), P(9:10))

; 線形continuum
;Y = Y * (P(9)*x + P(10))

return, Y
end

;---------------------------
function pyroxenes, x, p
;---------------------------
y1 = 1d - gauss1(x, [1.9d, 0.5d / (2 * sqrt(2*alog(2))), P(0)], /peak)
y2 = 1d - gauss1(x, [2.3d, 0.56d / (2 * sqrt(2*alog(2))), P(1)], /peak)
y = (y1 + y2)/2d * P(2)
return, y
end




Pro retrieval_pressure_fitting_v2


; define color platte
Set_Plot, 'x'
device, retain=1, decomposed=0
loadct, 39


; path 
path = '/data2/omega/sav/'
path2 = '/work1/LUT/SP/table/absorption/'
restore, path+'ORB0363_3.sav'
;restore, path+'ORB0030_1.sav'
;restore, path+'ORB1201_3.sav'
restore, path + 'specmars.sav'

;reference spectrum
ref = dblarr(2,3539)
openr, lun, '/work1/LUT/Common/psg_trn.txt', /get_lun
for i = 0, 3539-1 do begin
  readf, lun, a, b
  ref(0,i) = a
  ref(1,i) = b
endfor
free_lun,lun

; ind=76094,xind=62, yind=594, lati:22.705700 N , longi:311.76300 E (48.237 W)  [Forget+, retrievalすると852 Pa] ORB0363_3
ind=where_xyz(longi ge 311.73 and longi le 311.78 and lati ge 22.7 and lati le 22.72, xind=xind, yind=yind)

; ind=32015, xind=15, yind=1000, lati: 48.431797 S, longi: 60.808998 E [Forget+, retrievalすると1036 Pa] ORB0030_1
;ind=where_xyz(longi ge 60.79 and longi le 60.81 and lati ge -48.44 and lati le -48.43,xind=xind,yind=yind)

; ind=37135, xind=15, yind=1160, lati:7.764°S, longi:24.980°E　[Forget+, retrievalすると470 Pa] ORB1201_3
;ind=where_xyz(longi ge 24.979 and longi le 24.981 and lati ge -7.765 and lati le -7.763,xind=xind,yind=yind)

;test getting surface feature 
x0 = reform(wvl(0:127))
y0 = reform(jdat(xind,0:127,yind)/specmars(0:127))
nanserch=where(y0 ge 0 and y0 le 0.0001)
y0(nanserch)=!VALUES.F_NAN

ref_Mars = interpol(ref(1,*), ref(0,*), x0, /nan)
;good = where(x0 ge 1.2 and x0 le 2.6 and ref_Mars ge 0.995)
good = where(x0 ge 1.2 and x0 le 2.6 and ref_Mars ge 0.99)

pi = replicate({step:0d, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 3)
start = [0.06, 0.06, median(y0(good))]
Result_Fit0 = MPFITFUN('pyroxenes', x0(good), y0(good), y0(good)*1d-2, start, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit, /nan)
F0 = yfit

; CO2 absorption emission line
CO2=where(wvl ge 1.8 and wvl le 2.2)
wvl=wvl[CO2]
jdat=jdat(*,CO2,*)
specmars = specmars(CO2)


;    ; ----- MCD ------
lat = median(lati(xind,yind),/double)
lon = median(longi(xind,yind),/double)
Ls = SOLAR_LONGITUDE
Loct = 12. + (lon - SUB_SOLAR_LONGITUDE)/15. ;! TBC!!! !
dust = 2  ; our best guess MY24 scenario, with solar average conditions
hrkey = 1 ; set high resolution mode on (hrkey=0 to set high resolution off)
zkey = 3    ; specify that xz is the altitude above surface (m)
xz = 0. ; (array of altitude: step=2000 m start: 5m)
datekey = 1       ; <integer> type of input date (1=Mars date)
xdate = ls        ; <double precision> date (IF datekey = 1 : Value of Ls [deg.])
dset ='/work1/LUT/MCD/MCD5.3/data/' ;‘MCD_DATA/’  ; <character*50> data set
scena = 1         ; <integer> scenario (1 = Climatology ave solar)
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
cd, '/work1/LUT/MCD/MCD5.3/mcd/idl/'
a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
  dset,scena,perturkey,seedin,gwlength,extvarkeys, $
  pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
scaleH = extvar(13)
Z1 = scaleH * 0.1d
Z2 = scaleH * 4d
dust_opacity = extvar(36)*1.2090217E+00/4.6232791E+00
ice_opacity = 0.d;*4.2410289E-01/***  ;!TBD!
SP = PRES
xz=[Z1] ; (array of altitude with step=2000 m)
a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
  dset,scena,perturkey,seedin,gwlength,extvarkeys, $
  pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
TA = temp
xz=[Z2] ; (array of altitude with step=2000 m)
a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
  dset,scena,perturkey,seedin,gwlength,extvarkeys, $
  pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
TB = temp
;;-----------------

SZA = reform(geocube(xind,8,yind))*1.e-4
EA = reform(geocube(xind,9,yind))*1.e-4
PA = reform(geocube(xind,10,yind))*1.e-4
Albedo_input = jdat(xind,0,yind)/specmars(0) / cos(geocube(xind,8,yind)*1e-4*!DTOR)


jdat(*,0:3,*)= !VALUES.F_NAN
jdat(*,8,*)= !VALUES.F_NAN
jdat(*,17,*)= !VALUES.F_NAN
jdat(*,24,*)= !VALUES.F_NAN
jdat(*,27,*)= !VALUES.F_NAN

x = wvl
y = reform(jdat(xind(0), *, yind(0)))
y = y/specmars

y(0:3) = -0d/0d
y(8) = -0d/0d
y(17) = -0d/0d
y(24) = -0d/0d
y(27) = -0d/0d


;retrieval
;pi = replicate({step:0d, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 12)
;;start = [SP, TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input, Result_Fit0(0), Result_Fit0(1), Result_Fit0(2)]
;start = [SP, TA, TB, SZA, EA, PA, 0d, 0d, 0.29, Result_Fit0(0), Result_Fit0(1), Result_Fit0(2)]
;pi(1:8).fixed = 1

pi = replicate({step:0d, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 11)
;start = [SP, TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input, Result_Fit0(0), Result_Fit0(1), Result_Fit0(2)]
;start = [SP, TA, TB, SZA, EA, PA, 0d, 0d, 0.29, 1d, 0d]
;pi(1:11).fixed = 1

; Kazama LUT
; A = 0.15
start = [SP, TA, TB, SZA, EA, PA, 0.04d, 0d, 0.15, 1d, 0d]
pi(1:8).fixed = 1
err = y*1d-3
err(*) = median(y)*1d-3
Result_Fit = MPFITFUN('forward', x, y, err, start, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit, /nan)
F = yfit
con = poly(findgen(29), Result_Fit(9:10))

start1 = [SP, TA, TB, SZA, EA, PA, 0.44d, 0d, 0.15, 1d, 0d]
pi(1:8).fixed = 1
Result_Fit1 = MPFITFUN('forward', x, y, err, start1, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit1, /nan)
F1 = yfit1
con1 = poly(findgen(29), Result_Fit1(9:10))

start2 = [SP, TA, TB, SZA, EA, PA, 0.24d, 0d, 0.15, 1d, 0d]
pi(1:8).fixed = 1
Result_Fit2 = MPFITFUN('forward', x, y, err, start2, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit2, /nan)
F2 = yfit2
con2 = poly(findgen(29), Result_Fit2(9:10))

; A = 0.29
start3 = [SP, TA, TB, SZA, EA, PA, 0.04d, 0d, 0.29, 1d, 0d]
pi(1:8).fixed = 1
Result_Fit3 = MPFITFUN('forward', x, y, err, start3, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit3, /nan)
F3 = yfit3
con3 = poly(findgen(29), Result_Fit3(9:10))

start4 = [SP, TA, TB, SZA, EA, PA, 0.44d, 0d, 0.29, 1d, 0d]
pi(1:8).fixed = 1
Result_Fit4 = MPFITFUN('forward', x, y, err, start4, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit4, /nan)
F4 = yfit4
con4 = poly(findgen(29), Result_Fit4(9:10))

start5 = [SP, TA, TB, SZA, EA, PA, 0.24d, 0d, 0.29, 1d, 0d]
pi(1:8).fixed = 1
Result_Fit5 = MPFITFUN('forward', x, y, err, start5, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit5, /nan)
F5 = yfit5
con5 = poly(findgen(29), Result_Fit5(9:10))

; Forgetの値
; A = 0.15
start6 = [867.5d, TA, TB, SZA, EA, PA, 0.44d, 0d, 0.15, 1d, 0d]
pi(0:8).fixed = 1
Result_Fit6 = MPFITFUN('forward', x, y, err, start6, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit6, /nan)
F6 = yfit6
con6 = poly(findgen(29), Result_Fit6(9:10))

start7 = [822.5d, TA, TB, SZA, EA, PA, 0.24d, 0d, 0.15, 1d, 0d]
pi(0:8).fixed = 1
Result_Fit7 = MPFITFUN('forward', x, y, err, start7, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit7, /nan)
F7 = yfit7
con7 = poly(findgen(29), Result_Fit7(9:10))

start8 = [786.5d, TA, TB, SZA, EA, PA, 0.04d, 0d, 0.15, 1d, 0d]
pi(0:8).fixed = 1
Result_Fit8 = MPFITFUN('forward', x, y, err, start8, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit8, /nan)
F8 = yfit8
con8 = poly(findgen(29), Result_Fit8(9:10))

; A = 0.29
start9 = [837.5d, TA, TB, SZA, EA, PA, 0.44d, 0d, 0.29, 1d, 0d]
pi(0:8).fixed = 1
Result_Fit9 = MPFITFUN('forward', x, y, err, start9, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit9, /nan)
F9 = yfit9
con9 = poly(findgen(29), Result_Fit9(9:10))

start10 = [822.5d, TA, TB, SZA, EA, PA, 0.24d, 0d, 0.29, 1d, 0d]
pi(0:8).fixed = 1
Result_Fit10 = MPFITFUN('forward', x, y, err, start10, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit10, /nan)
F10 = yfit10
con10 = poly(findgen(29), Result_Fit10(9:10))

start11 = [815.5d, TA, TB, SZA, EA, PA, 0.04d, 0d, 0.29, 1d, 0d]
pi(0:8).fixed = 1
Result_Fit11 = MPFITFUN('forward', x, y, err, start11, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit11, /nan)
F11 = yfit11
con11 = poly(findgen(29), Result_Fit11(9:10))

;start2 = [852d, TA, TB, SZA, EA, PA, 0d, 0d, 0.29, 1d, 0d]
;pi(0:8).fixed = 1
;err = y*1d-1
;Result_Fit = MPFITFUN('forward', x, y, err, start2, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, yfit=yfit, /nan)
;F2 = yfit


f(0:3) = -0d/0d
f(8) = -0d/0d
f(17) = -0d/0d
f(24) = -0d/0d
f(27) = -0d/0d

con(0:3) = -0d/0d
con(8) = -0d/0d
con(17) = -0d/0d
con(24) = -0d/0d
con(27) = -0d/0d

f1(0:3) = -0d/0d
f1(8) = -0d/0d
f1(17) = -0d/0d
f1(24) = -0d/0d
f1(27) = -0d/0d

con1(0:3) = -0d/0d
con1(8) = -0d/0d
con1(17) = -0d/0d
con1(24) = -0d/0d
con1(27) = -0d/0d

f2(0:3) = -0d/0d
f2(8) = -0d/0d
f2(17) = -0d/0d
f2(24) = -0d/0d
f2(27) = -0d/0d

con2(0:3) = -0d/0d
con2(8) = -0d/0d
con2(17) = -0d/0d
con2(24) = -0d/0d
con2(27) = -0d/0d

f3(0:3) = -0d/0d
f3(8) = -0d/0d
f3(17) = -0d/0d
f3(24) = -0d/0d
f3(27) = -0d/0d

con3(0:3) = -0d/0d
con3(8) = -0d/0d
con3(17) = -0d/0d
con3(24) = -0d/0d
con3(27) = -0d/0d

f4(0:3) = -0d/0d
f4(8) = -0d/0d
f4(17) = -0d/0d
f4(24) = -0d/0d
f4(27) = -0d/0d

con4(0:3) = -0d/0d
con4(8) = -0d/0d
con4(17) = -0d/0d
con4(24) = -0d/0d
con4(27) = -0d/0d

f5(0:3) = -0d/0d
f5(8) = -0d/0d
f5(17) = -0d/0d
f5(24) = -0d/0d
f5(27) = -0d/0d

con5(0:3) = -0d/0d
con5(8) = -0d/0d
con5(17) = -0d/0d
con5(24) = -0d/0d
con5(27) = -0d/0d

f6(0:3) = -0d/0d
f6(8) = -0d/0d
f6(17) = -0d/0d
f6(24) = -0d/0d
f6(27) = -0d/0d

con6(0:3) = -0d/0d
con6(8) = -0d/0d
con6(17) = -0d/0d
con6(24) = -0d/0d
con6(27) = -0d/0d

f7(0:3) = -0d/0d
f7(8) = -0d/0d
f7(17) = -0d/0d
f7(24) = -0d/0d
f7(27) = -0d/0d

con7(0:3) = -0d/0d
con7(8) = -0d/0d
con7(17) = -0d/0d
con7(24) = -0d/0d
con7(27) = -0d/0d

F8(0:3) = -0d/0d
F8(8) = -0d/0d
F8(17) = -0d/0d
F8(24) = -0d/0d
F8(27) = -0d/0d

con8(0:3) = -0d/0d
con8(8) = -0d/0d
con8(17) = -0d/0d
con8(24) = -0d/0d
con8(27) = -0d/0d

F9(0:3) = -0d/0d
F9(8) = -0d/0d
F9(17) = -0d/0d
F9(24) = -0d/0d
F9(27) = -0d/0d

con9(0:3) = -0d/0d
con9(8) = -0d/0d
con9(17) = -0d/0d
con9(24) = -0d/0d
con9(27) = -0d/0d

F10(0:3) = -0d/0d
F10(8) = -0d/0d
F10(17) = -0d/0d
F10(24) = -0d/0d
F10(27) = -0d/0d

con10(0:3) = -0d/0d
con10(8) = -0d/0d
con10(17) = -0d/0d
con10(24) = -0d/0d
con10(27) = -0d/0d

F11(0:3) = -0d/0d
F11(8) = -0d/0d
F11(17) = -0d/0d
F11(24) = -0d/0d
F11(27) = -0d/0d

con11(0:3) = -0d/0d
con11(8) = -0d/0d
con11(17) = -0d/0d
con11(24) = -0d/0d
con11(27) = -0d/0d

;f2(0:3) = -0d/0d
;f2(8) = -0d/0d
;f2(17) = -0d/0d
;f2(24) = -0d/0d
;f2(27) = -0d/0d

good = where(FINITE(y) eq 1)
window,6, xs=800,ys=800
; plot, x(good), y(good), yr=[-0.1,0.3], back=255, color=0, thick=3, psym=-1, xr=[1.85, 2.2], xs=1
;plot, x(good), y(good), yr=[-0.1,0.4],back=255, color=0, thick=3, psym=-1, xr=[1.85, 2.2], xs=1
;oplot, x(good), f2(good), color=60, thick=3, psym=-1
plot, x(good), f(good),  back=255, color=0, thick=3, psym=-1, xr=[1.85, 2.2], xs=1, yr = [-0.1, 2.5], ys = 1
oplot, x(good), f1(good), color=60, thick=2, psym=-1, linestyle=2
oplot, x(good), f2(good), color=254, thick=2, psym=-1, linestyle=2
oplot, x(good), con(good), color=0, thick=2, psym=-1, linestyle=2
oplot, x(good), con1(good), color=60, thick=2, psym=-1, linestyle=2
oplot, x(good), con2(good), color=254, thick=2, psym=-1, linestyle=2
;oplot, x(good), (y(good)-f(good))*5, color=200, thick=2, psym=-1
;oplot, x(good), (y(good)-f2(good))*5, color=180, thick=2, psym=-1
;oplot, x(good), pyroxenes(x(good), Result_Fit(9:11)), color=100
;xyouts, 2.05, 0.1, 'SP='+strcompress(Result_Fit(0)), charsize=2, color=0
stop


end
