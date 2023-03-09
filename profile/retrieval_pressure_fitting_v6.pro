;------------------------------------------------------------------------------------------------------------------------
;
; ある1地点でのLUT + fitting をするプログラム
;
; create by Akira Kazama
; 
; retrieval_pressure_fitting_v1　 ::2022.11.16 Wed 10:52:00
; ver1: LUT update前に対応するプログラム
;
; retrieval_pressure_fitting_v2　 ::2022.11.16 Wed 11:52:00
; ver2: LUT update後に対応するプログラム
;
; retrieval_pressure_fitting_v3　 ::2022.11.16 Wed 12:55:00
; ver3: 高速化導入、多次元補間を外に出す
;
; retrieval_pressure_fitting_v4　 ::2022.11.25 Fri 10:07:00
; ver4: Albedoも同時にfittingをする
;
; retrieval_pressure_fitting_v5　 ::2022.11.25 Fri 13:45:00
; ver5: fitting (気圧だけFree) によるpressure mapを導入
; !TBD! MCDの読み出し分のみだけを回すようにする 
;
; retrieval_pressure_fitting_v6　 ::2022.12.12 Mon 15:38:00
; ver6: fitting (気圧+AlbedoがFree) によるpressure mapを導入
; ver6.1: continuum + best fit spectrumもsav fileに追加 [Edit 2022.12.13]
; !TBD! MCDの読み出し分のみだけを回すようにする 
;
; +++
;------------------------------------------------------------------------------------------------------------------------

;-------------------------------------------------------------------------------
function multifunction, TA, TB, SZA, EA, PA, Dust, Waterice
;-------------------------------------------------------------------------------

T1 = TA
T2 = TB

;note: y_tmp =n_elements(pressure, wavelength, Albedo)
y_tmp = dblarr(15, 29, 7)


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

X1 = (T1)^1.5d ;T1 as test
X2 = (T2)^1.5d ;T2 as test
x3 = exp(-cos(SZA/180d*!dpi))
x4 = exp(-cos(EA/180d*!dpi))
x5 = -cos(PA/180d*!dpi) + 1.d
x6 = Dust ;reduce_dust(loop)*1.86
x7 = Waterice; reduce_waterice(loop)*3.87

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

j1 = j1 -1
j2 = j2 -1
j3 = j3 -1
j4 = j4 -1
j5 = j5 -1
j6 = j6 -1
j7 = j7 -1

if j1 lt 0 then stop
if j2 lt 0 then stop
if j3 lt 0 then stop
if j4 lt 0 then stop
if j5 lt 0 then stop
if j6 lt 0 then stop
if j7 lt 0 then stop
if j1+1 ge 5 then stop
if j2+1 ge 3 then stop
if j3+1 ge 6 then stop
if j4+1 ge 5 then stop
if j5+1 ge 5 then stop
if j6+1 ge 6 then stop
if j7+1 ge 3 then stop

;mutli-dimensional interpolations
for i = 0, 15-1 do begin
  for j = 0, 29-1 do begin
  
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

  y11111111=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7,*)
  y21111111=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7,*)
  y12111111=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7,*)
  y11211111=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7,*)
  y11121111=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7,*)
  y11112111=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7,*)
  y11111211=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7,*)
  y22111111=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7,*)
  y21211111=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7,*)
  y21121111=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7,*)
  y21112111=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7,*)
  y21111211=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7,*)
  y12211111=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7,*)
  y12121111=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7,*)
  y12112111=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7,*)
  y12111211=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7,*)
  y11221111=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7,*)
  y11212111=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7,*)
  y11211211=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7,*)
  y11122111=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7,*)
  y11121211=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7,*)
  y11112211=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7,*)
  y22211111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7,*)
  y22121111=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7,*)
  y22112111=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7,*)
  y22111211=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7,*)
  y21221111=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7,*)
  y21212111=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7,*)
  y21211211=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7,*)
  y21122111=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7,*)
  y21121211=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7,*)
  y21112211=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7,*)
  y12221111=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7,*)
  y12212111=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7,*)
  y12211211=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7,*)
  y12122111=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7,*)
  y12121211=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7,*)
  y12112211=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7,*)
  y11222111=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7,*)
  y11221211=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7,*)
  y11212211=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7,*)
  y11122211=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7,*)
  y11222211=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7,*)
  y12122211=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7,*)
  y12212211=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7,*)
  y12221211=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7,*)
  y12222111=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7,*)
  y21122211=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7,*)
  y21212211=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7,*)
  y21221211=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7,*)
  y21222111=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7,*)
  y22112211=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7,*)
  y22121211=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7,*)
  y22122111=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7,*)
  y22211211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7,*)
  y22212111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7,*)
  y22221111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7,*)
  y12222211=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,*)
  y21222211=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7,*)
  y22122211=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7,*)
  y22212211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7,*)
  y22221211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7,*)
  y22222111=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7,*)
  y22222211=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7,*)
  y11111121=Table_Equivalent_width(j1,j2,j3,j4,j5,j6,j7+1,*)
  y21111121=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6,j7+1,*)
  y12111121=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6,j7+1,*)
  y11211121=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6,j7+1,*)
  y11121121=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6,j7+1,*)
  y11112121=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6,j7+1,*)
  y11111221=Table_Equivalent_width(j1,j2,j3,j4,j5,j6+1,j7+1,*)
  y22111121=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6,j7+1,*)
  y21211121=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6,j7+1,*)
  y21121121=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6,j7+1,*)
  y21112121=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6,j7+1,*)
  y21111221=Table_Equivalent_width(j1+1,j2,j3,j4,j5,j6+1,j7+1,*)
  y12211121=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6,j7+1,*)
  y12121121=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6,j7+1,*)
  y12112121=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6,j7+1,*)
  y12111221=Table_Equivalent_width(j1,j2+1,j3,j4,j5,j6+1,j7+1,*)
  y11221121=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6,j7+1,*)
  y11212121=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6,j7+1,*)
  y11211221=Table_Equivalent_width(j1,j2,j3+1,j4,j5,j6+1,j7+1,*)
  y11122121=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6,j7+1,*)
  y11121221=Table_Equivalent_width(j1,j2,j3,j4+1,j5,j6+1,j7+1,*)
  y11112221=Table_Equivalent_width(j1,j2,j3,j4,j5+1,j6+1,j7+1,*)
  y22211121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6,j7+1,*)
  y22121121=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6,j7+1,*)
  y22112121=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6,j7+1,*)
  y22111221=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5,j6+1,j7+1,*)
  y21221121=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6,j7+1,*)
  y21212121=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6,j7+1,*)
  y21211221=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5,j6+1,j7+1,*)
  y21122121=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6,j7+1,*)
  y21121221=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5,j6+1,j7+1,*)
  y21112221=Table_Equivalent_width(j1+1,j2,j3,j4,j5+1,j6+1,j7+1,*)
  y12221121=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6,j7+1,*)
  y12212121=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6,j7+1,*)
  y12211221=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5,j6+1,j7+1,*)
  y12122121=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6,j7+1,*)
  y12121221=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5,j6+1,j7+1,*)
  y12112221=Table_Equivalent_width(j1,j2+1,j3,j4,j5+1,j6+1,j7+1,*)
  y11222121=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6,j7+1,*)
  y11221221=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5,j6+1,j7+1,*)
  y11212221=Table_Equivalent_width(j1,j2,j3+1,j4,j5+1,j6+1,j7+1,*)
  y11122221=Table_Equivalent_width(j1,j2,j3,j4+1,j5+1,j6+1,j7+1,*)
  y11222221=Table_Equivalent_width(j1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,*)
  y12122221=Table_Equivalent_width(j1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,*)
  y12212221=Table_Equivalent_width(j1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,*)
  y12221221=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,*)
  y12222121=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,*)
  y21122221=Table_Equivalent_width(j1+1,j2,j3,j4+1,j5+1,j6+1,j7+1,*)
  y21212221=Table_Equivalent_width(j1+1,j2,j3+1,j4,j5+1,j6+1,j7+1,*)
  y21221221=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5,j6+1,j7+1,*)
  y21222121=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6,j7+1,*)
  y22112221=Table_Equivalent_width(j1+1,j2+1,j3,j4,j5+1,j6+1,j7+1,*)
  y22121221=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5,j6+1,j7+1,*)
  y22122121=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6,j7+1,*)
  y22211221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5,j6+1,j7+1,*)
  y22212121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6,j7+1,*)
  y22221121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6,j7+1,*)
  y12222221=Table_Equivalent_width(j1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,*)
  y21222221=Table_Equivalent_width(j1+1,j2,j3+1,j4+1,j5+1,j6+1,j7+1,*)
  y22122221=Table_Equivalent_width(j1+1,j2+1,j3,j4+1,j5+1,j6+1,j7+1,*)
  y22212221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4,j5+1,j6+1,j7+1,*)
  y22221221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5,j6+1,j7+1,*)
  y22222121=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6,j7+1,*)
  y22222221=Table_Equivalent_width(j1+1,j2+1,j3+1,j4+1,j5+1,j6+1,j7+1,*)

  t1 = (x1-x1a(j1)) / (x1a(j1+1)-x1a(j1))
  t1 = replicate(t1, 7)
  t2 = (x2-x2a(j2)) / (x2a(j2+1)-x2a(j2))
  t2 = replicate(t2, 7)
  t3 = (x3-x3a(j3)) / (x3a(j3+1)-x3a(j3))
  t3 = replicate(t3, 7)
  t4 = (x4-x4a(j4)) / (x4a(j4+1)-x4a(j4))
  t4 = replicate(t4, 7)
  t5 = (x5-x5a(j5)) / (x5a(j5+1)-x5a(j5))
  t5 = replicate(t5, 7)
  t6 = (x6-x6a(j6)) / (x6a(j6+1)-x6a(j6))
  t6 = replicate(t6, 7)
  t7 = (x7-x7a(j7)) / (x7a(j7+1)-x7a(j7))
  t7 = replicate(t7, 7)
  t8 = 0d;(x8-x8a(j8)) / (x8a(j8+1)-x8a(j8))
  t8 = replicate(t8, 7)

  y_tmp(i,j,*) = y11111111*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
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
    y22222121*t1*t2*t3*t4*t5*(1.d - t6)*t7*(1.d - t8)

  endfor
endfor

return, y_tmp
end


;---------------------------
function forward, x, p
;---------------------------
;x: wavenumber
;y: pressureによって変化する15個の放射輝度
;return: radiance
;---------------------------

restore, '/work1/LUT/Common/fitting_ytmp_albedo_1.sav'
nx = n_elements(x)
y = dblarr(nx)

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

Albedo_grid = dblarr(7) ; Surface Albedo
Albedo_grid(0) = 0.05d
Albedo_grid(1) = 0.1d
Albedo_grid(2) = 0.2d
Albedo_grid(3) = 0.3d
Albedo_grid(4) = 0.4d
Albedo_grid(5) = 0.5d
Albedo_grid(6) = 0.6d

for I = 0, 15-1 do begin
  F = Pressure_grid(I) - alog(P(0))
  if (F gt 0.0d0) then j1 = I
  if (F gt 0.0d0) then goto, skip01
endfor
skip01:

for I = 0, 7-1 do begin
  F = Albedo_grid(I) - P(1)
  if (F gt 0.0d0) then j2 = I
  if (F gt 0.0d0) then goto, skip02
endfor
skip02:

j1 = j1 -1
j2 = j2 -1

;if j1 lt 0 then stop
;if j2 lt 0 then stop

t1 = (alog(P(0))-Pressure_grid(j1)) / (Pressure_grid(j1+1)-Pressure_grid(j1))
t2 = (P(1)-Albedo_grid(j2)) / (Albedo_grid(j2+1)-Albedo_grid(j2))
Y(*) = Y_tmp(j1,*,j2)*(1.d - t1)*(1.d - t2) + Y_tmp(j1+1,*,j2)*t1*(1.d - t2) + Y_tmp(j1,*,j2+1)*(1.d - t1)*t2 + Y_tmp(j1+1,*,j2+1)*t1*t2

restore, '/work1/LUT/Common/specmars_CO2_update.sav'
;Y = Y/specmars

; 2 gaussian continuum
Y = Y * ( pyroxenes(x, P(2:4)) / mean(pyroxenes(x, P(2:4))) )

; simple continuum
;Y = Y * poly(findgen(nx)-float(nx/2), [P(2:3)])

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





Pro retrieval_pressure_fitting_v6

start_time = systime(1)
; define color platte
Set_Plot, 'x'
device, retain=1, decomposed=0
loadct, 39

; path 
path = '/data2/omega/sav/'
path2 = '/work1/LUT/SP/table/absorption/'

; 相対誤差の比較のためのORB [Forget+ 2007] Fig 13
restore, path+'ORB0920_3.sav'
;restore, path+'ORB0931_3.sav'

restore, path + 'specmars.sav'


;reference spectrum
; W guassiunのときに使う
; =========
ref = dblarr(2,3539)
openr, lun, '/work1/LUT/Common/psg_trn.txt', /get_lun
for i = 0, 3539-1 do begin
  readf, lun, a, b
  ref(0,i) = a
  ref(1,i) = b
endfor
free_lun,lun
; =========

;ind = where_xyz(longi ge 273 and longi le 277 and lati ge 53 and lati le 56, xind=xind, yind=yind)
; Block 1
ind = where_xyz(longi ge 274 and longi le 277 and lati ge 54 and lati le 56, xind=xind, yind=yind)

; 全部のやつ
ip = n_elements(LATI(*,0))
io = n_elements(LATI(0,*))

; loop変数
ip_b = min(xind)
ip_a = max(xind)
io_b = min(yind)
io_a = max(yind)

; create array
pressure = dblarr(ip,io)
MCDpressure = dblarr(ip,io)
Albedomap = dblarr(ip,io)
dustmap = dblarr(ip,io)
TAmap = dblarr(ip,io)
altitude = dblarr(ip,io)
bestfit = dblarr(ip,io,29)
continuum = dblarr(ip,io,29)

y_tmp = dblarr(15, 29)

for l = ip_b, ip_a -1 do begin ;loop for slit scan
  for k = io_b, io_a -1 do begin ;test getting surface feature

    start_time = systime(1)
    ;   ----- MCD ------
    lat = lati(l,k)
    lon = longi(l,k)
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

    ; W gaussiunで使う
    ; =========
    ; continuum の決定
    x0 = reform(wvl(0:127))
    y0 = reform(jdat(l,0:127,k)/specmars(0:127))
    nanserch=where(y0 ge 0 and y0 le 0.0001)
    y0(nanserch)=!VALUES.F_NAN

    ref_Mars = interpol(ref(1,*), ref(0,*), x0, /nan)
    good = where(x0 ge 1.2 and x0 le 2.6 and ref_Mars ge 0.9996)

    ; continuumの定数を決める
    pi = replicate({step:0d, fixed:0, limited:[0,0], limits:[0.D,0.D], MPPRINT:[0]}, 3)
    start = [0.06d, 0.06d, median(y0(good))]
    Result_Fit0 = MPFITFUN('pyroxenes', x0(good), y0(good), y0(good)*1d-2, start, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, NPRINT=0, status=status, yfit=yfit, /nan)
    pi(*).MPPRINT = 0
    F0 = yfit
    ; =========

    ; CO2 absorption emission line
    ; observation spectrum
    CO2=where(wvl ge 1.8 and wvl le 2.2)
    wvl_CO2=wvl[CO2]
    jdat_CO2=jdat(*,CO2,*)
    specmars_CO2 = specmars(CO2)

    ; Calc spectrum 計算
    SZA = reform(geocube(l,8,k))*1.e-4
    EA = reform(geocube(l,9,k))*1.e-4
    PA = reform(geocube(l,10,k))*1.e-4
    altitude(l,k) = reform(geocube(l,12,k))

    Albedo_input = jdat_CO2(l,0,k)/specmars_CO2(0) / cos(geocube(l,8,k)*1e-4*!DTOR)
    
    Waterice = ice_opacity
    Albedo = Albedo_input
    dust = dust_opacity

    MCDpressure(l,k) = SP
    dustmap(l,k) = dust 
    TAmap(l,k) = TA
    ;dust = 0.0744d ;test
    ;dust = 0d  ;CS

    ; not use spectrum
    jdat_CO2(*,0:3,*)= !VALUES.F_NAN
    jdat_CO2(*,8,*)= !VALUES.F_NAN
    jdat_CO2(*,17,*)= !VALUES.F_NAN
    jdat_CO2(*,24,*)= !VALUES.F_NAN
    jdat_CO2(*,27,*)= !VALUES.F_NAN

    x = wvl_CO2
    y_obs = reform(jdat_CO2(l, *, k))
    y_obs = y_obs/specmars_CO2

    y_tmp = multifunction(TA, TB, SZA, EA, PA, dust, ice_opacity)
    if y_tmp(0,0) eq 0.0d then goto, skip10

    save, y_tmp, filename='/work1/LUT/Common/fitting_ytmp_albedo_1.sav'

    ;retrieval
    ;pi = replicate({step:0d, fixed:0, limited:[0,0], limits:[0.D,0.D]}, 4)
    ;start = [SP, Albedo, 1d, 0d]

    ;W gaussiun fitting
    pi = replicate({step:0d, fixed:0, limited:[0,0], limits:[0.D,0.D], MPPRINT:[0]}, 5)
    start = [SP, Albedo, Result_Fit0(0), Result_Fit0(1), Result_Fit0(2)]

    err = y_obs*1d-3
    err(*) = median(y_obs)*1d-3

    ;--- added --->
    pi(0).limited(*) = 1
    pi(0).limits(0) = 51
    pi(0).limits(1) = 1499

    pi(1).limited(*) = 1
    pi(1).limits(0) = 0.051
    pi(1).limits(1) = 0.599

    ;simple continuum
    ;pi(2).fixed = 1 

    ;W gaussiun continuum
    pi(2:3).fixed = 1 

    ; iteration imformation print ;0はパラメータ表示させない
    pi(*).MPPRINT = 0

    ;<---

    Result_Fit = MPFITFUN('forward', x, y_obs, err, start, PARINFO=pi, MAXITER=20, BESTNORM=BESTNORM0, MPSIDE=2, status=status, NPRINT=0, yfit=yfit, /nan)
    F = yfit

    ; simple continuum
    ;cont = poly(findgen(29)-float(29/2), [Result_Fit(2:3)])

    ; WG continuum
    cont = pyroxenes(x, Result_Fit0(0:2)) / mean(pyroxenes(x, Result_Fit0(0:2)))

    f(0:3) = -0d/0d
    f(8) = -0d/0d
    f(17) = -0d/0d
    f(24) = -0d/0d
    f(27) = -0d/0d

    bestfit(l,k,*) = F
    continuum(l,k,*) = cont
    pressure(l,k) = Result_Fit(0)
    Albedomap(l,k) = Result_Fit(1)

    end_time = systime(1)
    print, "time:", end_time - start_time
    skip10:

    good = where(FINITE(y_obs) eq 1)
    window,6, xs=800,ys=800
    plot, x(good), y_obs(good), yr=[-0.1,0.4], back=255, color=0, thick=3, psym=-1, xr=[1.85, 2.2], xs=1
    ;oplot, x(good), f(good), color=254, thick=2, psym=-1, linestyle=2
    xyouts, 2.05, -0.08, 'SP='+strcompress(Result_Fit(0)), charsize=2, color=0
    xyouts, 2.05, -0.06, 'Albedo='+strcompress(Result_Fit(1)), charsize=2, color=0
    xyouts, 2.05, -0.04, 'Dust='+strcompress(dust), charsize=2, color=0


    stop

    print, "l LOOP: ", l - ip_b, "/", ip_a - ip_b
    print, "k LOOP: ", k - io_b, "/", io_a - io_b

    ;save, lati, longi, altitude, pressure, Albedomap, dustmap, TAmap, MCDpressure, bestfit, continuum, filename='/work1/LUT/SP/table/absorption/SPmap_ORB0920_3_albedo_WG.sav'

  endfor
endfor


stop
end
