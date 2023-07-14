;------------------------------------------------------------------------------------------------------------------------
;
; Equivalent width methodを用いて気圧mapを作成するプログラム
;
; create by Akira Kazama
; 
; retrieval_pres_map_ver2　 ::2022.10.31 Mon 12:08:00
; ver2: LUT update後に対応するプログラム
;
; retrieval_pres_map_ver3　 ::2023.02.21 Tue 16:21:00
; ver3: 等価幅法でQLを作成するためのプログラム
; add -> 取り敢えずは、大まかな大きな分布が全部ほしいから、観測条件が良いモノだけをpick upさせている。
; limited: Albedo, EA, SZA
;
; retrieval_pres_map_ver4　 ::2023.03.01 Wed 18:43:00
; ver4: 等価幅法でQLを作成するためのプログラム
; add -> 条件が悪いところを除いて、全部のORBに適応するようにプログラムを改変する
; limited: EA<30, SZA<55, 極域(70°以南/以北)、LUTの軸を超えた値ならNaNが入るようにする
; *TBD* dust factorを追加、T1を追加
;
; retrieval_pres_map_ver5　 ::2023.04.12 Wed 11:17:00
; ver5: 等価幅法でQLを作成するためのプログラム
; add -> 条件が悪いところを除いて、全部のORBに適応するようにプログラムを改変する
; 信頼性の低い素子をORBごとに変更する仕様を追加
; !TBD! MY28以前/以降も適用可能に変更
; 
;------------------------------------------------------------------------------------------------------------------------

;------------------------------------------------------------------------------------------------------------------------
function ret_pressure, trans, TA, TB, SZA, EA, PA, Dust, Waterice, Albedo, fileorbit
;------------------------------------------------------------------------------------------------------------------------

T1 = TA 
T2 = TB

pressure_CD = 0d

; ここでorbごとにrestoreするfileを変更する仕様を追加
; restore, into equivalent width table
if fileorbit le 2000 then begin 
  restore,'/work1/LUT/SP/table/absorption/density/Table_SP_calc_ver2_orb0000.sav'
endif

if fileorbit gt 2000 and fileorbit le 2431 then begin
  restore,'/work1/LUT/SP/table/absorption/density/Table_SP_calc_ver2_orb2000.sav'
endif

if fileorbit gt 2431 and fileorbit le 2618 then begin
  restore,'/work1/LUT/SP/table/absorption/density/Table_SP_calc_ver2_orb2432.sav'
endif

if fileorbit gt 2618 and fileorbit le 5306 then begin
  restore,'/work1/LUT/SP/table/absorption/density/Table_SP_calc_ver2_orb2618.sav'
endif

if fileorbit gt 5306 and fileorbit le 5692 then begin
  restore,'/work1/LUT/SP/table/absorption/density/Table_SP_calc_ver2_orb5306.sav'
endif

if fileorbit gt 5692 and fileorbit le 6008 then begin
  restore,'/work1/LUT/SP/table/absorption/density/Table_SP_calc_ver2_orb5692.sav'
endif

if fileorbit gt 6008 then begin
  restore,'/work1/LUT/SP/table/absorption/density/Table_SP_calc_ver2_orb6008.sav'
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

y = dblarr(15)
X1 = (T1)^1.5d 
X2 = (T2)^1.5d 
x3 = exp(-cos(SZA/180d*!dpi))
x4 = exp(-cos(EA/180d*!dpi))
x5 = -cos(PA/180d*!dpi) + 1.d
x6 = Dust
x7 = Waterice
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

for I = 0, 5-1 do begin
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

for I = 0, 15-1 do begin

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
  
  y(i) = y11111111*(1.d - t1)*(1.d - t2)*(1.d - t3)*(1.d - t4)*(1.d - t5)*(1.d - t6)*(1.d - t7)*(1.d - t8) +   $
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

pressure_CD=interpol(Pressure_grid,Y,trans)
return, pressure_CD

end



Pro retrieval_pres_map_ver5

; =========== path ===========
pathfile = '/data2/omega/sav/'
path_work = '/work1/LUT/SP/EWwork/'
path_ql = '/work1/LUT/SP/QL/'

files=file_search(pathfile+'*.sav',count=count) ; sav fileすべてを持ってくる

; どの軌道番号から始めたいかを指定する
;loop_b = where(files eq '/data2/omega/sav/ORB0920_3.sav')
loop_b = where(files eq '/data2/omega/sav/ORB0923_3.sav')

; どの軌道番号まで回したいかを指定する
loop_a = where(files eq '/data2/omega/sav/ORB7697_3.sav')


; ======== > 参考
; MY27
;loop_b = where(files eq '/data2/omega/sav/ORB0181_0.sav') ; MY27はじまりのORB
;loop_a = where(files eq '/data2/omega/sav/ORB2595_0.sav') ; MY27おわりのORB

; MY28
;loop_b = where(files eq '/data2/omega/sav/ORB2607_0.sav') ; MY28はじまりのORB
;loop_a = where(files eq '/data2/omega/sav/ORB5055_5.sav') ; MY28おわりのORB

; MY29
;loop_b = where(files eq '/data2/omega/sav/ORB5062_0.sav') ; MY29はじまりのORB
;loop_a = where(files eq '/data2/omega/sav/ORB7454_4.sav') ; MY29おわりのORB

; < ======== 

loop_b = fix(loop_b(0))
loop_a = fix(loop_a(0))

; === can not restore file memo ===
; ORB1991_3.sav, ORB2079_6.sav, ORB2138_5-7.sav, ORB2139_0.sav, ORB2141_0.sav, ORB2456_0.sav, ORB2533_0.sav,  
; ORB3194_4.sav, ORB3354_3.sav, ORB3354_4.sav, ORB4254_4-6.sav, ORB4260_0-6.sav, ORB4262_0-6.sav, ORB4265_0-6.sav, ORB4266_0-6.sav
; ORB4267_0-6.sav, ORB4270_0-2.sav, ORB4287_5.sav
; ORB6602_3.sav, ORB7393_0.sav, ORB7401_0.sav, ORB7407_1-4.sav, ORB7410_0-2.sav


; =========== file loop start ===========
for loop = loop_b, loop_b do begin
  start_time = systime(1)
  
; ++++++ restore file ++++++
  ; retrieval するfileをrestore
  file = files(loop)
  restore, file

  ; getting orb file name
  fileorbit = strupcase(strmid(file,12,9,/REVERSE_OFFSET))
  input_fileorbit = strupcase(strmid(file,9,4,/REVERSE_OFFSET))

; ++++++++ retrieval準備 ++++++++++
  ; 全部indexを持ってくる
  ip = n_elements(LATI(*,0))
  io = n_elements(LATI(0,*))
  
  ; CO2スペクトル範囲を選択
  CO2=where(wvl gt 1.8 and wvl lt 2.2) ; 波長範囲はここで変更可能
  wvl=wvl[CO2]
  jdat_ind = jdat(*,CO2,*)
  jdat=jdat(*,CO2,*)

  ; search for broken spectrel [下参照]
  ; 経年劣化分を探す
  nan = where(jdat le 0.0000001)
  jdat(nan) = !VALUES.F_NAN

  ; SOFT10に記載されている場所を探す
  ; 2 um吸収帯(1.8から2.2 um特論)
  bad1 = (78 - CO2(0))
  bad2 = (69 - CO2(0))
  bad3 = (88 - CO2(0))
  bad4 = (66 - CO2(0))
  bad5 = (79 - CO2(0))
  bad6 = (85 - CO2(0))

  input_bad = [bad1, bad2, bad3, bad4, bad5, bad6]

  jdat(*,input_bad(0:2),*)= !VALUES.F_NAN
  ; ORB2000から使用しないスペクトル
  if input_fileorbit ge 2000 then jdat(*,input_bad(3:5),*)= !VALUES.F_NAN

  ; ========= > reference@OMEGA SOFT10
  ; 使用しないスペクトル

  ; since beginning
  ; 1.4142500(34)   2.0392201(78)   3.1515501(158)
  ; 1.9131700(69)   2.1776299(88)   4.4997802(224)

  ; since orb1147
  ; 3.7730401(188)

  ; since orb1990
  ; 3.0893900(155)

  ; since orb2000
  ; 1.7149000(55)   1.8708900(66)   2.0531399(79)   2.1363101(85)   2.6180999(121)   2.6948700(127)   4.0179400(200)   4.4563298(222)
  ; < ==========


  ; 太陽輝度fileをrestore
  restore, pathfile + 'specmars.sav'
  specmars = specmars(CO2)

  ; continuumの波長点を選択
  if input_fileorbit le 5306 then begin
    x = [wvl(0),wvl(1),wvl(2),wvl(23),wvl(25),wvl(26)]
  endif

  if input_fileorbit gt 5306 and input_fileorbit le 6008 then begin
    x = [wvl(0),wvl(2),wvl(23),wvl(25),wvl(26)]
  endif

  if input_fileorbit gt 6008 then begin
    x = [wvl(0),wvl(23),wvl(25),wvl(26)]
  endif

; +++++++ 使用するband幅を指定 ++++++++
  ; band幅2: 全CO2吸収帯を使用
  band = where(wvl gt 1.94 and wvl lt 2.09)

  ; band幅3: ある1つのCO2吸収帯を使用
  ;band = where(wvl gt 1.94 and wvl lt 1.99)

; +++++++ create array +++++++
  pressure = dblarr(ip,io)
  trans = dblarr(ip,io)
  MCDpressure = dblarr(ip,io)
  dustmap = dblarr(ip,io)
  TAmap = dblarr(ip,io)
  altitude = dblarr(ip,io)
  InputAlbedo = dblarr(ip,io)

  ; observation geometry save
  SZA_all = reform(geocube(*,8,*))*1.e-4
  EA_all = reform(geocube(*,9,*))*1.e-4
  PA_all = reform(geocube(*,10,*))*1.e-4

; +++++++ MCD v5.3 data input file +++++++
  ; local timeの計算
  timei=reform(geocube(0:6,1,*))
  time=timei(*,*)
  local_time=(longi-sub_solar_longitude)*24/360+12
  indLT1=where(local_time lt 0)
  if n_elements(indLT1) gt 1 then local_time(indLT1)=24+local_time(indLT1)
  indLT2=where(local_time ge 24)
  if n_elements(indLT2) gt 1 then local_time(indLT2)=local_time(indLT2)-24

  ; ========== retrieval loop start ========== 
  ;for l = 0, ip -1 do begin ;loop for slit scan
  ;  for k = 0, io -1 do begin ;test getting surface feature

  ; ========== retrieval loop start ========== 
  for l = 1, 1 do begin ;loop for slit scan
    for k = 71, 71 do begin ;test getting surface feature

; =========== bad spectrum skip =========== 
      ; 29個の素子の中で、放射輝度の値がおかしいところを検索
      notgood = where(jdat_ind(l,*,k) lt 0.000001)
      notgood_ind = n_elements(notgood)

      if input_fileorbit lt 2432 then begin
        if notgood_ind gt 2 then begin
          ; print, 'BAD INDEX: ', notgood_ind
          ; print, 'WARNING: This spectrum is bad data'
          goto, skip20
        endif
      endif

      if input_fileorbit ge 2432 and input_fileorbit lt 2618 then begin
        if notgood_ind gt 3 then begin
          ; print, 'BAD INDEX: ', notgood_ind
          ; print, 'WARNING: This spectrum is bad data'
          goto, skip20
        endif
      endif

      if input_fileorbit ge 2618  and input_fileorbit lt 5306 then begin
        if notgood_ind gt 4 then begin
          ; print, 'BAD INDEX: ', notgood_ind
          ; print, 'WARNING: This spectrum is bad data'
          goto, skip20
        endif
      endif

      if input_fileorbit ge 5306 and input_fileorbit lt 5408 then begin
        if notgood_ind gt 5 then begin
          ; print, 'BAD INDEX: ', notgood_ind
          ; print, 'WARNING: This spectrum is bad data'
          goto, skip20
        endif
      endif

      if input_fileorbit ge 5408 and input_fileorbit lt 5692 then begin
        if notgood_ind gt 6 then begin
          ; print, 'BAD INDEX: ', notgood_ind
          ; print, 'WARNING: This spectrum is bad data'
          goto, skip20
        endif
      endif

      if input_fileorbit ge 5692 and input_fileorbit lt 6008 then begin
        if notgood_ind gt 7 then begin
          ; print, 'BAD INDEX: ', notgood_ind
          ; print, 'WARNING: This spectrum is bad data'
          goto, skip20
        endif
      endif

      if input_fileorbit ge 6008 then begin
        if notgood_ind gt 9 then begin
          ; print, 'BAD INDEX: ', notgood_ind
          ; print, 'WARNING: This spectrum is bad data'
          goto, skip20
        endif
      endif

      ; define geometry
      lat = lati(l,k)
      lon = longi(l,k)
      SZA = reform(geocube(l,8,k))*1.e-4
      EA = reform(geocube(l,9,k))*1.e-4
      PA = reform(geocube(l,10,k))*1.e-4
      altitude(l,k) = reform(geocube(l,12,k))

      ; reflectance factor ( I/F / cos(SZA) )
      Albedo_input = jdat(l,0,k)/specmars(0) / cos(geocube(l,8,k)*1e-4*!DTOR)

; =========== retrieval skip ===========
      ; SZA limited: cos(SZA)<0.6 ~ 55
      if SZA gt 55. then begin
        ;print, 'SZA Warning:', SZA
        goto, skip20
      endif
      ; EA limited: EA < 30° (only nadir obeservation)
      if EA gt 30. then begin
        ;print, 'EA Warning:', EA
        goto, skip20
      endif

      ; geometry limited: lati > ±70° (exception polar region)
      if lat gt 70. then begin 
        ;print, 'latitude Warning:', lat
        goto, skip20
      endif

      if lat lt -70. then begin 
        ;print, 'latitude Warning:', lat
        goto, skip20
      endif

      ; =========== MCD version 6.1 ===========
      Ls = SOLAR_LONGITUDE
      Loct = local_time
      hrkey = 1 ; set high resolution mode on (hrkey=0 to set high resolution off)
      zkey = 3    ; specify that xz is the altitude above surface (m)
      xz = 0. ; (array of altitude: step=2000 m start: 5m)
      datekey = 1       ; <integer> type of input date (1=Mars date)
      xdate = Ls       ; <double precision> date (IF datekey = 1 : Value of Ls [deg.])
      dset ='/work1/LUT/MCDv6.1/MCD_6.1_queleclim/data/'
      scena = 27         ; <integer> scenario (1 = Climatology ave solar)
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

      cd, '/work1/LUT/MCDv6.1/MCD_6.1_queleclim/mcd/interfaces/idl/'
      a = call_mcd(zkey,xz,lon,lat,hrkey,datekey,xdate,Loct, $
        dset,scena,perturkey,seedin,gwlength,extvarkeys, $
        pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)

      scaleH = extvar(63)
      Z1 = scaleH * 0.1d
      Z2 = scaleH * 4d
      dust_opacity = extvar(40)*1.2090217E+00/4.6232791E+00
      ice_opacity = 0.d;*4.2410289E-01/***  ;!TBD!
      SP = pres

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

; +++++++++ save MCD data +++++++++++
      InputAlbedo(l,k) = Albedo_input
      MCDpressure(l,k) = SP
      dustmap(l,k) = dust_opacity
      TAmap(l,k) = TA

; =========== retrieval skip ===========
      ; LUT parameter limited
      if PA ge 180. then begin
        ;print, 'PA Warning:', PA
        goto, skip20
      endif

      if TA ge 285. then begin
        ;print, 'TA Warning:', TA
        goto, skip20
      endif

      if TA le 135. then begin
        ;print, 'TA Warning:', TA
        goto, skip20
      endif

      if TB ge 200. then begin
        ;print, 'TB Warning:', TB
        goto, skip20
      endif

      if TB le 80. then begin
        ;print, 'TB Warning:', TB
        goto, skip20
      endif

      if dust_opacity ge 1.5 then begin
        print, 'Dust Warning:', dust_opacity
        goto, skip20
      endif

      if Albedo_input le 0.05 then begin
        print, 'Albedo Warning:', Albedo_input
        goto, skip20
      endif

      if Albedo_input ge 0.6 then begin
        print, 'Albedo Warning:', Albedo_input
        goto, skip20
      endif

; =========== retrieval start ===========
      if input_fileorbit le 5306 then begin
        Y = [jdat(l,0,k), jdat(l,1,k), jdat(l,2,k), jdat(l,23,k), jdat(l,25,k), jdat(l,26,k)]
      endif

      if input_fileorbit gt 5306 and input_fileorbit le 6008 then begin
        Y = [jdat(l,0,k), jdat(l,2,k), jdat(l,23,k), jdat(l,25,k), jdat(l,26,k)]
      endif

      if input_fileorbit gt 6008 then begin
        Y = [jdat(l,0,k), jdat(l,23,k), jdat(l,25,k), jdat(l,26,k)]
      endif

      coef = linfit(X,Y)
      cont = coef(0) + coef(1)*wvl

; +++++++++ Equivalent width計算 ++++++++++++
      width = 1.0 - jdat(l,*,k)/cont
      trans(l,k) = (total(width[band], /nan))

; +++++++++ retrieval ++++++++++++
      pressure(l,k) = ret_pressure(trans(l,k), TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input, input_fileorbit)
      
      skip20:
    
    endfor
    if l mod 5 eq 0 then print, 'LOOP: ', l, '/', ip
  endfor

;  +++++++++ 各種indexを計算 ++++++++++++
  ; retrievalされた圧力値から地形情報除去
  g = 3.72 ;m/s2
  Gconst = 6.67430d-11
  MMars = 6.42d23
  RMars = 3.4d6
  g = -Gconst*MMars/(-RMars*RMars)
  goodT = where(TAmap gt 0)
  R = 192d ;J/K/kg

  p1_slev = pressure
  p1 = pressure
  p1 = exp(p1)
  for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) eq 1 then p1(i,j) = !VALUES.F_NAN
  
  p5 = altitude
  for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 0 then p1_slev(i,j) = p1(i,j) * exp( (p5(i,j)) / (R*TAmap(i,j)/(-Gconst*MMars/(-1d*(RMars+p5(i,J))*(RMars+p5(i,J)) ))))

  ; =========== 使えるデータ数を記入 ===========
  ; bad data skip
  bad_skip = n_elements(nan)
  bad_ratio = bad_skip *100d / (ip * io * 29)
  ; 素子29個中平均何個がbad dataかわかる指標
  bad_index = bad_skip / (ip * io)

  ; sza skip
  where_SZA = where(SZA_all gt 55)
  SZA_skip = n_elements(where_SZA)
  SZA_ratio = SZA_skip *100d / (ip * io)
  
  ; albedo low skip
  where_albedo_low = where(InputAlbedo lt 0.05)
  albedo_low_skip = n_elements(where_albedo_low)
  albedo_low_ratio = albedo_low_skip *100d / (ip * io)

  ; albedo high skip
  where_albedo_high = where(InputAlbedo lt 0.6)
  albedo_high_skip = n_elements(where_albedo_high)
  albedo_high_ratio = albedo_high_skip *100d / (ip * io)

  ; EA skip
  where_EA = where(EA_all gt 30)
  EA_skip = n_elements(where_EA)
  EA_ratio = EA_skip *100d / (ip * io)

  ; dust skip
  where_dust = where(dustmap gt 0.4d)
  dust_skip = n_elements(where_dust)
  dust_ratio = dust_skip *100d / (ip * io)

  print, 'bad data ratio: ', bad_ratio
  print, 'bad sza ratio: ', SZA_ratio

  ; =========== save work file ===========
  ; create work file
  save, lati, longi, altitude, trans, pressure, InputAlbedo, dustmap, TAmap, MCDpressure, SZA_all, EA_all, PA_all, local_time, SOLAR_LONGITUDE, p1_slev, bad_index, bad_ratio, SZA_ratio, EA_ratio, dust_ratio, albedo_high_ratio, albedo_low_skip, filename = path_work + 'EW_work_' + fileorbit + '.sav'

  ; なんらかのskipに引っかかってimageの97%以上Retrievalされなかった場合は、そのQLを作らない
  skip_ind = where(pressure eq 0)
  n_skip = n_elements(skip_ind)
  skip_ratio = n_skip *100d / (ip * io)
  if skip_ratio gt 97d  then begin
    print, 'ERROR: This orbit could not retrieval. "' + fileorbit + '"'
    print, 'WARNING: Could not create Quick Look'
    goto, skip10
  endif

  ; 軌道情報 other parameters
  localtime = strmid(mean(local_time),5,8)
  Ls=strmid(SOLAR_LONGITUDE,5,7)
  Latimean=strmid(mean(lati),5,6)
  Longimean=strmid(mean(longi),5,6)
  year=strmid(time(0),8,4)
  month=strmid(time(1),10,2)
  day=strmid(time(2),10,2)
  hour=strmid(time(3),10,2)
  minit=strmid(time(4),10,2)
  

  ; =========== plot pressure sav ===========
  device, decomposed = 0, retain = 2
  loadct,39

  lat1 = lati
  minlati = min(lati)
  maxlati = max(lati)
  lon1 = longi
  minlon = min(longi)
  maxlon = max(longi)

  ; retrieval pressure
  minp1 = min(p1)
  maxp1 = max(p1)
  medip1 = median(p1)

  ; MCD pressure
  p2 = MCDpressure
  for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) eq 0d then p2(i,j) = !VALUES.F_NAN
  minp2 = min(p2)
  maxp2 = max(p2)

  ; derivetion pressure
  p3 = p1 - p2
  medip3 = median(p3)

  ; 何%変動までみたいかによってcolor rangeを変更
  ; 現在は5%変動までをみられるように設定
  dev_range = 0.05*medip1
  minp3 = medip3 - dev_range
  maxp3 = medip3 + dev_range

  ; dust opacity
  p4 = dustmap
  for i = 0, ip-1 do for j = 0, io-1 do if p4(i,j) eq 0d then p4(i,j) = !VALUES.F_NAN
  minp4 = min(p4)
  maxp4 = max(p4)

  ; MOLA altitude
  minp5 = min(p5)
  maxp5 = max(p5)

  ; albedo input
  p6 = inputalbedo
  for i = 0, ip-1 do for j = 0, io-1 do if p6(i,j) eq 0d then p6(i,j) = !VALUES.F_NAN
  minp6 = min(p6)
  maxp6 = max(p6)

  ; topography remove
  p8 = p1_slev
  minp8 = min(p8)
  maxp8 = max(p8)


  ; プロットの枠組み
  A = FINDGEN(17) * (!PI*2/16.)
  USERSYM, COS(A), SIN(A), /FILL

  window, 22, xs=1800, ys=1600
  !P.Multi = [0, 3, 2]

  ; Retrieved pressure plot
  loadct,39
  plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.55,0.31,0.92], /nodata, charsize=3, title='Retrieved pressure', xtitle='Lon', ytitle='Lat'
  cgLOADCT, 16
  cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp1,maxp1],position=[0.04,0.58,0.046,0.86],title='Pa', TCHARSIZE=10, /vertical 
  loadct, 16
  for i = 0, ip-1 do for j = 0, io-1 do if p1(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p1(i,j)-minp1)/(maxp1-minp1)*254., psym=8, symsize=1

  ; MCD prediction pressure
  loadct,39
  plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.55,0.65,0.92], /nodata, charsize=3, title='MCD pressure', xtitle='Lon', ytitle='Lat'
  cgLOADCT, 16
  cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp2,maxp2],position=[0.38,0.58,0.386,0.86],title='Pa', TCHARSIZE=10, /vertical 
  loadct, 16
  for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p2(i,j)-minp2)/(maxp2-minp2)*254., psym=8, symsize=1

  ; retrieved - MCD
  loadct,39
  plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.55,0.99,0.92], /nodata, charsize=3, title='Retrieval - MCD', xtitle='Lon', ytitle='Lat'
  cgLOADCT, 39
  cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp3,maxp3],position=[0.72,0.58,0.726,0.86],title='Pa', TCHARSIZE=10, /vertical 
  for i = 0, ip-1 do for j = 0, io-1 do if p2(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p3(i,j)-minp3)/(maxp3-minp3)*254., psym=8, symsize=1

  ; input Albedo map 
  loadct,39
  plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.1,0.08,0.31,0.45], /nodata, charsize=3, title='Albedo map', xtitle='Lon', ytitle='Lat'
  cgLOADCT, 16
  cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp6,maxp6],position=[0.04,0.12,0.046,0.4],title='Albedo', TCHARSIZE=10, /vertical 
  loadct, 16
  for i = 0, ip-1 do for j = 0, io-1 do if p6(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p6(i,j)-minp6)/(maxp6-minp6)*254., psym=8, symsize=1

  ; dust opacity map
  loadct,39
  plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.44,0.08,0.65,0.45], /nodata, charsize=3, title='dust opacity', xtitle='Lon', ytitle='Lat'
  cgLOADCT, 16
  cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp4,maxp4],position=[0.38,0.12,0.386,0.4],title='Dust opacity', TCHARSIZE=10, /vertical 
  loadct, 16
  for i = 0, ip-1 do for j = 0, io-1 do if p4(i,j) gt 0 then plots, lon1(i,j), lat1(i,j), color=(p4(i,j)-minp4)/(maxp4-minp4)*254., psym=8, symsize=1

  ; MOLA altitude map
  loadct,39
  plot, findgen(10),xs=1, ys=1, yr=[minlati, maxlati], xr=[minlon, maxlon], back=255, color=0, position=[0.78,0.08,0.99,0.45], /nodata, charsize=3, title='MOLA altitude', xtitle='Lon', ytitle='Lat'
  cgLOADCT, 16
  cgColorbar, Divisions=4, Minor=5, charsize=2, Range=[minp5,maxp5],position=[0.72,0.12,0.726,0.4],title='m', TCHARSIZE=10, /vertical 
  loadct, 16
  for i = 0, ip-1 do for j = 0, io-1 do plots, lon1(i,j), lat1(i,j), color=(p5(i,j)-minp5)/(maxp5-minp5)*254., psym=8, symsize=1

  loadct, 39
  xyouts,0.06,0.975,fileorbit,charsize=2.5,color=0,/normal
  xyouts,0.2,0.975,'Time : '+year+'/'+month+'/'+day+' '+hour+':'+minit,color=0,/normal,charsize=2
  xyouts,0.45,0.975,'Ls : '+Ls,color=0,/normal,charsize=2
  xyouts,0.57,0.975,'Local Time : '+localtime,color=0,/normal,charsize=2

  snapshot = TVRD(True=1)
  Write_JPEG, path_ql + 'SP_' + fileorbit +'_test.jpg', snapshot, True=1, Quality=75
  
  skip10:
  print, 'ORB LOOP: ', (loop - loop_b) + 1 , '/', (loop_a - loop_b) + 1
  print, 'ORB number: ', fileorbit

  end_time = systime(1)
  print, 'time: ', (end_time - start_time)/60 , 'min'

endfor
stop

end
