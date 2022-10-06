; まずは1点で波長3個分とかでLMS計算を行ってみる

function LMS_calc, TA, TB, SZA, EA, PA, Dust, Waterice, Albedo, wave, Intensity

T1 = TA 
T2 = TB

LMS = dblarr(15)

; for wave = 0,1 do begin
;restore ARS calc data
  if wave eq 0 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave0.sav'
  if wave eq 1 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave1.sav'
  if wave eq 2 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave2.sav'
  if wave eq 3 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave3.sav'
  if wave eq 4 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave4.sav'
  if wave eq 5 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave5.sav'
  if wave eq 6 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave6.sav'
  if wave eq 7 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave7.sav'
  if wave eq 8 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave8.sav'
  if wave eq 9 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave9.sav'
  if wave eq 10 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave10.sav'
  if wave eq 11 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave11.sav'
  if wave eq 12 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave12.sav'
  if wave eq 13 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave13.sav'
  if wave eq 14 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave14.sav'
  if wave eq 15 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave15.sav'
  if wave eq 16 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave16.sav'
  if wave eq 17 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave17.sav'
  if wave eq 18 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave18.sav'
  if wave eq 19 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave19.sav'
  if wave eq 20 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave20.sav'
  if wave eq 21 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave21.sav'
  if wave eq 22 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave22.sav'
  if wave eq 23 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave23.sav'
  if wave eq 24 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave24.sav'
  if wave eq 25 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave25.sav'
  if wave eq 26 then restore,'/work1/LUT/SP/table/LUT_fitting/Table_calc_wave26.sav'


  ;skip
  if T1 gt 285. then begin
    print,'Warning: T1 > limit'
    goto, skip
  endif
  if T1 lt 135. then begin
    print,'Warning: T1 < limit'
    goto, skip
  endif
  if T2 gt 200. then begin
    print,'Warning: T2 > limit'
    goto, skip
  endif
  if T2 lt 80. then begin
    print,'Warning: T2 < limit'
    goto, skip
  endif
  if SZA gt 75. then begin
    print,'Warning: SZA > limit'
    goto, skip
  endif
  if EA gt 50. then begin 
    print,'Warning: EA > limit'
    goto, skip
  endif
  if PA gt 180. then begin
    print,'Warning: PA > limit'
    goto, skip
  endif
  if Dust gt 1.5 then begin
    print,'Warning: DUST > limit'
    goto, skip
  endif
  if Waterice gt 1.0 then begin
    print,'Warning: ICE > limit'
    goto, skip
  endif
  if Albedo gt 0.6 then begin
    print,'Warning: Albedo > limit'
    goto, skip
  endif
  if Albedo lt 0.05 then begin
    print,'Warning: Albedo < limit'
    goto, skip
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

  y = dblarr(15)
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

  for I = 0, 6-1 do begin
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
      t8 = (x8-x8a(j7)) / (x8a(j8+1)-x8a(j8))
      
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

        ;----------------------------------------------------------------------------------
        bad = where(Table_Equivalent_width(j1:j1+1,j2:j2+1,j3:j3+1,j4:j4+1,j5:j5+1,j6:j6+1,j7:j7+1,j8:j8+1) eq 0, count)
        if count ge 1 then Y(i) = -0d/0d
        ;----------------------------------------------------------------------------------

  endfor
   ; print, "y:" ,y

  LMS = (y - Intensity)^2

; endfor

return, LMS

skip:

end


Pro retrieval_pressure_fitting_test
; fitting test

path = '/data2/omega/sav/'
path2 = '/work1/LUT/SP/table/absorption/'
restore, path+'ORB0363_3.sav'
restore, path + 'specmars.sav'

; CO2 absorption emission line   
CO2=where(wvl gt 1.81 and wvl lt 2.19)
wvl=wvl[CO2]

LMS0 = dblarr(15)
LMS1 = dblarr(15)

jdat=jdat(*,CO2,*)
specmars = specmars(CO2)

; ind=32015, xind=15, yind=1000, lati: 48.431797 S, longi: 60.808998 E [Forget+, retrievalすると1036 Pa] ORB0030_1
; ind=where_xyz(longi ge 60.79 and longi le 60.81 and lati ge -48.44 and lati le -48.43,xind=xind,yind=yind)

; ind=76094,xind=62, yind=594, lati:22.705700 N , longi:311.76300 E (48.237 W)  [Forget+, retrievalすると852 Pa] ORB0363_3
ind=where_xyz(longi ge 311.73 and longi le 311.78 and lati ge 22.7 and lati le 22.72,xind=xind,yind=yind)

; ind=37135, xind=15, yind=1160, lati:7.764°S, longi:24.980°E　[Forget+, retrievalすると470 Pa] ORB1201_3
; ind=where_xyz(longi ge 24.979 and longi le 24.981 and lati ge -7.765 and lati le -7.763,xind=xind,yind=yind)

; ind = 68287, xind=63, yind=533, lati:51.068897 N, longi:276.50281 E, 青木さんの結果を再現 ORB0931_3
; ind=where_xyz(longi ge 276.50 and longi le 276.51 and lati ge 51.05 and lati le 51.1,xind=xind,yind=yind)
; ind = where_xyz(longi eq 275.82278 and lati eq 52.079899,xind=xind,yind=yind)


nanserch=where(jdat ge 0 and jdat le 0.0001)
jdat(nanserch)=!VALUES.F_NAN

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

x = [wvl(0),wvl(1),wvl(2),wvl(23),wvl(24),wvl(25)]
Y = [jdat(xind,0,yind), jdat(xind,1,yind), jdat(xind,2,yind), jdat(xind,23,yind), jdat(xind,24,yind), jdat(xind,25,yind)]
coef = linfit(X,Y)
cont = coef(0) + coef(1)*wvl

SZA = reform(geocube(xind,8,yind))*1.e-4
EA = reform(geocube(xind,9,yind))*1.e-4
PA = reform(geocube(xind,10,yind))*1.e-4

Albedo_input = jdat(xind,0,yind)/specmars(0) / cos(geocube(xind,8,yind)*1e-4*!DTOR)

width = 1.0 - jdat(xind,*,yind)/cont
Intensity = width

for wave = 0,26 do begin

  if wave eq 0 then LMS0 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 1 then LMS1 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 2 then LMS2 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 3 then LMS3 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 4 then LMS4 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 5 then LMS5 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 6 then LMS6 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 7 then LMS7 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 8 then LMS8 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 9 then LMS9 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 10 then LMS10 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 11 then LMS11 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 12 then LMS12 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 13 then LMS13 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 14 then LMS14 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 15 then LMS15 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 16 then LMS16 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 17 then LMS17 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 18 then LMS18 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 19 then LMS19 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 20 then LMS20 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 21 then LMS21 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 22 then LMS22 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 23 then LMS23 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 24 then LMS24 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 25 then LMS25 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)
  if wave eq 26 then LMS26 = LMS_calc(TA, TB, SZA, EA, PA, dust_opacity, ice_opacity, Albedo_input,wave,Intensity)

endfor

; (1) 最小二乗和を計算する
; total_LMS = LMS0 + LMS1 + LMS2 + LMS3 + LMS4 + LMS5 + LMS6 + LMS7 + LMS8 + LMS9 + LMS10 + LMS11 + LMS12 + LMS13 + LMS14 + LMS15  + LMS17 + LMS18 + LMS19 + LMS20 + LMS21 + LMS22 + LMS23 + LMS24 + LMS25 + LMS26
total_LMS = LMS9 + LMS10 + LMS11 + LMS12 + LMS13 + LMS14 + LMS15

; (2) 最小二乗和の中で最小値を探す
min_LMS = min(total_LMS)

; debug →→

; min_ind = (where(total_LMS eq min_LMS)) + 5
; print, "total_LMS", total_LMS
; print, "total_min", min_LMS
; print, "min_ind", min_ind
; help, total_LMS
; total_LMS_1 = LMS0
; total_LMS_2 = LMS26
; total_LMS_3 = LMS23
; total_LMS_4 = LMS24
; total_LMS_5 = LMS25
; print, "LMS0:", LMS0(0)
; print, "LMS26:", LMS26(0)
; print, "total:", total_LMS
; print, "before_pres", exp(before_pres)
; print, "after_pres", exp(after_pres)
; print, "before:", before_interpol
; print, "after:", after_interpol
; print, "result:", result_compare

; ←← debug

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

; TBD! indの与え方が違う。初期値とかの部分を変えないといけない。あとは、
; (3) 極小値を探す
; 最小二乗和の中で最小値を持つindex前後の値を持ってきて、その中間地点での値を出す
; repeat begin
min_ind = where(total_LMS eq min_LMS) + 8b
print, "intial:" ,min_ind

; repeat begin
for n=0, 5 do begin
  pup = 1d
  pdown = -1d
  d = (1/2d)^n

  before_pres = (Pressure_grid(min_ind - d) + Pressure_grid(min_ind))/2
  after_pres = (Pressure_grid(min_ind + d) + Pressure_grid(min_ind))/2

; 最小値からどっちが最小二乗和が小さいかを判断する
  before_interpol = interpol(total_LMS, Pressure_grid, before_pres)
  after_interpol = interpol(total_LMS, Pressure_grid, after_pres)

  result_compare = before_interpol < after_interpol

  if before_interpol eq result_compare then begin
    min_ind =  min_ind + (pdown*(1/2d)^(n+1))
    print, "if:" , min_ind
  endif else  begin
    min_ind =  min_ind + (pup*(1/2d)^(n+1))
    print, "else:" ,min_ind
  Endelse

  ; print, "min_ind:", min_ind
  ; print, "回数" , n

endfor

print, "LMS_result",result_compare

pressure = interpol(Pressure_grid, total_LMS, result_compare)
print, "pressure:", exp(pressure)

; endrep until result_compare lt 0.0001
stop

repeat begin
if before_interpol eq result_compare then begin

  loop_ind_1 = (Pressure_grid(min_ind - 1) + before_pres)/2
  loop_ind_2 = (Pressure_grid(min_ind) + before_pres)/2

  print, "loop1:", loop_ind_1
  print, "loop2:", loop_ind_2

  before_interpol = interpol(total_LMS, Pressure_grid, loop_ind_1)
  after_interpol = interpol(total_LMS, Pressure_grid, loop_ind_2)

  result_compare = before_interpol < after_interpol
  print, "loop_result:", result_compare

endif else  begin

  loop_ind_1 = (Pressure_grid(min_ind + 1) + after_pres)/2
  loop_ind_2 = (Pressure_grid(min_ind) + after_pres)/2

  print, "loop1:", loop_ind_1
  print, "loop2:", loop_ind_2

  before_interpol = interpol(total_LMS, Pressure_grid, loop_ind_1)
  after_interpol = interpol(total_LMS, Pressure_grid, loop_ind_2)

  result_compare = before_interpol < after_interpol
  print, "loop_result:", result_compare

Endelse
endrep until result_compare gt 0.0001


; *** TBD! ***
; pressure1 = interpol(Pressure_grid,total_LMS_1,Intensity[0])
; pressure2 = interpol(Pressure_grid,total_LMS_2,Intensity[26])
; pressure3 = interpol(Pressure_grid,total_LMS_3,Intensity[23])
; pressure4 = interpol(Pressure_grid,total_LMS_4,Intensity[24])
; pressure5 = interpol(Pressure_grid,total_LMS_5,Intensity[25])

; ORB0363_3でチェック中 retrievalすると852 Pa
; wn(0):5.9526977e-228      
; wn(1):6.0212888e+92
; wn(2):0.0000000
; wn(3):0.0000000
; wn(4):22033.758
; wn(5):2239.2206
; wn(6):89674.872      
; wn(7):7.3313055e+11
; wn(8):5.6501814e+13
; wn(9):1858.5543
; wn(10):1027.2930
; wn(11):1217.2015      
; wn(12):8661.1304
; wn(13):1326.0085
; wn(14):1078.4356
; wn(15):1380.4640
; wn(16):NaN   
; wn(17):1359.7269
; wn(18):1156.3858
; wn(19):1411.8142
; wn(20):2294.0319
; wn(21):4820.2817 
; wn(22):5.5723976e+10
; wn(23):265.76018
; wn(24):6.5631814e-13
; wn(25):1.8564746e-06
; wn(26):NaN
 
end
