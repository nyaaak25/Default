; ARSのデータを観測データとMCDから補間するfunction
function ARS_calc, TA, TB, SZA, EA, PA, Dust, Waterice, Albedo, wave, Intensity

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

; endfor

return, y

skip:

end
