; Arsで計算されたスペクトルを観測データとMCDから補間し、それをさらにinterpolで振るための関数
function ARS_calc_function, LMS0,LMS1,LMS2,LMS3,LMS4,LMS5,LMS6,LMS7,LMS8,LMS9,LMS10,LMS11,LMS12,LMS13,LMS14,LMS15,LMS16,LMS17,LMS18,LMS19,LMS20,LMS21,LMS22,LMS23,LMS24,LMS25,LMS26,input_beta

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

; inputされたpressureによって補間の更新がされていく
wav0 = interpol(Pressure_grid,LMS0,input_beta)
wav1 = interpol(Pressure_grid,LMS1,input_beta)
wav2 = interpol(Pressure_grid,LMS2,input_beta)
wav3 = interpol(Pressure_grid,LMS3,input_beta)
wav4 = interpol(Pressure_grid,LMS4,input_beta)
wav5 = interpol(Pressure_grid,LMS5,input_beta)
wav6 = interpol(Pressure_grid,LMS6,input_beta)
wav7 = interpol(Pressure_grid,LMS7,input_beta)
wav8 = interpol(Pressure_grid,LMS8,input_beta)
wav9 = interpol(Pressure_grid,LMS9,input_beta)
wav10 = interpol(Pressure_grid,LMS10,input_beta)
wav11 = interpol(Pressure_grid,LMS11,input_beta)
wav12 = interpol(Pressure_grid,LMS12,input_beta)
wav13 = interpol(Pressure_grid,LMS13,input_beta)
wav14 = interpol(Pressure_grid,LMS14,input_beta)
wav15 = interpol(Pressure_grid,LMS15,input_beta)
wav16 = interpol(Pressure_grid,LMS16,input_beta)
wav17 = interpol(Pressure_grid,LMS17,input_beta)
wav18 = interpol(Pressure_grid,LMS18,input_beta)
wav19 = interpol(Pressure_grid,LMS19,input_beta)
wav20 = interpol(Pressure_grid,LMS20,input_beta)
wav21 = interpol(Pressure_grid,LMS21,input_beta)
wav22 = interpol(Pressure_grid,LMS22,input_beta)
wav23 = interpol(Pressure_grid,LMS23,input_beta)
wav24 = interpol(Pressure_grid,LMS24,input_beta)
wav25 = interpol(Pressure_grid,LMS25,input_beta)
wav26 = interpol(Pressure_grid,LMS26,input_beta)

; それを1次元データに落とし込む
spectrum = [wav0,wav1,wav2,wav3,wav4,wav5,wav6,wav7,wav8,wav9,wav10,wav11,wav12,wav13,wav14,wav15,wav16,wav17,wav18,wav19,wav20,wav21,wav22,wav23,wav24,wav25,wav26]

return, spectrum
 
end
