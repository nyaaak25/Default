;------------------------------------------------------------------------------------------------------------------------
Pro Cal_LUtable_ver2_nyonn
; +++
; to read outputs from ARS and calculate equivalent width (=absorption depth) for each spectrum
; create by Shohei Aoki
; 
; edit by Akira Kazama
; 
; cal_lutable_ver2_nyonn　 ::2022.6.1 Wed 10:49:00
; ver2: 吸収の深さ量を計算するプログラム
; 
; +++
;------------------------------------------------------------------------------------------------------------------------
path1 = '/work1/LUT/SP/table/output_SP1/'
path2 = '/work1/LUT/SP/table/output_SP2/'
path3 = '/work1/LUT/SP/table/output_SP3/'
path4 = '/work1/LUT/SP/table/output_SP4/'
path5 = '/work1/LUT/SP/table/output_SP5/'
path6 = '/work1/LUT/SP/table/output_SP6/'
path7 = '/work1/LUT/SP/table/output_SP7/'
path8 = '/work1/LUT/SP/table/output_SP8/'
path9 = '/work1/LUT/SP/table/output_SP9/'
path10 = '/work1/LUT/SP/table/output_SP10/'
path11 = '/work1/LUT/SP/table/output_SP11/'
path12 = '/work1/LUT/SP/table/output_SP12/'
path13 = '/work1/LUT/SP/table/output_SP13/'
path14 = '/work1/LUT/SP/table/output_SP14/'
path15 = '/work1/LUT/SP/table/output_SP15/'

Table_Equivalent_pressure1=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure2=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure3=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure4=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure5=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure6=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure7=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure8=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure9=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure10=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure11=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure12=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure13=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure14=dblarr(5,3,6,5,5,6,3,7)
Table_Equivalent_pressure15=dblarr(5,3,6,5,5,6,3,7)

wn=dblarr(27)
rad1=dblarr(27)
rad2=dblarr(27)
rad3=dblarr(27)
rad4=dblarr(27)
rad5=dblarr(27)
rad6=dblarr(27)
rad7=dblarr(27)
rad8=dblarr(27)
rad9=dblarr(27)
rad10=dblarr(27)
rad11=dblarr(27)
rad12=dblarr(27)
rad13=dblarr(27)
rad14=dblarr(27)
rad15=dblarr(27)

x = dblarr(6)
y1 = dblarr(6)
y2 = dblarr(6)
y3 = dblarr(6)
y4 = dblarr(6)
y5 = dblarr(6)
y6 = dblarr(6)
y7 = dblarr(6)
y8 = dblarr(6)
y9 = dblarr(6)
y10 = dblarr(6)
y11 = dblarr(6)
y12 = dblarr(6)
y13 = dblarr(6)
y14 = dblarr(6)
y15 = dblarr(6)

for IT1 = 1, 5 do begin
  for IT2 = 1, 3 do begin
    for ISZA = 1, 6 do begin
      for IEA = 1, 5 do begin
        for IPA = 1, 5 do begin
          for ID = 1, 6 do begin
            for IWI = 1, 3 do begin
              for IAB = 1, 7 do begin
                
                ;file names
                file1 = path1 + 'SP1' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file2 = path2 + 'SP2' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file3 = path3 + 'SP3' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                
                file4 = path4 + 'SP4' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file5 = path5 + 'SP5' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file6 = path6 + 'SP6' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'    
                
                file7 = path7 + 'SP7' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file8 = path8 + 'SP8' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file9 = path9 + 'SP9' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file10 = path10 + 'SP10' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file11 = path11 + 'SP11' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file12 = path12 + 'SP12' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file13 = path13 + 'SP13' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat' 
                                     
                file14 = path14 + 'SP14' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'                                                                                                            

                file15 = path15 + 'SP15' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'


                ;read files
                openr,lun,file1,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  wn(i)=a
                  rad1(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file2,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad2(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file3,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad3(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file4,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad4(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file5,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad5(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file6,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad6(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file7,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad7(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file8,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad8(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file9,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad9(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file10,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad10(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file11,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad11(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file12,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad12(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file13,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad13(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file14,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad14(i)=b
                endfor
                
                free_lun,lun
                openr,lun,file15,/get_lun
                for i = 0, 27-1 do begin
                  readf,lun,a,b
                  rad15(i)=b
                endfor
                free_lun,lun
                
                wn = (1/wn)*10000
                wn = reverse(wn)

                ; ver1
                ; band=where(wn gt 1.85 and wn lt 2.10)
                
                ; ver2
                ; band=where(wn gt 1.94 and wn lt 2.09)

                ; ver3
                ; band=where(wn gt 1.94 and wn lt 1.99)

                ; ver4
                band = where(wn gt 1.93 and wn lt 2.04)
                
                rad1 = reverse(rad1)
                rad2 = reverse(rad2)
                rad3 = reverse(rad3)
                rad4 = reverse(rad4)
                rad5 = reverse(rad5)
                rad6 = reverse(rad6)
                rad7 = reverse(rad7)
                rad8 = reverse(rad8)
                rad9 = reverse(rad9)
                rad10 = reverse(rad10)
                rad11 = reverse(rad11)
                rad12 = reverse(rad12)
                rad13 = reverse(rad13)
                rad14 = reverse(rad14)
                rad15 = reverse(rad15)
                
                ;Cal equivalent width (absorption depth)
                x = [wn(0), wn(1), wn(2), wn(23), wn(24), wn(25)]
                y1 = [rad1(0), rad1(1), rad1(2), rad1(23), rad1(24), rad1(25)]
                y2 = [rad2(0), rad2(1), rad2(2), rad2(23), rad2(24), rad2(25)]
                y3 = [rad3(0), rad3(1), rad3(2), rad3(23), rad3(24), rad3(25)]
                y4 = [rad4(0), rad4(1), rad4(2), rad4(23), rad4(24), rad4(25)]
                y5 = [rad5(0), rad5(1), rad5(2), rad5(23), rad5(24), rad5(25)]
                y6 = [rad6(0), rad6(1), rad6(2), rad6(23), rad6(24), rad6(25)]
                y7 = [rad7(0), rad7(1), rad7(2), rad7(23), rad7(24), rad7(25)]
                y8 = [rad8(0), rad8(1), rad8(2), rad8(23), rad8(24), rad8(25)]
                y9 = [rad9(0), rad9(1), rad9(2), rad9(23), rad9(24), rad9(25)]
                y10 = [rad10(0), rad10(1), rad10(2), rad10(23), rad10(24), rad10(25)]
                y11 = [rad11(0), rad11(1), rad11(2), rad11(23), rad11(24), rad11(25)]
                y12 = [rad12(0), rad12(1), rad12(2), rad12(23), rad12(24), rad12(25)]
                y13 = [rad13(0), rad13(1), rad13(2), rad13(23), rad13(24), rad13(25)]
                y14 = [rad14(0), rad14(1), rad14(2), rad14(23), rad14(24), rad14(25)]
                y15 = [rad15(0), rad15(1), rad15(2), rad15(23), rad15(24), rad15(25)]
                
                coef1 = linfit(X,Y1)
                coef2 = linfit(X,Y2)
                coef3 = linfit(X,Y3)
                coef4 = linfit(X,Y4)
                coef5 = linfit(X,Y5)
                coef6 = linfit(X,Y6)
                coef7 = linfit(X,Y7)
                coef8 = linfit(X,Y8)
                coef9 = linfit(X,Y9)
                coef10 = linfit(X,Y10)
                coef11 = linfit(X,Y11)
                coef12 = linfit(X,Y12)
                coef13 = linfit(X,Y13)
                coef14 = linfit(X,Y14)
                coef15 = linfit(X,Y15)
                
                cont1 = coef1(0) + coef1(1)*wn
                width1 = 1.0 - rad1/cont1
                width1[16] = !VALUES.F_NAN
                                
                cont2 = coef2(0) + coef2(1)*wn
                width2 = 1.0 - rad2/cont2
                width2[16] = !VALUES.F_NAN
                
                cont3 = coef3(0) + coef3(1)*wn
                width3 = 1.0 - rad3/cont3
                width3[16] = !VALUES.F_NAN
                
                cont4 = coef4(0) + coef4(1)*wn
                width4 = 1.0 - rad4/cont4
                width4[16] = !VALUES.F_NAN
                
                cont5 = coef5(0) + coef5(1)*wn
                width5 = 1.0 - rad5/cont5
                width5[16] = !VALUES.F_NAN
                
                cont6 = coef6(0) + coef6(1)*wn
                width6 = 1.0 - rad6/cont6
                width6[16] = !VALUES.F_NAN
                
                cont7 = coef7(0) + coef7(1)*wn
                width7 = 1.0 - rad7/cont7
                width7[16] = !VALUES.F_NAN
                
                cont8 = coef8(0) + coef8(1)*wn
                width8 = 1.0 - rad8/cont8
                width8[16] = !VALUES.F_NAN
                
                cont9 = coef9(0) + coef9(1)*wn
                width9 = 1.0 - rad9/cont9
                width9[16] = !VALUES.F_NAN
                
                cont10 = coef10(0) + coef10(1)*wn
                width10 = 1.0 - rad10/cont10
                width10[16] = !VALUES.F_NAN
                
                cont11 = coef11(0) + coef11(1)*wn
                width11 = 1.0 - rad11/cont11
                width11[16] = !VALUES.F_NAN
                
                cont12 = coef12(0) + coef12(1)*wn
                width12 = 1.0 - rad12/cont12
                width12[16] = !VALUES.F_NAN
                
                cont13 = coef13(0) + coef13(1)*wn
                width13 = 1.0 - rad13/cont13
                width13[16] = !VALUES.F_NAN
                
                cont14 = coef14(0) + coef14(1)*wn
                width14 = 1.0 - rad14/cont14
                width14[16] = !VALUES.F_NAN
                
                cont15 = coef15(0) + coef15(1)*wn
                width15 = 1.0 - rad15/cont15
                width15[16] = !VALUES.F_NAN


;                stop
;                window,0
;                plot, wn, rad7 / (coef7(0) + coef7(1)*wn), xs=1, ys=1
;                stop

                Table_Equivalent_pressure1(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width1[band],/nan)
                Table_Equivalent_pressure2(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width2[band],/nan)
                Table_Equivalent_pressure3(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width3[band],/nan)
                Table_Equivalent_pressure4(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width4[band],/nan)
                Table_Equivalent_pressure5(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width5[band],/nan)
                Table_Equivalent_pressure6(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width6[band],/nan)
                Table_Equivalent_pressure7(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width7[band],/nan)
                Table_Equivalent_pressure8(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width8[band],/nan)
                Table_Equivalent_pressure9(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width9[band],/nan)
                Table_Equivalent_pressure10(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width10[band],/nan)
                Table_Equivalent_pressure11(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width11[band],/nan)
                Table_Equivalent_pressure12(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width12[band],/nan)
                Table_Equivalent_pressure13(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width13[band],/nan)
                Table_Equivalent_pressure14(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width14[band],/nan)
                Table_Equivalent_pressure15(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width15[band],/nan)

                if y1(0) ge 100 then stop
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor

save,Table_Equivalent_pressure1,$
     Table_Equivalent_pressure2,$
     Table_Equivalent_pressure3,$
     Table_Equivalent_pressure4,$
     Table_Equivalent_pressure5,$
     Table_Equivalent_pressure6,$
     Table_Equivalent_pressure7,$
     Table_Equivalent_pressure8,$
     Table_Equivalent_pressure9,$
     Table_Equivalent_pressure10,$
     Table_Equivalent_pressure11,$
     Table_Equivalent_pressure12,$
     Table_Equivalent_pressure13,$
     Table_Equivalent_pressure14,$
     Table_Equivalent_pressure15,$
     
     filename='/work1/LUT/SP/table/absorption/Table_SP_calc_ver4_add_albedo.sav'
stop
END