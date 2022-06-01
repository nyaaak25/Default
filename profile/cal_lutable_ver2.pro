;------------------------------------------------------------------------------------------------------------------------
Pro Cal_LUtable
; +++
; to read outputs from ARS and calculate equivalent width (=absorption depth) for each spectrum
; +++
;------------------------------------------------------------------------------------------------------------------------
path = '/work1/LUT/SP/table/output/'

Table_Equivalent_pressure1=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure2=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure3=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure4=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure5=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure6=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure7=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure8=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure9=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure10=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure11=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure12=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure13=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure14=dblarr(5,3,6,3,5,6,3,6)
Table_Equivalent_pressure15=dblarr(5,3,6,3,5,6,3,6)

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

x = dblarr(2)
y1 = dblarr(2)
y2 = dblarr(2)
y3 = dblarr(2)
y4 = dblarr(2)
y5 = dblarr(2)
y6 = dblarr(2)
y7 = dblarr(2)
y8 = dblarr(2)
y9 = dblarr(2)
y10 = dblarr(2)
y11 = dblarr(2)
y12 = dblarr(2)
y13 = dblarr(2)
y14 = dblarr(2)
y15 = dblarr(2)

for IT1 = 1, 5 do begin
  for IT2 = 1, 3 do begin
    for ISZA = 1, 6 do begin
      for IEA = 1, 3 do begin
        for IPA = 1, 5 do begin
          for ID = 1, 6 do begin
            for IWI = 1, 3 do begin
              for IAB = 1, 6 do begin
                
                ;file names
                file1 = path + 'SP1' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file2 = path + 'SP2' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file3 = path + 'SP3' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                
                file4 = path + 'SP4' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file5 = path + 'SP5' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file6 = path + 'SP6' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'    
                
                file7 = path + 'SP7' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file8 = path + 'SP8' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file9 = path + 'SP9' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file10 = path + 'SP10' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'

                file11 = path + 'SP11' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file12 = path + 'SP12' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                     
                file13 = path + 'SP13' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat' 
                                     
                file14 = path + 'SP14' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                     + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                     + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                     + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                     + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                     + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                     + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                     + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'                                                                                                            

                file15 = path + 'SP15' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
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
                

                ;Cal equivalent width (absorption depth)
                x = [wn(0), wn(3), wn(5), wn(23), wn(24), wn(25)]
                y1 = [rad1(0), rad1(3), rad1(5), rad1(23), rad1(24), rad1(25)]
                y2 = [rad2(0), rad2(3), rad2(5), rad2(23), rad2(24), rad2(25)]
                y3 = [rad3(0), rad3(3), rad3(5), rad3(23), rad3(24), rad3(25)]
                y4 = [rad4(0), rad4(3), rad4(5), rad4(23), rad4(24), rad4(25)]
                y5 = [rad5(0), rad5(3), rad5(5), rad5(23), rad5(24), rad5(25)]
                y6 = [rad6(0), rad6(3), rad6(5), rad6(23), rad6(24), rad6(25)]
                y7 = [rad7(0), rad7(3), rad7(5), rad7(23), rad7(24), rad7(25)]
                y8 = [rad8(0), rad8(3), rad8(5), rad8(23), rad8(24), rad8(25)]
                y9 = [rad9(0), rad9(3), rad9(5), rad9(23), rad9(24), rad9(25)]
                y10 = [rad10(0), rad10(3), rad10(5), rad10(23), rad10(24), rad10(25)]
                y11 = [rad11(0), rad11(3), rad11(5), rad11(23), rad11(24), rad11(25)]
                y12 = [rad12(0), rad12(3), rad12(5), rad12(23), rad12(24), rad12(25)]
                y13 = [rad13(0), rad13(3), rad13(5), rad13(23), rad13(24), rad13(25)]
                y14 = [rad14(0), rad14(3), rad14(5), rad14(23), rad14(24), rad14(25)]
                y15 = [rad15(0), rad15(3), rad15(5), rad15(23), rad15(24), rad15(25)]
                
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
                
                cont1_1 = coef1(0) + coef1(1)*wn(8)
                cont1_2 = coef1(0) + coef1(1)*wn(12)
                cont1_3 = coef1(0) + coef1(1)*wn(16)
                width1_1 = 1.0 - rad1(8)/cont1_1
                width1_2 = 1.0 - rad1(12)/cont1_2
                width1_3 = 1.0 - rad1(16)/cont1_3
                
                cont2_1 = coef2(0) + coef2(1)*wn(8)
                cont2_2 = coef2(0) + coef2(1)*wn(12)
                cont2_3 = coef2(0) + coef2(1)*wn(16)
                width2_1 = 1.0 - rad2(8)/cont2_1
                width2_2 = 1.0 - rad2(12)/cont2_2
                width2_3 = 1.0 - rad2(16)/cont2_3
                
                cont3_1 = coef3(0) + coef3(1)*wn(8)
                cont3_2 = coef3(0) + coef3(1)*wn(12)
                cont3_3 = coef3(0) + coef3(1)*wn(16)
                width3_1 = 1.0 - rad3(8)/cont3_1
                width3_2 = 1.0 - rad3(12)/cont3_2
                width3_3 = 1.0 - rad3(16)/cont3_3
                
                cont4_1 = coef4(0) + coef4(1)*wn(8)
                cont4_2 = coef4(0) + coef4(1)*wn(12)
                cont4_3 = coef4(0) + coef4(1)*wn(16)
                width4_1 = 1.0 - rad4(8)/cont4_1
                width4_2 = 1.0 - rad4(12)/cont4_2
                width4_3 = 1.0 - rad4(16)/cont4_3
                
                cont5_1 = coef5(0) + coef5(1)*wn(8)
                cont5_2 = coef5(0) + coef5(1)*wn(12)
                cont5_3 = coef5(0) + coef5(1)*wn(16)
                width5_1 = 1.0 - rad5(8)/cont5_1
                width5_2 = 1.0 - rad5(12)/cont5_2
                width5_3 = 1.0 - rad5(16)/cont5_3
                
                cont6_1 = coef6(0) + coef6(1)*wn(8)
                cont6_2 = coef6(0) + coef6(1)*wn(12)
                cont6_3 = coef6(0) + coef6(1)*wn(16)
                width6_1 = 1.0 - rad6(8)/cont6_1
                width6_2 = 1.0 - rad6(12)/cont6_2
                width6_3 = 1.0 - rad6(16)/cont6_3
                
                cont7_1 = coef7(0) + coef7(1)*wn(8)
                cont7_2 = coef7(0) + coef7(1)*wn(12)
                cont7_3 = coef7(0) + coef7(1)*wn(16)
                width7_1 = 1.0 - rad7(8)/cont7_1
                width7_2 = 1.0 - rad7(12)/cont7_2
                width7_3 = 1.0 - rad7(16)/cont7_3
                
                cont8_1 = coef8(0) + coef8(1)*wn(8)
                cont8_2 = coef8(0) + coef8(1)*wn(12)
                cont8_3 = coef8(0) + coef8(1)*wn(16)
                width8_1 = 1.0 - rad8(8)/cont8_1
                width8_2 = 1.0 - rad8(12)/cont8_2
                width8_3 = 1.0 - rad8(16)/cont8_3
                
                cont9_1 = coef9(0) + coef9(1)*wn(8)
                cont9_2 = coef9(0) + coef9(1)*wn(12)
                cont9_3 = coef9(0) + coef9(1)*wn(16)
                width9_1 = 1.0 - rad9(8)/cont9_1
                width9_2 = 1.0 - rad9(12)/cont9_2
                width9_3 = 1.0 - rad9(16)/cont9_3
                
                cont10_1 = coef10(0) + coef10(1)*wn(8)
                cont10_2 = coef10(0) + coef10(1)*wn(12)
                cont10_3 = coef10(0) + coef10(1)*wn(16)
                width10_1 = 1.0 - rad10(8)/cont10_1
                width10_2 = 1.0 - rad10(12)/cont10_2
                width10_3 = 1.0 - rad10(16)/cont10_3
                
                cont11_1 = coef11(0) + coef11(1)*wn(8)
                cont11_2 = coef11(0) + coef11(1)*wn(12)
                cont11_3 = coef11(0) + coef11(1)*wn(16)
                width11_1 = 1.0 - rad11(8)/cont11_1
                width11_2 = 1.0 - rad11(12)/cont11_2
                width11_3 = 1.0 - rad11(16)/cont11_3
                
                cont12_1 = coef12(0) + coef12(1)*wn(8)
                cont12_2 = coef12(0) + coef12(1)*wn(12)
                cont12_3 = coef12(0) + coef12(1)*wn(16)
                width12_1 = 1.0 - rad12(8)/cont12_1
                width12_2 = 1.0 - rad12(12)/cont12_2
                width12_3 = 1.0 - rad12(16)/cont12_3
                
                cont13_1 = coef13(0) + coef13(1)*wn(8)
                cont13_2 = coef13(0) + coef13(1)*wn(12)
                cont13_3 = coef13(0) + coef13(1)*wn(16)
                width13_1 = 1.0 - rad13(8)/cont13_1
                width13_2 = 1.0 - rad13(12)/cont13_2
                width13_3 = 1.0 - rad13(16)/cont13_3
                
                cont14_1 = coef14(0) + coef14(1)*wn(8)
                cont14_2 = coef14(0) + coef14(1)*wn(12)
                cont14_3 = coef14(0) + coef14(1)*wn(16)
                width14_1 = 1.0 - rad14(8)/cont14_1
                width14_2 = 1.0 - rad14(12)/cont14_2
                width14_3 = 1.0 - rad14(16)/cont14_3
                
                cont15_1 = coef15(0) + coef15(1)*wn(8)
                cont15_2 = coef15(0) + coef15(1)*wn(12)
                cont15_3 = coef15(0) + coef15(1)*wn(16)
                width15_1 = 1.0 - rad15(8)/cont15_1
                width15_2 = 1.0 - rad15(12)/cont15_2
                width15_3 = 1.0 - rad15(16)/cont15_3

;                stop
;                window,0
;                plot, wn, rad7 / (coef7(0) + coef7(1)*wn), xs=1, ys=1
;                stop

                Table_Equivalent_pressure1(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width1_1 + width1_2 + width1_3
                Table_Equivalent_pressure2(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width2_1 + width2_2 + width2_3
                Table_Equivalent_pressure3(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width3_1 + width3_2 + width3_3
                Table_Equivalent_pressure4(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width4_1 + width4_2 + width4_3
                Table_Equivalent_pressure5(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width5_1 + width5_2 + width5_3
                Table_Equivalent_pressure6(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width6_1 + width6_2 + width6_3
                Table_Equivalent_pressure7(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width7_1 + width7_2 + width7_3
                Table_Equivalent_pressure8(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width8_1 + width8_2 + width8_3
                Table_Equivalent_pressure9(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width9_1 + width9_2 + width9_3
                Table_Equivalent_pressure10(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width10_1 + width10_2 + width10_3
                Table_Equivalent_pressure11(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width11_1 + width11_2 + width11_3
                Table_Equivalent_pressure12(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width12_1 + width12_2 + width12_3
                Table_Equivalent_pressure13(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width13_1 + width13_2 + width13_3
                Table_Equivalent_pressure14(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width14_1 + width14_2 + width14_3
                Table_Equivalent_pressure15(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = width15_1 + width15_2 + width15_3

                if y0(0) ge 100 then stop
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
     
     filename='/work1/LUT/SP/table/absorption/Table_SP_Trans.sav'
stop
END