;------------------------------------------------------------------------------------------------------------------------
Pro Cal_LUtable_ver2_nyonn
; +++
; to read outputs from ARS and calculate equivalent width (=absorption depth) for each spectrum
; create by Shohei Aoki
; 
; edit by Akira Kazama
; 
; cal_lutable_obsspec　 ::2022.6.1 Wed 11:43:00
; OMEGAで観測されたスペクトルの吸収深さ量を計算するプログラム
; 
; +++
;------------------------------------------------------------------------------------------------------------------------
path = '/Users/nyonn/IDLWorkspace/Default/savfile/'


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
                
                cont1 = coef1(0) + coef1(1)*wn
                width1 = 1.0 - rad1/cont1
                                
                cont2 = coef2(0) + coef2(1)*wn
                width2 = 1.0 - rad2/cont2
                
                cont3 = coef3(0) + coef3(1)*wn
                width3 = 1.0 - rad3/cont3
                
                cont4 = coef4(0) + coef4(1)*wn
                width4 = 1.0 - rad4/cont4
                
                cont5 = coef5(0) + coef5(1)*wn
                width5 = 1.0 - rad5/cont5
                
                cont6 = coef6(0) + coef6(1)*wn
                width6 = 1.0 - rad6/cont6
                
                cont7 = coef7(0) + coef7(1)*wn
                width7 = 1.0 - rad7/cont7
                
                cont8 = coef8(0) + coef8(1)*wn
                width8 = 1.0 - rad8/cont8
                
                cont9 = coef9(0) + coef9(1)*wn
                width9 = 1.0 - rad9/cont9
                
                cont10 = coef10(0) + coef10(1)*wn
                width10 = 1.0 - rad10/cont10
                
                cont11 = coef11(0) + coef11(1)*wn
                width11 = 1.0 - rad11/cont11
                
                cont12 = coef12(0) + coef12(1)*wn
                width12 = 1.0 - rad12/cont12
                
                cont13 = coef13(0) + coef13(1)*wn
                width13 = 1.0 - rad13/cont13
                
                cont14 = coef14(0) + coef14(1)*wn
                width14 = 1.0 - rad14/cont14
                
                cont15 = coef15(0) + coef15(1)*wn
                width15 = 1.0 - rad15/cont15


;                stop
;                window,0
;                plot, wn, rad7 / (coef7(0) + coef7(1)*wn), xs=1, ys=1
;                stop

                Table_Equivalent_pressure1(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = total(width1[7:18])

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
     
     filename='/work1/LUT/SP/table/absorption/Table_SP_Trans_calc.sav'
stop
END