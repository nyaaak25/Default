;------------------------------------------------------------------------------------------------------------------------
Pro Cal_LUtable_water_v1
; +++
; to read outputs from ARS and calculate equivalent width (=absorption depth) for each spectrum
; +++
;------------------------------------------------------------------------------------------------------------------------
path_LUT_Output1 = '/work1/LUT/H2O/table/output/'
path_LUT_Output2 = '/work1/LUT/H2O/table/output2/'


Table_Equivalent_width1 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width2 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width3 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width4 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width5 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width6 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width7 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width8 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width9 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width10 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width11 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width12 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width13 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width14 = dblarr(15,5,3,6,5,5,6,6)
Table_Equivalent_width15 = dblarr(15,5,3,6,5,5,6,6)

x = dblarr(11)
y = dblarr(11)

for ISP = 1, 15 do begin
  for IT1 = 2, 5 do begin
    print, 'ISP:', ISP, ' IT1:', IT1
    for IT2 = 1, 3 do begin
;      print, 'IT2:', IT2
      for ISZA = 1, 6 do begin
        for IEA = 1, 5 do begin
          for IPA = 1, 5 do begin
            for ID = 1, 6 do begin
              for IAB = 1, 6 do begin
                
               output_name =  '_SP' + STRCOMPRESS(fix(ISP),/REMOVE_AL) $
                      + '_TA' + STRCOMPRESS(IT1,/REMOVE_AL) $
                      + '_TB' + STRCOMPRESS(IT2,/REMOVE_AL) $
                      + '_SZA' + STRCOMPRESS(ISZA,/REMOVE_AL) $
                      + '_EA' + STRCOMPRESS(IEA,/REMOVE_AL) $
                      + '_PA' + STRCOMPRESS(IPA,/REMOVE_AL) $
                      + '_Dust' + STRCOMPRESS(ID,/REMOVE_AL) $
                      + '_WaterI' + STRCOMPRESS(1,/REMOVE_AL) $
                      + '_SurfaceA' + STRCOMPRESS(IAB,/REMOVE_AL) + '_rad.dat'

                ;--------------
                file = path_LUT_Output1 + 'Water1' + output_name          
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width1(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip1
                endif

                file = path_LUT_Output2 + 'Water1' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width1(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip1:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water2' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width2(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip2
                endif

                file = path_LUT_Output2 + 'Water2' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width1(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip2:
                ;--------------

                ;--------------
                file = path_LUT_Output1 + 'Water3' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width3(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip3
                endif

                file = path_LUT_Output2 + 'Water3' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width3(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip3:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water4' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width4(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip4
                endif

                file = path_LUT_Output2 + 'Water4' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width4(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip4:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water5' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width5(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip5
                endif

                file = path_LUT_Output2 + 'Water5' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width5(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip5:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water6' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width6(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip6
                endif

                file = path_LUT_Output2 + 'Water6' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width6(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip6:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water7' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width7(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip7
                endif

                file = path_LUT_Output2 + 'Water7' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width7(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip7:
                ;--------------                
                ;--------------
                file = path_LUT_Output1 + 'Water8' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width8(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip8
                endif

                file = path_LUT_Output2 + 'Water8' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width8(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip8:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water9' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width9(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip9
                endif

                file = path_LUT_Output2 + 'Water9' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width9(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip9:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water10' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width10(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip10
                endif

                file = path_LUT_Output2 + 'Water10' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width10(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip10:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water11' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width11(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip11
                endif

                file = path_LUT_Output2 + 'Water11' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width11(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip11:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water12' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width12(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip12
                endif

                file = path_LUT_Output2 + 'Water12' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width12(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip12:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water13' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width13(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip13
                endif

                file = path_LUT_Output2 + 'Water13' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width13(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip13:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water14' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width14(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                  goto, skip14
                endif

                file = path_LUT_Output2 + 'Water14' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width14(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip14:
                ;--------------
                ;--------------
                file = path_LUT_Output1 + 'Water15' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width15(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
;                  stop
                  goto, skip15
                endif

                file = path_LUT_Output2 + 'Water15' + output_name
                if file_test(file) eq 1 then begin
                  openr,lun,file,/get_lun
                  for i = 0, 11-1 do begin
                    readf,lun,a,b
                    x(i)=a
                    y(i)=b
                  endfor
                  free_lun,lun
                  x_local = reform([x(10),x(9),x(7),x(6),x(5),x(4),x(3),x(2),x(1),x(0)])
                  y_local = reform([y(10),y(9),y(7),y(6),y(5),y(4),y(3),y(2),y(1),y(0)])
                  x_cont = reform([x_local(0), x_local(1), x_local(2), x_local(9)])
                  y_cont = reform([y_local(0), y_local(1), y_local(2), y_local(9)])
                  coef = linfit(X_cont, Y_cont)
                  cont = coef(0) + coef(1)*x_local
                  Table_Equivalent_width15(ISP-1,IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IAB-1) = total(1.0 - Y_local/cont)
                endif
                skip15:
                ;--------------
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor

save,Table_Equivalent_width1,$
     Table_Equivalent_width2,$
     Table_Equivalent_width3,$
     Table_Equivalent_width4,$
     Table_Equivalent_width5,$
     Table_Equivalent_width6,$
     Table_Equivalent_width7,$
     Table_Equivalent_width8,$
     Table_Equivalent_width9,$
     Table_Equivalent_width10,$
     Table_Equivalent_width11,$
     Table_Equivalent_width12,$
     Table_Equivalent_width13,$
     Table_Equivalent_width14,$
     Table_Equivalent_width15,$
     filename='/work1/LUT/H2O/table/Table_water_EW_v1.sav'
stop
END