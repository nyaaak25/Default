;------------------------------------------------------------------------------------------------------------------------
Pro Cal_LUtable_fitting
; +++
; to read outputs from ARS and calculate equivalent width (=absorption depth) for each spectrum
; create by Shohei Aoki
; 
; edit by Akira Kazama
; 
; cal_lutable_fitting　 ::2022.9.1 Thu 16:21:00
; ver1: 波長1つ1つに対して吸収の深さを計算させて、sav fileを作成するプログラム
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

wn=dblarr(29)
rad1=dblarr(29)
rad2=dblarr(29)
rad3=dblarr(29)
rad4=dblarr(29)
rad5=dblarr(29)
rad6=dblarr(29)
rad7=dblarr(29)
rad8=dblarr(29)
rad9=dblarr(29)
rad10=dblarr(29)
rad11=dblarr(29)
rad12=dblarr(29)
rad13=dblarr(29)
rad14=dblarr(29)
rad15=dblarr(29)

for w = 0, 29-1 do begin
  for IT1 = 1, 5 do begin
    for IT2 = 1, 3 do begin
      for ISZA = 1, 6 do begin
        for IEA = 1, 5 do begin
          for IPA = 1, 5 do begin
            for ID = 1, 6 do begin;6 do begin
              for IWI = 1, 1 do begin;3 do begin
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
;  
                  file15 = path15 + 'SP15' + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                       + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                       + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                       + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                       + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                       + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                       + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                       + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad_test.dat'
  
  
                  ;read files
                  openr,lun,file1,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    wn(i)=a
                    rad1(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file2,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad2(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file3,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad3(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file4,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad4(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file5,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad5(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file6,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad6(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file7,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad7(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file8,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad8(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file9,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad9(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file10,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad10(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file11,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad11(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file12,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad12(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file13,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad13(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file14,/get_lun
                  for i = 0, 29-1 do begin
                    readf,lun,a,b
                    rad14(i)=b
                  endfor
                  
                  free_lun,lun
                  openr,lun,file15,/get_lun
                  for i = 0, 27-1 do begin
                    readf,lun,a,b
                    wn(i)=a
                    rad15(i)=b
                  endfor
                  free_lun,lun

                  wn = reverse(wn)
                  wav = 1/wn
                  wn = (1/wn)*10000
                  
                  rad1 = reverse(rad1)
                  rad1 = (rad1/wav^2)*1e-7  ; 波数から波長換算 (erg = 10-7W, cm^2 = 10-4m^2, μm = 10^4cm)

                  rad2 = reverse(rad2)
                  rad2 = (rad2/wav^2)*1e-7

                  rad3 = reverse(rad3)
                  rad3 = (rad3/wav^2)*1e-7

                  rad4 = reverse(rad4)
                  rad4 = (rad4/wav^2)*1e-7

                  rad5 = reverse(rad5)
                  rad5 = (rad5/wav^2)*1e-7

                  rad6 = reverse(rad6)
                  rad6 = (rad6/wav^2)*1e-7

                  rad7 = reverse(rad7)
                  rad7 = (rad7/wav^2)*1e-7

                  rad8 = reverse(rad8)
                  rad8 = (rad8/wav^2)*1e-7

                  rad9 = reverse(rad9)
                  rad9 = (rad9/wav^2)*1e-7

                  rad10 = reverse(rad10)
                  rad10 = (rad10/wav^2)*1e-7

                  rad11 = reverse(rad11)
                  rad11 = (rad11/wav^2)*1e-7

                  rad12 = reverse(rad12)
                  rad12 = (rad12/wav^2)*1e-7

                  rad13 = reverse(rad13)
                  rad13 = (rad13/wav^2)*1e-7

                  rad14 = reverse(rad14)
                  rad14 = (rad14/wav^2)*1e-7

                  rad15 = reverse(rad15)
                  rad15 = (rad15/wav^2)*1e-7

                  Table_Equivalent_pressure1(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad1[w]
                  Table_Equivalent_pressure2(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad2[w]
                  Table_Equivalent_pressure3(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad3[w]
                  Table_Equivalent_pressure4(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad4[w]
                  Table_Equivalent_pressure5(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad5[w]
                  Table_Equivalent_pressure6(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad6[w]
                  Table_Equivalent_pressure7(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad7[w]
                  Table_Equivalent_pressure8(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad8[w]
                  Table_Equivalent_pressure9(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad9[w]
                  Table_Equivalent_pressure10(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad10[w]
                  Table_Equivalent_pressure11(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad11[w]
                  Table_Equivalent_pressure12(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad12[w]
                  Table_Equivalent_pressure13(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad13[w]
                  Table_Equivalent_pressure14(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad14[w]
                  Table_Equivalent_pressure15(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = rad15[w]
                  
                  ;check
;                  if IAB eq 1 then ref = rad15[2]
;                  if IAB ne 1 then print, rad15[2]/ref, (IAB-1)*0.1/0.05d;, '  ', file15
;                  if IAB ne 1 then spawn, 'ls -l '+file15
                  
                endfor
;                stop
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
       
       filename='/work1/LUT/SP/table/LUT_fitting/Table_calc_wave_new_update_rad'+ STRCOMPRESS(fix(w),/REMOVE_AL) + '.sav'
       
endfor
END