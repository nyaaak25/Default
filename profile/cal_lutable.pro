;------------------------------------------------------------------------------------------------------------------------
Pro Cal_LUtable
; +++
; to read outputs from ARS and calculate equivalent width (=absorption depth) for each spectrum
; +++
;------------------------------------------------------------------------------------------------------------------------
path = '/Users/Shohei/tmp/OMEGA_Nadir/Output_ver2/'

Table_Equivalent_width0=dblarr(4,3,5,3,4,5,3,5)
Table_Equivalent_width1=dblarr(4,3,5,3,4,5,3,5)
Table_Equivalent_width2=dblarr(4,3,5,3,4,5,3,5)
Table_Equivalent_width3=dblarr(4,3,5,3,4,5,3,5)
Table_Equivalent_width4=dblarr(4,3,5,3,4,5,3,5)
Table_Equivalent_width5=dblarr(4,3,5,3,4,5,3,5)
Table_Equivalent_width6=dblarr(4,3,5,3,4,5,3,5)
Table_Equivalent_width7=dblarr(4,3,5,3,4,5,3,5)
Table_Equivalent_width10=dblarr(4,3,5,3,4,5,3,5)
Table_Equivalent_width20=dblarr(4,3,5,3,4,5,3,5)

wn=dblarr(3)
rad0=dblarr(3)
rad1=dblarr(3)
rad2=dblarr(3)
rad3=dblarr(3)
rad4=dblarr(3)
rad5=dblarr(3)
rad6=dblarr(3)
rad7=dblarr(3)

x = dblarr(2)
y0 = dblarr(2)
y1 = dblarr(2)
y2 = dblarr(2)
y3 = dblarr(2)
y4 = dblarr(2)
y5 = dblarr(2)
y6 = dblarr(2)
y7 = dblarr(2)

for IT1 = 1, 4 do begin
  for IT2 = 1, 3 do begin
    for ISZA = 1, 5 do begin
      for IEA = 1, 3 do begin
        for IPA = 1, 4 do begin
          for ID = 1, 5 do begin
            for IWI = 1, 3 do begin
              for IAB = 1, 5 do begin
                
                ;file names
                file0 = path+'Water0'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'
                file1 = path+'Water1'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'
                file2 = path+'Water2'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'
                file3 = path+'Water3'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'
                file4 = path+'Water4'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'
                file5 = path+'Water5'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'
                file6 = path+'Water6'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'
                file7 = path+'Water7'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'
                 file10 = path+'Water10'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'
                 file20 = path+'Water20'+'_TA'+STRCOMPRESS(IT1,/REMOVE_AL)+'_TB'+STRCOMPRESS(IT2,/REMOVE_AL)+'_SZA'+STRCOMPRESS(ISZA,/REMOVE_AL)+'_EA'+STRCOMPRESS(IEA,/REMOVE_AL) $
                                +'_PA'+STRCOMPRESS(IPA,/REMOVE_AL)+'_Dust' + STRCOMPRESS(ID,/REMOVE_AL)+'_WaterI' + STRCOMPRESS(IWI,/REMOVE_AL)+'_SurfaceA'+STRCOMPRESS(IAB,/REMOVE_AL)+'.dat'

                ;read files
                openr,lun,file10,/get_lun
                for i = 0, 3-1 do begin
                  readf,lun,a,b
                  wn(i)=a
                  rad0(i)=b
                endfor
                free_lun,lun
                openr,lun,file20,/get_lun
                for i = 0, 3-1 do begin
                  readf,lun,a,b
                  rad1(i)=b
                endfor
                free_lun,lun
                openr,lun,file2,/get_lun
                for i = 0, 3-1 do begin
                  readf,lun,a,b
                  rad2(i)=b
                endfor
                free_lun,lun
                openr,lun,file3,/get_lun
                for i = 0, 3-1 do begin
                  readf,lun,a,b
                  rad3(i)=b
                endfor
                free_lun,lun
                openr,lun,file4,/get_lun
                for i = 0, 3-1 do begin
                  readf,lun,a,b
                  rad4(i)=b
                endfor
                free_lun,lun
                openr,lun,file5,/get_lun
                for i = 0, 3-1 do begin
                  readf,lun,a,b
                  rad5(i)=b
                endfor
                free_lun,lun
                openr,lun,file6,/get_lun
                for i = 0, 3-1 do begin
                  readf,lun,a,b
                  rad6(i)=b
                endfor
                free_lun,lun
                openr,lun,file7,/get_lun
                for i = 0, 3-1 do begin
                  readf,lun,a,b
                  rad7(i)=b
                endfor
                free_lun,lun

                ;Cal equivalent width (absorption depth)
                x = [wn(0), wn(2)]
                y0 = [rad0(0), rad0(2)]
                y1 = [rad1(0), rad1(2)]
                y2 = [rad2(0), rad2(2)]
                y3 = [rad3(0), rad3(2)]
                y4 = [rad4(0), rad4(2)]
                y5 = [rad5(0), rad5(2)]
                y6 = [rad6(0), rad6(2)]
                y7 = [rad7(0), rad7(2)]
                coef0 = linfit(X,Y0)
                coef1 = linfit(X,Y1)
                coef2 = linfit(X,Y2)
                coef3 = linfit(X,Y3)
                coef4 = linfit(X,Y4)
                coef5 = linfit(X,Y5)
                coef6 = linfit(X,Y6)
                coef7 = linfit(X,Y7)
                cont0 = coef0(0) + coef0(1)*wn(1)
                cont1 = coef1(0) + coef1(1)*wn(1)
                cont2 = coef2(0) + coef2(1)*wn(1)
                cont3 = coef3(0) + coef3(1)*wn(1)
                cont4 = coef4(0) + coef4(1)*wn(1)
                cont5 = coef5(0) + coef5(1)*wn(1)
                cont6 = coef6(0) + coef6(1)*wn(1)
                cont7 = coef7(0) + coef7(1)*wn(1)

;                stop
;                window,0
;                plot, wn, rad7 / (coef7(0) + coef7(1)*wn), xs=1, ys=1
;                stop

                Table_Equivalent_width0(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad0(1)/cont0
                Table_Equivalent_width1(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad1(1)/cont1
                Table_Equivalent_width2(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad2(1)/cont2
                Table_Equivalent_width3(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad3(1)/cont3
                Table_Equivalent_width4(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad4(1)/cont4
                Table_Equivalent_width5(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad5(1)/cont5
                Table_Equivalent_width6(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad6(1)/cont6
                Table_Equivalent_width7(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad7(1)/cont7
                Table_Equivalent_width10(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad0(1)/cont0
                Table_Equivalent_width20(IT1-1,IT2-1,ISZA-1,IEA-1,IPA-1,ID-1,IWI-1,IAB-1) = 1.0 - rad1(1)/cont1

                if y0(0) ge 100 then stop
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor
  endfor
endfor

save,Table_Equivalent_width10,$
     Table_Equivalent_width1,$
     Table_Equivalent_width2,$
     Table_Equivalent_width3,$
     Table_Equivalent_width4,$
     Table_Equivalent_width5,$
     Table_Equivalent_width6,$
     Table_Equivalent_width20,$
     filename='/Users/Shohei/tmp/OMEGA_Nadir/Table_water_Trans_ver2.sav'
stop
END