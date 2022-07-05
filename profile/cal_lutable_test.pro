;------------------------------------------------------------------------------------------------------------------------
Pro Cal_lutable_test
; +++
; to read outputs from ARS and calculate equivalent width (=absorption depth) for each spectrum
; +++
;------------------------------------------------------------------------------------------------------------------------
file1 = '/Users/nyonn/IDLWorkspace/Default/LUT/SP1_TA1_TB1_SZA1_EA1_PA1_Dust1_WaterI1_SurfaceA1_rad.dat'
file2 = '/Users/nyonn/IDLWorkspace/Default/LUT/SP1_TA1_TB1_SZA1_EA3_PA1_Dust5_WaterI2_SurfaceA6_rad.dat'

Table_Equivalent_pressure1=dblarr(5,3,6,3,5,6,3,6)

wn=dblarr(27)
rad0=dblarr(27)
rad1=dblarr(27)

x = dblarr(2)
y0 = dblarr(2)
y1 = dblarr(2)

  
;read files
openr,lun,file1,/get_lun
for i = 0, 27-1 do begin
  readf,lun,a,b
  wn(i)=a
  rad0(i)=b
endfor
  free_lun,lun
  
  openr,lun,file2,/get_lun
  for i = 0, 27-1 do begin
    readf,lun,a,b
    rad1(i)=b
  endfor
  free_lun,lun
               
wn = (1/wn)*10000
wn = reverse(wn)
band=where(wn gt 1.9 and wn lt 2.1)
rad0 = reverse(rad0)  

;Cal equivalent width (absorption depth) by Aoki
x = [wn(0),wn(3),wn(5),wn(23),wn(24), wn(25)]
y0 = [rad0(0), rad0(3), rad0(5), rad0(23), rad0(24), rad0(25)]
;y1 = [rad0(0), rad0(25)]
;y1 = [rad1(0), rad1(2)]

coef0 = linfit(X,Y0)
; coef1 = linfit(X,Y1)

cont0 = coef0(0) + coef0(1)*wn
y_calc = 1 - rad0/cont0
y_calc[16] = !VALUES.F_NAN
y_total = total(y_calc[7:18],/nan)


cont0 = coef0(0) + coef0(1)*wn(8)  ; 4837.81250 cm-1  1個目の吸収ライン
cont1 = coef0(0) + coef0(1)*wn(12)  ; 4971.85888 cm-1  2個目の吸収ライン
cont2 = coef0(0) + coef0(1)*wn(16)  ; 5114.22607 cm-1  3個目の吸収ライン
; cont1 = coef1(0) + coef1(1)*wn(1)

    
;                stop
;                window,0
;                plot, wn, rad7 / (coef7(0) + coef7(1)*wn), xs=1, ys=1
;                stop

Table1 = 1.0 - rad0(8)/cont0
Table2 = 1.0 - rad0(12)/cont1
Table3 = 1.0 - rad0(16)/cont2


; Cal absorption depth by Kazama
;y1 = rad0(0)
;y2 = rad0(3)
;y3 = rad0(5)
;y4 = rad0(23)
;y5 = rad0(24)
;y6 = rad0(25)
;yave = (y1+y2+y3+y4+y5+y6)/6
;
;y_nor = rad0/yave
;
;y_total = total(y_nor[7:18])


save,Table_Equivalent_pressure1,$
     filename='/Users/nyonn/IDLWorkspace/Default/LUT/Table_SP_Trans_ver2.sav'
stop
END