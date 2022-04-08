pro nadiaspectrum

  restore, '/Users/nyonn/IDLWorkspace/Default/savfile/orb0278_3.sav'
;;specsol_0403:太陽輝度情報
  openr,2,'/Users/nyonn/IDLWorkspace/Default/profile/'+'specsol_0403.dat'
  specmars=0B
  ;;spcmarsに入れるOMEGAの波長分
  ;格納する場所：specmars
  specmars=fltarr(352)
  readf,2,specmars
  close,2
  
  device, retain=2, decomposed=0,SET_FONT='DejaVuSans', /TT_FONT
  !p.background = 255
  !p.color = 0
  
  specmars=specmars/dmars/dmars
  
;;where関数にxyの２次元が指定できる
;;for文で回して、画像をたくさん表示させることを次回までの宿題で

 ;longitude=where_xyz(longi ge 15.482 and longi le 15.49 and lati ge 9.9 and lati le 10.0,xind=xind,yind=yind)
 longitude2=where_xyz(longi ge 15.485 and longi le 15.49 and lati ge 10.56 and lati le 10.58, xind=xind,yind=yind)
  

  CO2=where(wvl gt 1.8 and wvl lt 2.2)
  co2wvl=wvl(CO2)
  
  ;どうやってjdatの全部を読み込むか。longiとlatiのCO2(75)分に対応する全部
  ;for文で回して全部？それとも格納庫を作って再格納？
  
  longi=longi(co2)
  lati=lati(co2)
  
  ;result=indgen(longi.lati)
  
  ;for i=40, 111 do began
    flux=jdat(*,co2,*)/specmars(co2)
    ;flux(38)=!VALUES.F_NAN
    ;flux(48)=!VALUES.F_NAN
    
    ;y1=flux(0)
    ;y2=flux(5)
    ;y3=flux(10)
    ;y4=flux(15)
    ;y5=flux(55)
    ;y6=flux(60)
    ;y7=flux(65)
    ;y8=flux(70)
    ;yave1=(y1+y2+y3+y4+y5+y6+y7+y8)/8

   ;ConFLUX=flux/yave1
   ;absol=total(flux(band))
   ;
   ;result + = absol
   
  
  flux = jdat(32,CO2,163)/specmars(CO2)
;  flux(38)=!VALUES.F_NAN
 ; flux(48)=!VALUES.F_NAN
  
  flux1 = jdat(31,CO2,199)/specmars(CO2)
;  flux1(38)=!VALUES.F_NAN
 ; flux1(48)=!VALUES.F_NAN
  
  
  
  
 plot, co2wvl,flux1, /nodata;,yrange=[0.2,1.2]
 oplot, co2wvl,flux1, color=1,psym=2
 ;oplot, x, y, color=230, psym=2
 
  
  ;規格化
  y1=flux(0)
  y2=flux(20)
  y3=flux(10)
  y4=flux(15)
  y5=flux(55)
  y6=flux(60)
  y7=flux(65)
  y8=flux(70)
  yave1=(y1+y2+y3+y4+y5+y6+y7+y8)/8
  
  ConFLUX=flux/yave1
  
  y21=flux1(0)
  y22=flux1(20)
  y23=flux1(10)
  y24=flux1(15)
  y25=flux1(55)
  y26=flux1(60)
  y27=flux1(65)
  y28=flux1(70)
  yave2=(y21+y22+y23+y24+y25+y26+y27+y28)/8

  ConFLUX1=flux1/yave2
  
  ;積分
  ;band=where(wvl gt 1.9 and wvl lt 2.1)
  ;band1=wvl(band)
  
  ;absol=total(flux(band))
  
  ;img=image(longi,lati,absol)
  
  ;小暮さんに質問
  ;どうやるの？
  
  x1=co2wvl(0)
  x2=co2wvl(20)
  x3=co2wvl(10)
  x4=co2wvl(15)
  x5=co2wvl(55)
  x6=co2wvl(60)
  x7=co2wvl(65)
  x8=co2wvl(70)
 
x=[x1,x2,x3,x4,x5,x6,x7,x8] 
y=[y21,y22,y23,y24,y25,y26,y27,y28]

loadct,39
;set symbol color  0=black,255=white,254=red
sycolor=0

;pl=plot(co2wvl,ConFLUX,xtitle='wav', ytitle='flux')
;pl2=plot(co2wvl,ConFLUX1,/overplot)
;pl.color="red"
;pl2.color="blue"
  
  
;  
;  plot, co2wvl,flux1, /nodata;,yrange=[0.2,1.2]
;  oplot, co2wvl,flux1, color=1,psym=2
;  ;oplot, x, y, color=230, psym=2
  
  ;snapshot = TVRD(True=1)
  ;Write_JPEG, 'onaji2.jpg', snapshot, true=1

  stop


end


