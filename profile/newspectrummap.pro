pro NEWspectrummap

  restore, '/Users/nyonn/IDLWorkspace/Default/orb0313_4.sav'
  ;;specsol_0403:太陽輝度情報
  openr,2,'/Users/nyonn/IDLWorkspace/Default/'+'specsol_0403.dat'
  specmars=0B
  ;;spcmarsに入れるOMEGAの波長分
  ;格納する場所：specmars
  specmars=fltarr(352)
  readf,2,specmars
  close,2

  ;カラーの入れ方
  device, retain=2, decomposed=0,SET_FONT='DejaVuSans', /TT_FONT
  !p.background = 255
  !p.color = 0
  loadct,39


  specmars=specmars/dmars/dmars

  CO2=where(wvl gt 1.8 and wvl lt 2.2)
  co2wvl=wvl(CO2)

  specmars=specmars[61:89]
  wvl=wvl[61:89]
  nwvl=n_elements(wvl)
  io=n_elements(jdat(*,1,1))
  ip=n_elements(jdat(1,1,*))
  flux=dblarr(io,nwvl,ip)
  for i=0,io-1 do begin
    for o=0,ip-1 do begin
      flux(i,*,o)=jdat(i,co2,o)/specmars
      ;flux(i,38,o)=!VALUES.F_NAN
      ;flux(i,48,o)=!VALUES.F_NAN
    endfor
  endfor

  maxlongi = max(longi)
  minlongi = min(longi)
  maxlati = max(lati)
  minlati = min(lati)
  
  ind=where_xyz(longi ge minlongi and longi le maxlongi and lati ge minlati and lati le maxlati,xind=xind,yind=yind)
  band=where(wvl gt 1.9 and wvl lt 2.1)
  band2=wvl(band)

  ;where_xyzで探してきたindを入れると、その場所のlongi,latiを取り出してくることができる
  longi1=longi(ind)
  lon2=n_elements(longi1)
  lati1=lati(ind)

  ;データ番号と波長のデータセットを作成する
  radiance=dblarr(n_elements(ind),nwvl)
  
  ;reformは次元を落としてくれるよ
  for i=0,nwvl-1 do begin
    newflux1=reform(flux(*,i,*))
    radiance(*,i)=double(newflux1(ind))
  endfor
  
  nanserch=where(radiance ge 0 and radiance le 0.0001)
  radiance(nanserch)=!VALUES.F_NAN
; plot, radiance(0,*)
 
 y1=radiance(*,0)
 y2=radiance(*,3)
 y3=radiance(*,5)
 y4=radiance(*,23)
 y5=radiance(*,24)
 y6=radiance(*,25)
 yave=(y1+y2+y3+y4+y5+y6)/6
 
 conflux=dblarr(n_elements(ind),nwvl)
 
 for m=0,nwvl-1 do begin
   conflux1=radiance(*,m)/yave
   conflux(*,m) = conflux1
 endfor
 
 absol=total(conflux(*,band),2,/nan)
 
 maxabsol=max(absol)
 minabsol=min(absol)
 
 ;カラーの値に補間している
 abs_for_color=((absol-minabsol)/(maxabsol-minabsol))*255
 color2=abs_for_color
; plot,longi(ind),lati(ind), xstyle=1,ystyle=1,title='CO2_absorption',xtitle='latitude'$
;   ,color=0,position=[0.05,0.13,0.5,0.95],xticks=2,/nodata
; cgColorbar,Range=[10,10.5],position=[0.05,0.05,0.5,0.07]
; plots, longi(ind),lati(ind),color=color2,psym=2 

 snapshot = TVRD(True=1)
; Write_JPEG, 'Spiga_case1_CO2absorption.jpg', snapshot, true=1

;ran=where_xyz(longi ge 15.45 and longi le 15.55 and lati ge 9.5 and lati le 13.0,xind=xind,yind=yind)
;range=where(longi ge 15.4 and longi le 15.6)
;plot,lati(ind),color2,xrange=[9.5,13],yrange=[0,250],xstyle=1,title='CO2absorption longi:15.48-15.52',xtitle='latitude'
;xstyleは軸を消したり、色々できる。x=1は強制的にxrangeを表す

MOLAaltitude=reform(geocube(*,12,*))
MOLA=double(MOLAaltitude(ind))+2000d

aveMOLA=mean(MOLA)
avecolor=mean(color2)

maxMOLA=max(MOLA)
minMOLA=min(MOLA)

scalling=aveMOLA/avecolor
newcolor=color2*scalling

;----MOLA+CO2absopyion altitude--------
;plot,lati(ind),MOLA,/nodata,title='MOLA&CO2_absoption',xrange=[9.5,13],yrange=[-2800,10],xstyle=1,ystyle=1,xtitle='latitude'
;plot,lati(ind),MOLA,/nodata,title='MOLA&CO2_absoption',xstyle=1,ystyle=1,xtitle='latitude'
;oplot,lati(ind),MOLA,color=1
;oplot,lati(ind),newcolor,color=230
;xyouts,0.8,0.88,'Red=CO2absorption',charsize=1.5,color=230,/normal
;xyouts,0.8,0.9,'MOLAaltitude',charsize=1.5,color=1,/normal

;---MOLA topography---------
MOLA_for_color=double((MOLA-minMOLA)/(maxMOLA-minMOLA))*255
colorMOLA=MOLA_for_color
;plot,longi(ind),lati(ind),xrange=[24.75,25.2],yrange=[-10.1,-4.8],xstyle=1,ystyle=1,title='MOLA altitude',xtitle='latitude'$
;,color=0,position=[0.05,0.13,0.5,0.95],xticks=2,/nodata
;plot,longi(ind),lati(ind),xstyle=1,ystyle=1,title='MOLA altitude',xtitle='latitude'$
; ,color=0,position=[0.05,0.13,0.5,0.95],xticks=2,/nodata
;cgColorbar,Range=[10,10.5],position=[0.05,0.05,0.5,0.07]
;plots, longi(ind),lati(ind),color=colorMOLA,psym=2
;
;snapshot = TVRD(True=1)
;Write_JPEG, 'Spiga_case1_MOLAto.jpg', snapshot, true=1


;------plot差異------
newcolor=color2-colorMOLA
maxvalue=max(newcolor)
minvalue=min(newcolor)
plotdev=double((newcolor-minvalue)/(maxvalue-minvalue))*255
colordev=plotdev
;
plot,longi(ind),lati(ind),xstyle=1,ystyle=1,title='CO2 - MOLA altitude',xtitle='latitude'$
 ,color=0,position=[0.05,0.13,0.7,0.95],xticks=2,/nodata
cgColorbar,Range=[10,10.5],position=[0.05,0.05,0.7,0.07]
plots, longi(ind),lati(ind),color=colordev,psym=2

;--------OMEGA data time------
time=timei(*,yind)

;local_time=(geocube(*,6,*)*1.e-4-sub_solar_longitude)*24/360+12
;indLT1=where_xyz(local_time lt 0,xind=xind,yind=yind,zind=zind)
;if n_elements(indLT1) gt 1 then local_time(xind,yind,zind)=24+local_time(xind,yind,zind)
;indLT2=where_xyz(local_time ge 24,xind=xind,yind=yind,zind=zind)
;if n_elements(indLT2) gt 1 then local_time(xind,yind,zind)=local_time(xind,yind,zind)-24

local_time=(longi-sub_solar_longitude)*24/360+12
indLT1=where(local_time lt 0)
if n_elements(indLT1) gt 1 then local_time(indLT1)=24+local_time(indLT1)
indLT2=where(local_time ge 24)
if n_elements(indLT2) gt 1 then local_time(indLT2)=local_time(indLT2)-24

; [other parameters]
Ls=strmid(SOLAR_LONGITUDE,5,7)
Latimean=strmid(mean(lati),5,6)
Longimean=strmid(mean(longi),5,6)
year=strmid(time(0),8,4)
month=strmid(time(1),10,2)
day=strmid(time(2),10,2)
hour=strmid(time(3),10,2)
minit=strmid(time(4),10,2)
;
;plot,lati,longi,position=[0,0.83,0.95,1],color=255
;;xyouts,0.03,0.975,fileorbit,charsize=2.5,color=0,/normal
;xyouts,0.03,0.955,'Time : '+year+'/'+month+'/'+day+' '+hour+':'+minit,color=0,/normal,charsize=1.5
;xyouts,0.03,0.94,'Ls : '+Ls,color=0,/normal,charsize=1.5
;xyouts,0.03,0.89,'Local Time : '+local_time,color=0,/normal,charsize=1.5
;    xyouts,0.03,0.84,'Latitude : '+Latimean,color=0,/normal,charsize=2
;    xyouts,0.03,0.81,'Longitude : '+Longimean,color=0,/normal,charsize=2

;snapshot = TVRD(True=1)
;Write_JPEG, 'CO2absorption.jpg', snapshot, true=1

  stop



end
