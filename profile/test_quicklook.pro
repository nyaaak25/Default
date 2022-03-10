pro test_quicklook

  restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0006_0.sav'
  ;;specsol_0403:太陽輝度情報
  openr,2,'/Users/nyonn/IDLWorkspace/Default/profile/'+'specsol_0403.dat'
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
  
  file = '/Users/nyonn/IDLWorkspace/Default/orb0313_4.sav'
  fileorbit=strupcase(strmid(file,12,9,/REVERSE_OFFSET))
  path_save='/Users/nyonn/IDLWorkspace/Default/'
  
  specmars=specmars/dmars/dmars

  CO2=where(wvl gt 1.8 and wvl lt 2.2)
  co2wvl=wvl(CO2)

  ;61:89はCO2の範囲内
  specmars=specmars[61:89]
  wvl=wvl[61:89]
  nwvl=n_elements(wvl)
  io=n_elements(jdat(*,1,1))
  ip=n_elements(jdat(1,1,*))
  flux=dblarr(io,nwvl,ip)

  for i=0,io-1 do begin
    for o=0,ip-1 do begin
      flux(i,*,o)=jdat(i,co2,o)/specmars
    endfor
  endfor

  ; 全ての緯度・経度を持ってきて、plotさせる
  maxlongi = max(longi)
  minlongi = min(longi)
  maxlati = max(lati)
  minlati = min(lati)

  ind=where_xyz(longi ge minlongi and longi le maxlongi and lati ge minlati and lati le maxlati,xind=xind,yind=yind)

  ; CO2の吸収線の中
  band=where(wvl gt 1.9 and wvl lt 2.1)
  band2=wvl(band)

  ;where_xyzで探してきたindを入れると、その場所のlongi,latiを取り出してくることができる。今は全範囲。
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

  ;OMEGAの38番目と48番目の素子が死んでた
  nanserch=where(radiance ge 0 and radiance le 0.0001)
  radiance(nanserch)=!VALUES.F_NAN

  ; ======== CO2吸収量 mapping =============
  ; continuumで規格化を行っています
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

  ; カラー値に補間している
  abs_for_color=((absol-minabsol)/(maxabsol-minabsol))*255
  color2=abs_for_color

  ; ====== MOLA高度 mapping =======

  MOLAaltitude=reform(geocube(*,12,*))
  MOLA=double(MOLAaltitude(ind))+2000d

  maxMOLA = max(MOLA)
  minMOLA = min(MOLA)

  MOLA_for_color = double((MOLA-minMOLA)/(maxMOLA-minMOLA))*255
  colorMOLA = MOLA_for_color

  ; ======= MOLA高度 - CO2吸収量 mapping =========

  newcolor = color2-colorMOLA
  maxvalue = max(newcolor)
  minvalue = min(newcolor)
  plotdev = double((newcolor-minvalue)/(maxvalue-minvalue))*255
  colordev = plotdev
  
  ; ========= 他の軌道データについて =============
  
  timei=reform(geocube(0:6,1,*))
  time=timei(*,yind)
  
  local_time=(longi-sub_solar_longitude)*24/360+12
  indLT1=where(local_time lt 0)
  if n_elements(indLT1) gt 1 then local_time(indLT1)=24+local_time(indLT1)
  indLT2=where(local_time ge 24)
  if n_elements(indLT2) gt 1 then local_time(indLT2)=local_time(indLT2)-24
  
  localtime = mean(local_time)

  ; other parameters
  Ls=strmid(SOLAR_LONGITUDE,5,7)
  Latimean=strmid(mean(lati),5,6)
  Longimean=strmid(mean(longi),5,6)
  year=strmid(time(0),8,4)
  month=strmid(time(1),10,2)
  day=strmid(time(2),10,2)
  hour=strmid(time(3),10,2)
  minit=strmid(time(4),10,2)
  
  ; ======== True Color画像 ==============
  
  blue=where(wvl gt 0.476 and wvl lt 0.484)
  red=where(wvl gt 0.645 and wvl lt 0.655)
  
  radred=reform(jdat(*,red,*))
  radblue = jdat(*,blue,*)
  
  redimage =  double(radred(ind))
  
  maxredimage = max(redimage)
  minredimage = min(redimage)

  redimage_for_color = double((redimage-minredimage)/(maxredimage-minredimage))*255
  colorredimage = redimage_for_color
  
  yraw=n_elements(jdat(0,0,*))
  xraw=n_elements(jdat(*,0,0))
  
;  radred=radred(*,*,0:round(yraw/10)*10-1)
;  radred=reform(radred,1,xraw,round(yraw/10)*10)
;  red_image=rebin(radred,1,64,round(yraw/10))
  
  radblue=radblue(*,*,0:round(yraw/10)*10-1)
  radblue=reform(radblue,1,xraw,round(yraw/10)*10)
  blue_image=rebin(radblue,1,64,round(yraw/10))
  
;  image_mars = red_image + blue_image
  
  ; ========= plotする場所 ================
  
  !p.multi = [0,3,1,0,0]

  ; CO2吸収量 mapping
;  plot,longi(ind),lati(ind), xstyle=1,ystyle=1,title='CO2_absorption',xtitle='latitude'$
;    ,color=0,position=[0.05,0.13,0.3,0.92],xticks=2,/nodata, charsize=2
;  ; cgColorbar,Range=[10,10.5],position=[0.01,0.02,0.3,0.03]
;  plots, longi(ind),lati(ind),color=color2,psym=2
;  
;  ; 日付と軌道番号
;  xyouts,0.06,0.975,fileorbit,charsize=1.5,color=0,/normal
;  xyouts,0.25,0.975,'Time : '+year+'/'+month+'/'+day+' '+hour+':'+minit,color=0,/normal,charsize=1
;  xyouts,0.45,0.975,'Ls : '+Ls,color=0,/normal,charsize=1
;  xyouts,0.65,0.975,'Local Time : '+localtime,color=0,/normal,charsize=1
;  
;
;  ; MOLA高度 mapping
;   plot,longi(ind),lati(ind),xstyle=1,ystyle=1,title='MOLA altitude',xtitle='latitude'$
;    ,color=0,position=[0.35,0.13,0.6,0.92],xticks=2,/nodata, charsize=2
;  cgColorbar,Range=[10,10.5],position=[0.2,0.03,0.7,0.04]
;  plots, longi(ind),lati(ind),color=colorMOLA,psym=2
;
;  ; CO2吸収量 - MOLA高度
;  plot,longi(ind),lati(ind),xstyle=1,ystyle=1,title='CO2 - MOLA altitude',xtitle='latitude'$
;    ,color=0,position=[0.65,0.13,0.9,0.92],xticks=2,/nodata, charsize=2
;  ; cgColorbar,Range=[10,10.5],position=[0.65,0.01,0.9,0.03]
;  plots, longi(ind),lati(ind),color=colordev,psym=2
  
  
  plot, longi(ind),lati(ind),xstyle=1,ystyle=1,title='CO2 - MOLA altitude',xtitle='latitude'$
      ,color=0,position=[0.65,0.13,0.9,0.92],xticks=2,/nodata, charsize=2
  
  plots, longi(ind),lati(ind),color=colorredimage,psym=2
  
;  snapshot = TVRD(True=1)
;  Write_JPEG, path_save+fileorbit+'.jpg', snapshot, True=1, Quality=100
 
 stop


end