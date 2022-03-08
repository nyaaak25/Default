pro spectrummap

  restore, '/Users/kazama/IDLWorkspace/orb0278_3.sav'
  ;;specsol_0403:太陽輝度情報
  openr,2,'/Users/kazama/IDLWorkspace/'+'specsol_0403.dat'
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
  
  specmars=specmars[40:111]
  wvl=wvl[40:111]
  nwvl=n_elements(wvl)
  io=n_elements(jdat(*,1,1))
  ip=n_elements(jdat(1,1,*))
  flux=dblarr(io,nwvl,ip)
  for i=0,io-1 do begin
    for o=0,ip-1 do begin
      flux(i,*,o)=jdat(i,*,o)/specmars
      flux(i,38,o)=!VALUES.F_NAN
      flux(i,48,o)=!VALUES.F_NAN
    endfor
  endfor
  
  
  
  ind=where_xyz(longi ge 14.0 and longi le 16.0 and lati ge 9.5 and lati le 13.0,xind=xind,yind=yind)
  band=where(wvl gt 1.9 and wvl lt 2.1)
  band2=wvl(band)
  
  ;where_xyzで探してきたindを入れると、その場所のlongi,latiを取り出してくることができる
  longi1=longi(ind)
  lon2=n_elements(longi1)
  lati1=lati(ind)
  
  ;データ番号と波長のデータセットを作成する
  radiance=dblarr(n_elements(ind),nwvl)
  
  
  io1=n_elements(jdat(xind,1,1))
  ip2=n_elements(jdat(1,1,yind))
  
  wvl1=wvl[29:42]
  band=n_elements(wvl1)
  
  ;配列の作り方
  ;absol=dblarr(io1,band,ip2)
  ;y1=dblarr(io,nwvl,ip)
  ;conflux=dblarr(io1,nwvl,ip2)
  
  y1=flux(7:63,0,168:334)
  y2=flux(7:63,5,168:334)
  y3=flux(7:63,10,168:334)
  y4=flux(7:63,15,168:334)
  y5=flux(7:63,55,168:334)
  y6=flux(7:63,60,168:334)
  y7=flux(7:63,65,168:334)
  y8=flux(7:63,70,168:334)

  yave=(y1+y2+y3+y4+y5+y6+y7+y8)/8
  
  conflux=dblarr(63-7+1,nwvl,334-168+1)
  
  for m=0,nwvl-1 do begin
    conflux1=flux(7:63,m,168:334)/yave
    conflux(*,m,*) = conflux1
  endfor

 absol=total(conflux(*,29:42,*),2,/nan)
 
 
 ;maxabsol=max(absol)
 ;minabsol=min(absol)
 
 ;カラーの値に補間している
 abs_for_color=((absol-minabsol)/(maxabsol-minabsol))*255
 color2=abs_for_color
 
 ;plot, lati(41,168:334),color2(41,*)
 
 ;img=image(color2,rgb_table=33)
 ;cb=colorbar(target=img)
 
 
 ;plots,longi(7:63,168:334),lati(7:63,168:334),color=color2
 
 
 stop
  

  
 end
 