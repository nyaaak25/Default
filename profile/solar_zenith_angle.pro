pro solar_zenith_angle
;
;観測地点の太陽天頂角を知るために、太陽緯度経度を計算し、太陽天頂角を計算するプログラム。 Kogure 2021/3/19
;

pathfile='/Users/nyonn/IDLWorkspace/Default/orb0072_2.sav'
pathKernels='/Users/nyonn/IDLWorkspace/Default/kernels'

cspice_kclear;カーネルをloadする空きスペースを消去して確保
cspice_furnsh, '/Users/nyonn/IDLWorkspace/Default/kernels/OMEGA.mk';loadする
;
;cspice_furnsh, 'c:¥SPICE¥kernel¥lsk¥naif0010.tls' 
;if getkernelsspan
;tocheck=(pathKernels+'spk');
;for ii=1,n_elements(tocheck) do begin
;  ligne=strcmpi(kernelslist(*,1),tocheck(ii));
;  for i=1,n_elements(kernelslist(ligne,2)) do begin
;    ids=cspice_spkobj([pathKernels+tocheck(ii)+filesep+kernelslist(ligne,2)(i)],10000);
;    fprintf('File %s - %s\n',tocheck(ii),kernelslist{ligne,2}{i});
;    for j=1,n_elements(ids) do begin
;      cover=cspice_spkcov([pathKernels+tocheck(ii)+filesep+kernelslist(ligne,2)(i)],ids(j),10000);
;      fprintf('ID %i:\n',ids(j))
;      for k=1:2:n_elements(cover) do begin
;        timestrstart=cspice_timout(cover(k),'YYYY MON DD HR:MN:SC.### (TDB) ::TDB');
;        timestrend=cspice_timout(cover(k+1),'YYYY MON DD HR:MN:SC.### (TDB) ::TDB');
;        fprintf('%s to %s\n',timestrstart,timestrend);
;      endfor
;    endfor
;  endfor
;endfor

;'pathfiles'の下に存在する.savファイルの検索
;files=file_search(pathfile+'*.sav',count=count)
;for loop=185,count-1 do begin;;;;;;;;loop=185, orb4483
;  print,loop,'/',count
;  loadct,39


;ファイルのリストア
  ;file=files(loop)
  restore,pathfile
 
  
;求めたい高度域のデータの抽出
  alt=alt-65.536
  maxaltitude=150
  ind=where_xyz(longi ge 178 and longi le 202 and lati ge 10 and lati le 45,xind=xind,yind=yind)
  latipi=lati(ind)/180*!pi
  longipi=longi(ind)/180*!pi
  nind=n_elements(ind)


;太陽緯度経度を求めるための準備として、観測時間の取得  
  time=timei(*,yind)
  month=string(time(1,*))
  month=month.Replace('1','JAN')
  month=month.Replace('2','FEB')
  month=month.Replace('3','MAR')
  month=month.Replace('4','APR')
  month=month.Replace('5','MAY')
  month=month.Replace('6','JUN')
  month=month.Replace('7','JUL')
  month=month.Replace('8','AUG')
  month=month.Replace('9','SEP')
  month=month.Replace('10','OCT')
  month=month.Replace('11','NOV')
  month=month.Replace('12','DEC')
  
  month=strmid(month,11,3)
  day=time(2,*)
  year=time(0,*)
  time=time(3:6,*)
  day=string(day,format='%02i')
  time=string([year,time(0,*),time(1,*),time(2,*),time(3,*)],format='%04i %02i:%02i:%02i.%03i')
  data=day+' '+month+' '+time;cspice_str2etに突っ込む時間情報のstring
  

;準備した観測時間から太陽緯度経度を取得。関数は検索すると詳細が見れる。  
  etcube=dblarr(nind)
  spoint1=dblarr(3,nind)
  trgepc1=dblarr(nind)
  srfvec1=dblarr(3,nind)
  radius1=dblarr(nind)
  lonsol=dblarr(nind)
  latsol=dblarr(nind)
  spoint_mex=dblarr(3,nind)
  trgepc_mex=dblarr(nind)
  srfvec_mex=dblarr(3,nind)
  radius_mex=dblarr(nind)
  lon_mex=dblarr(nind)
  lat_mex=dblarr(nind)
  
  for i=0,nind-1 do begin
    cspice_str2et, data(i), et;epocデータに変換(時間型のデータ)
    etcube(i)=et
    cspice_subpnt, 'NEAR POINT: ELLIPSOID','MARS',etcube(i),'IAU_MARS','NONE','SUN',spoint, trgepc, srfvec;火星からみた太陽の位置を直交座標で
    cspice_subpnt, 'NEAR POINT: ELLIPSOID','MARS',etcube(i),'IAU_MARS','NONE','MEX',spoint_mex, trgepc_mex, srfvec_mex;火星からみたMExの位置を直交座標で
    spoint1(*,i)=spoint
    trgepc1(i)=trgepc
    srfvec1(*,i)=srfvec
    spoint_mex(*,i)=spoint_mex
    trgepc_mex(i)=trgepc_mex
    srfvec_mex(*,i)=srfvec_mex
    cspice_reclat, spoint, radius, longitude, latitude;曲座標に変換
    cspice_reclat, spoint_mex, radius_mex, longitude_mex, latitude_mex;曲座標に変換
    radius1(i)=radius
    lonsol(i)=longitude
    latsol(i)=latitude
    radius_mex(i)=radius_mex
    lon_mex(i)=longitude_mex
    lat_mex(i)=latitude_mex
  
;求めた太陽緯度経度と観測点の緯度経度から、太陽天頂角を計算。
  ;cspice_subpnt('NEAR POINT: ELLIPSOID','MARS',etcube(i),'IAU_MARS','NONE','SUN');
  cosd_sol=sin(latsol)*sin(latipi)+cos(latsol)*cos(latipi)*cos(longipi-lonsol);%太陽天頂角
  sza=acos(cosd_sol)
  sza=sza*180/!pi;[degree]
  
  print, '太陽天頂角', sza[10]
 endfor  
end