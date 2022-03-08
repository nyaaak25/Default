FUNCTION READOMEGA, FORCE = force

IF ~KEYWORD_SET(force) THEN force = 0L

;++++++++++++++++++++++++
;Written by SHOHEI AOKI
;++++++++++++++++++++++++
;
;================================================
;path
;================================================
;pathnames = PATH_MANAGEMENT()

dirsoft='/Users/nyonn/IDLWorkspace/Default/'
dirdata ='/Users/nyonn/IDLWorkspace/Default/'

datapath = dirdata
geompath = dirdata

;================================================
;search file
;================================================

files_nav = FILE_SEARCH(dirdata + '*.nav', COUNT = count)
i_nav = count

files_qub = FILE_SEARCH(dirdata + '*.qub', COUNT = count)
i_qub = count

if i_nav ne i_qub then begin
  print, '.nav and .qub is not consistent'
  stop
endif

;================================================
;Loop start
;================================================
FOR Loop = 0, i_qub - 1 DO BEGIN
  
  IF ~force THEN BEGIN
    sdir = FILE_SEARCH(geompath + FILE_BASENAME(files_qub(loop),'.qub') + '*.sav', COUNT = tnf)
    IF tnf GT 0L THEN CONTINUE
  ENDIF
  
  file_qub = files_qub(loop)
  file_nav = files_nav(loop)
  
  PRINT, 'Loop', STRING(Loop, "(I04)"), '/', STRING(i_nav, "(I04)"), ', file: ', file_qub

  dirsave = FILE_DIRNAME(file_qub)

  nomfic0 = file_basename(file_qub,'.qub')
  orbnum=fix(strmid(nomfic0,3,4))
  nomgeo=nomfic0+'.nav'
  nomfic0=nomfic0+'.qub'
  nomfic=datapath+nomfic0
  openr,2,nomfic,ERROR=errflag
  close,2
;  if (errflag ne 0) then begin
;    print,'file ',nomfic,' not found'
;    goto, start
;  endif
  readcube, nomfic, idat, sdat0, sdat1,info
  if(size(idat))(0) ne 3 then begin
    print,'***** only one line in cube ',nomfic0,' ******'
    CONTINUE
  endif
  exposure=0
  exposure=info(0:2)
  summation=info(3)
  bits_per_data=info(4)
  nomgeo=geompath+nomgeo

  openr,1,nomgeo,ERROR=errflag
  if(errflag ne 0) then begin
    print,' no corresponding nav cube'
    close,1
    dmars=1.52
    specmars=fltarr(352)
    openr,2,dirsoft + 'specsol_0403.dat'
    readf,2,specmars
    close,2
    specmars=specmars/dmars/dmars
    goto, nogeom
  endif
  trans= !pi/180.*1.e-4

  geocube=0B
  ilec=0B

  read_geolbl,nomgeo, '', values, data

  lrec   = data(0)
  nrec   = data(1)
  npixel = data(2)
  npara  = data(3)
  nscan  = data(4)
  nau    = data(5)
  solar_longitude = data(8)
  sub_solar_longitude = data(6)
  dmars=nau*1.e-4
  close,2
  openr,2,dirsoft + 'specsol_0403.dat'
  specmars=0B
  specmars=fltarr(352)
  readf,2,specmars
  close,2
  specmars=specmars/dmars/dmars

  data = 0B
  geocube=0B

  geocube = lonarr(npixel, npara, nscan)
  ilec = lonarr(npixel)

  point_lun,1,lrec*nrec

  for  k=0,nscan-1 do begin
    for j=0,npara-1 do begin
      readu, 1, ilec
      geocube(*,j,k) = ilec
    endfor
  endfor
  close, 1
  if(geocube(1,1,0) gt 13) then geocube=swap_endian(temporary(geocube))
  mirror_position=reform(geocube(*,0,*))
  timei=reform(geocube(0:6,1,*))
  longi=reform(geocube(*,6,*)*1.e-4)
  lati=reform(geocube(*,7,*)*1.e-4)
  slant=reform(geocube(*,11,*)*1.e-3)
  alt=reform(geocube(*,12,*)*1.e-3)
  nogeom:
  close,/all
  ; preliminary pipeline tool

  a=size(idat)
  nbal=a(3)
  npix=a(1)
  jdat=0
  jdat=float(idat)

  fond2=intarr(256)
  openr,2,dirsoft + 'fond2.dat'
  readf,2,fond2
  close,2
  if(exposure(0) gt 4.) then fond2=2*fond2
  fondcur=fond2(128:255)#(0.*indgen(nbal)+1.)

  jdat=0
  jdat=float(idat)
  i=where(idat le 1)
  if(i(0) ne -1) then jdat(i)=1.e-5

  pix0IR=0
  i=where(idat(*,0:255,*) le 0)
  if(i(0) ne -1) then pix0IR=(size(i))(3)
  print,'       0 or less  IR: ',pix0IR

  hkmin=min(sdat1(14,1,*))
  if(hkmin lt 6) then hkmin=6
  indj=where(sdat1(14,1,*) ge hkmin)
  i6=indj(0)
  indj=where(sdat1(14,1,*) ge hkmin+1)
  balHK=indj(0)-i6
  indi=where(sdat1(14,1,*) ge 6 and sdat0(10,*) gt 100)
  balsm=balHK*8

  if(balHK gt 0 and nbal gt indi(0)+balsm) then begin
    a=0.
    b=0.
    c=0.
    d=0.
    ndeb=indi(0)
    nf=(size(indi))(1)-1
    nfin=indi(nf)
    for k=0,255 do begin
      b=reform(float(sdat0(k,indi)))
      a=[2*b(0)-rotate(b(0:balsm-1),2),b,2*b(nf)-rotate(b(nf-balsm+1:nf),2)]
      c=(smooth(a,balsm))(balsm:balsm+nf)
      d= -sdat0(k,ndeb:nbal-1)+spline(indi,c,ndeb+indgen(nbal-ndeb))
      for i=0,npix-1 do jdat(i,k,ndeb:nbal-1)= $
        jdat(i,k,ndeb:nbal-1)+d
    endfor
  endif

  for i=0,npix-1 do jdat(i,128:255,*)=jdat(i,128:255,*)-fondcur

  i=where(jdat lt 1.e-5)
  if(i(0) ne -1) then jdat(i)=1.e-5

  if(summation ne 1) then jdat=jdat/summation
  linearC=0
  linearC=fltarr(4096)
  openr,2,dirsoft + 'linearC.dat'
  readf,2,linearC
  close,2
  jdat(*,0:127,*)=linearC(fix(jdat(*,0:127,*)+0.5))
  wvl=0
  wvl=fltarr(352)
  openr,2,dirsoft + 'lambda_0403.dat'
  readf,2,wvl
  close,2
  mtf=fltarr(352)
  rap=fltarr(256)
  bound=intarr(256)
  openr,2,dirsoft + 'bound060110.dat'
  readf,2,bound
  close,2
  if(exposure(0) lt 4.) then begin
    openr,2,dirsoft + 'mtf060110_25.dat'
    readf,2,mtf
    close,2
    openr,2,dirsoft + 'rap060110_25.dat'
    readf,2,rap
    close,2
  endif else begin
    openr,2,dirsoft + 'mtf060110_50.dat'
    readf,2,mtf
    close,2
    openr,2,dirsoft + 'rap060110_50.dat'
    readf,2,rap
    close,2
  endelse
  ib=where(orbnum gt bound)
  mtf(ib)=mtf(ib)*rap(ib)
  for n=0,255 do jdat(*,n,*)=jdat(*,n,*)/mtf(n)
  ic=[where(rap lt 10.),256+indgen(96)]

  ;**************************************************************************************************************************
  ;bit error correction
  ;**************************************************************************************************************************

  vis=float(idat(*,256:351,*))
  siz = size(idat)
  lines = siz(3)
  pixels = siz(1)
  exptime = info(2)/1000.
  image_ratio = fltarr(pixels, lines)

  level = 4095.*info(3)
  if pixels EQ 128 then level = info(3)*4095*2.

  i=where(vis le 0)
  counter_neg=0
  if(i(0) ne -1) then begin
    counter_neg=(size(i))(3)
    vis(i)=0.001
  endif

  i=where(vis gt level)
  counter_pos=0
  if(i(0) ne -1) then begin
    counter_pos=(size(i))(3)
    vis(i)=level
  endif
  i=where(vis le level and vis gt 0.8*level)
  counter_sat=0
  if(i(0) ne -1) then begin
    counter_sat=(size(i))(3)
    vis(i)=level*0.8
  endif

  plan3=reform(vis(*,3,*))
  counter_spike3=0
  i=where(plan3 gt 0.5*(vis(*,2,*)+vis(*,3,*))+50)
  if(i(0) ne -1) then begin
    counter_spike3=(size(i))(3)
    plan3(i)=0.5*(reform(vis(*,2,*)+vis(*,3,*)))(i)
    vis(*,3,*)=plan3
  endif

  median_line = fltarr(lines)
  column = fltarr(lines)

  counter_spikes=long(0)

  if(lines lt 10) then goto, nodespike

  for k = 0, 95 do begin
    for m = 0, pixels - 1 do begin
      column(*) = vis(m, k, *)
      median_line = median(column, 7)
      for i = 0, lines-1 do begin
        if (abs(median_line(i) - vis(m, k, i)) GT 200) then begin
          vis(m, k, i) = median_line(i)
          counter_spikes = counter_spikes + 1
        endif
      endfor
    endfor
  endfor
  nodespike:

  print, ' negative pixels VIS:', counter_neg
  print, 'anomalous pixels VIS:', counter_pos
  print, 'saturated pixels VIS:', counter_sat
  print, '          spikes VIS:', counter_spikes

  vis = float(fix(vis))

  ;******************************************************************************************************************************
  ; bias correction
  ;******************************************************************************************************************************


  bias = fltarr(2, 20000)
  openr, 1, dirsoft + 'bias010705.txt'
  readf, 1, bias
  close, 1

  mean_ch_mat = fltarr(lines)
  mean_sl_mat = fltarr(lines)
  for n = 0, lines - 1 do begin
    mean_sl = total(vis(*, *, n), 1) / pixels
    mean_sl_mat(n) = mean_sl(21)
    mean_ch_mat(n) = mean_sl(58)
    area_sl = total(mean_sl(0:84), 1) / 85.

    for i = 0, pixels-1 do begin
      vis(i, *, n) = offset_corr_050701( vis(i, *, n), mean_sl, info(3), pixels, bias)
    endfor
  endfor

  ;******************************************************************************************************************************
  ; smear correction
  ;******************************************************************************************************************************

  for i = 0, pixels-1 do begin
    for n = 0, lines-1 do begin
      vis(i, *, n) = smear_corr_050701( vis(i, *, n), info(2))

    endfor
  endfor


  ;******************************************************************************************************************************
  ;flat cube
  ;******************************************************************************************************************************

  slice = ulonarr(128, 96)

  openr, 1, dirsoft + 'flatVIS050701.bin'
  readu, 1, slice
  close, 1
  if(slice(48,0) gt 40000) then slice=swap_endian(slice)


  istart= (128 - pixels)/2.
  iend  = istart + pixels - 1

  for i = 0, lines-1 do begin

    vis(*, *, i) = vis(*, *, i) * 2.^15 / slice(istart:iend, *)

  endfor

  ;***************************************************************************************************************************
  ;radiometric calibration
  ;***************************************************************************************************************************

  f = 1.
  if pixels eq 128 then  f = 2; internal summation
  f2 = 1.
  if exptime eq 0.1 then f2 = 2.
  if exptime eq 0.2 then f2 = 4.

  spectrum = fltarr(96)

  for i = 0, lines - 1 do begin
    for m = 0, pixels - 1 do begin
      spectrum = vis(m, *, i)
      jdat(m, 256:351, i) = ordcorr(wvl(256:351), spectrum)/ ( f2 * mtf(256:351)) / f / summation
    endfor
  endfor
  
  ;================================================
  ; save data
  ;================================================
  save, wvl, jdat, lati, longi, geocube, dmars, Solar_longitude, SUB_SOLAR_LONGITUDE, timei, mirror_position, alt, slant, filename=dirsave+'/'+file_basename(file_qub,'.qub')+'.sav'

endfor

RETURN, 0

end
