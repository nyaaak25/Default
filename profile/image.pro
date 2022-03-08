pro image

device, retain=2, decomposed=0,SET_FONT='DejaVuSans', /TT_FONT

restore, '/Users/kazama/IDLWorkspace/orb0044_1.sav'
openr,2,'/Users/kazama/IDLWorkspace/specsol_0403.dat'
specmars=0B
specmars=fltarr(352)
readf,2,specmars
close, 2

blue=where(wvl gt 0.476 and wvl lt 0.484)
red=where(wvl gt 0.645 and wvl lt 0.655)
SWIR=wvl(65);(0:127)
LWIR=wvl(150);(128:255)


radred=jdat(*,red,*)
radSW=jdat(*,SWIR,*)
radLW=jdat(*,LWIR,*)

yraw=n_elements(jdat(0,0,*))
xraw=n_elements(jdat(*,0,0))

radred=radred(*,*,0:round(yraw/10)*10-1)
radred=reform(radred,1,xraw,round(yraw/10)*10)
image=rebin(radred,1,64,round(yraw/10))

radSW=radSW(*,*,0:round(yraw/10)*10-1)
radSW=reform(radSW,1,xraw,round(yraw/10)*10)
image1=rebin(radSW,1,64,round(yraw/10))

radLW=radLW(*,*,0:round(yraw/10)*10-1)
radLW=reform(radLW,1,xraw,round(yraw/10)*10)
image2=rebin(radLW,1,64,round(yraw/10))


loadct,0

tvscl, image
tv, image, 65,0

tvscl, image1, 130,0
tv, image1, 195,0

tvscl, image2, 260,0
tv, image2, 325,0


;写真の保存がうまく行かないから、そこをなんとかしよう！
snapshot = TVRD(True=1)
Write_JPEG, 'omegaimage.jpg', snapshot, true=1


stop



end
