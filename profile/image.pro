pro image

device, retain=2, decomposed=0, SET_FONT='DejaVuSans', /TT_FONT
!p.background = 255

restore, '/Users/nyonn/IDLWorkspace/Default/savfile/ORB0024_1.sav'
openr,2,'/Users/nyonn/IDLWorkspace/Default/profile/'+'specsol_0403.dat'
specmars=0B
specmars=fltarr(352)
readf,2,specmars
close, 2

;blue=where(wvl gt 0.476 and wvl lt 0.484)
;red=where(wvl gt 0.645 and wvl lt 0.655)
;green=where(wvl gt 0.569 and wvl lt 0.575)
blue=where(wvl gt 0.432 and wvl lt 0.44)
red=where(wvl gt 0.695 and wvl lt 0.7)
green=where(wvl gt 0.54 and wvl lt 0.55)
SWIR=wvl(65);(0:127)
LWIR=wvl(150);(128:255)


radred=jdat(*,red,*)
radblue=jdat(*,blue,*)
radgreen=jdat(*,green,*)
radSW=jdat(*,SWIR,*)
radLW=jdat(*,LWIR,*)

yraw=n_elements(jdat(0,0,*))
xraw=n_elements(jdat(*,0,0))

;radred=radred(*,*,0:round(yraw/10)*10-1)
;radred=reform(radred,1,xraw,round(yraw/10)*10)
;red_image=rebin(radred,1,64,round(yraw/10))

radred=radred(*,*,0:yraw-1)
radred=reform(radred,1,xraw,yraw)
; red_image=rebin(radred,1,64,yraw)

radblue=radblue(*,*,0:round(yraw/10)*10-1)
radblue=reform(radblue,1,xraw,round(yraw/10)*10)
blue_image=rebin(radblue,1,64,round(yraw/10))

radgreen=radgreen(*,*,0:round(yraw/10)*10-1)
radgreen=reform(radgreen,1,xraw,round(yraw/10)*10)
green_image=rebin(radgreen,1,64,round(yraw/10))

radSW=radSW(*,*,0:round(yraw/10)*10-1)
radSW=reform(radSW,1,xraw,round(yraw/10)*10)
image1=rebin(radSW,1,64,round(yraw/10))

radLW=radLW(*,*,0:round(yraw/10)*10-1)
radLW=reform(radLW,1,xraw,round(yraw/10)*10)
image2=rebin(radLW,1,64,round(yraw/10))

;image1 = red_image ;+ blue_image  + green_image
image1 = radred ;+ (radblue*0.7) + (radgreen*0.9)

loadct,39

;xsize= 1000
;ysize= 1000
;window,xsize=xsize,ysize=ysize

tvscl, image1, 0,0
tv, image1, 300,0

;tvscl, image1, 0,100
;tv, image1, 195,0
;
;tvscl, image2, 260,0
;tv, image2, 325,0


; 赤、青、緑の配色をうまくいれていきたい！

;写真の保存がうまく行かないから、そこをなんとかしよう！
snapshot = TVRD(True=1)
; Write_JPEG, 'omegaimage.jpg', snapshot, true=1


stop



end
