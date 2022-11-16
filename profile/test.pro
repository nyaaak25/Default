pro test 


ind = 0l
h = 0.
t = 0.
p = 0.
kw0 = fltarr(101101)
kw_SA = fltarr(31,101101)
x = findgen(101101)*0.01 + 4545

file = '/Users/nyonn/IDLWorkspace/Default/savfile/CO2_SP15_TA5_TB3.k'
openr, 1, file
for j = 0, 30 do begin
  readu, 1, ind, h, t, p, kw0
  KW_SA(j,*) = KW0
endfor
close,1
help, ind, h, t, p

kw2 = fltarr(101100)
kw2_K = fltarr(31,101100)
file = '/Users/nyonn/IDLWorkspace/Default/savfile/LUTable_T1_285_T2_200_PRS1500.k'
openr, 1, file
for j = 0, 30 do begin
  readu, 1, ind, h, t, p, kw2
;  help, ind, h, t, p
  KW2_K(j,*) = KW2
endfor
help, ind, h, t, p
close,1

kw3 = fltarr(101100)
kw3_K = fltarr(31,101100)
file = '/Users/nyonn/IDLWorkspace/Default/savfile/test_T1_285_T2_200_PRS1500.k'
openr, 1, file
for j = 0, 30 do begin
  readu, 1, ind, h, t, p, kw3
  ;  help, ind, h, t, p
  KW3_K(j,*) = KW3
endfor
help, ind, h, t, p
close,1

kw4 = fltarr(101100)
kw4_K = fltarr(31,101100)
file = '/Users/nyonn/IDLWorkspace/Default/savfile/test3_T1_285_T2_200_PRS1500.k'
openr, 1, file
for j = 0, 30 do begin
  readu, 1, ind, h, t, p, kw4
  ;  help, ind, h, t, p
  KW4_K(j,*) = KW4
endfor
help, ind, h, t, p
close,1


loadct, 39
Set_plot, 'x'
Device, retain=1, decomposed=0

window, 1
plot, x, KW_SA(30, *), /ylog, back=255, color=0, xs=1
oplot, x, KW_SA(30, *), color=254, thick=2
oplot, x, KW2_K(30, *), color=0
stop

window, 0
plot, x, KW_SA(30, *), /ylog, back=255, color=0, xs=1, xr=[5200,5410]
oplot, x, KW3_K(30, *), color=60, linestyle=2 ;青：風間cutoff120
; oplot, x, KW_SA(30, *), color=254;, thick=3 ; 赤：青木
oplot, x, KW2_K(30, *), color=0, linestyle=2;, thick=3, linestyle=2 ; 黒：風間cutoff300
oplot, x, KW4_K(30, *), color=120, linestyle=2
oplot, x, KW_SA(30, *), color=254 
stop

window, 0
plot, x, KW2_K(30, *), /ylog, back=255, color=0, xs=1, xr=[5400,5410]
oplot, x, KW_SA(30, *), color=254, thick=3
oplot, x, KW2_K(30, *), color=0, thick=3, linestyle=2
stop

window, 1
plot, KW_SA(30, *), /ylog, back=255, color=0, xs=1, xr=[8d4,1d5]
oplot, KW_SA(30, *), color=254, thick=3
oplot, KW2_K(30, *), color=0, thick=3, linestyle=2
stop

window, 1
plot, x, smooth(KW_SA(30, *),1000), /ylog, back=255, color=0, xs=1, /nodata
oplot, x, smooth(KW_SA(30, *),1000), color=254, thick=2
oplot, x, smooth(KW2_K(30, *),1000), color=0
stop

window, 1
plot, KW2_K(30, *)/KW_SA(30, *), back=255, color=0, xs=1

window, 1
plot, smooth(KW2_K(30, *), 1000)/smooth(KW_SA(30, *), 1000), back=255, color=0, xs=1
stop

;
;window, 1
;plot, KW_SA(20, *), /ylog, back=255, color=0, xs=1
;oplot, KW2_K(20, *), color=0
;oplot, KW_SA(20, *), color=254

end