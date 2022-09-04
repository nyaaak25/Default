Pro Plo_pressure_results

  restore, '/work1/LUT/SP/table/absorption/work_ORB0920_3.sav'

loadct, 39
Set_Plot, 'x'
device, retain=1, decomposed=0
!Y.OMargin = [10,0]
window,0,xs=800,ys=1400

good = where(pressure ne -999, count)
print, count

pressure = exp(pressure)

min_pressure = min(pressure(good),/nan)
max_pressure = max(pressure(good),/nan)
rn_pressure = max_pressure - min_pressure
plot,longi,lati,xs=1,ys=1,psym=1,back=255,color=0,/nodata,thick=3,ytitle='Latitude',xtitle='Longitude',$
  title='Ls:'+STRCOMPRESS(Solar_longitude)+' Max:'+string(max_pressure,format='(f4.0)' ),charsize=1.;, yr=[50,61], xr=[272,278]
  
for i = 0, io-1 do begin
  for j = 0, ip-1 do begin
    color = (pressure(j,i)-min_pressure)/rn_pressure*254.
    if color gt 254 then color = 254
    if pressure(j,i) eq 0 then color = 255
    if color lt 0 then color = 0
    plots, longi(j,i),lati(j,i),color=color,psym=6,symsize=0.5,thick=3
  endfor
endfor

colorbar_u2,format='(f6.1)',charsize=1.5,charthick=1,range=[min_pressure,max_pressure],DIVISIONS=4,position=[0.8/21., 5/29.7, 1.25/21., 25/29.7],title='Pa',color=0
stop
end
