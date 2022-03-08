pro newOMEGA

restore, '/Users/kazama/IDLWorkspace/orb0044_1.sav'

openr,2,'/Users/kazama/IDLWorkspace/'+'specsol_0403.dat'
specmars=0B
specmars=fltarr(352)
readf,2,specmars
close,2
specmars=specmars/dmars/dmars

maxaltitude=60
alt2 = alt-65.536
ind=where_xyz(alt2 ge 0 and alt2 le maxaltitude+0.1,xind=xind,yind=yind)

CO2=where(wvl gt 1.8 and wvl lt 2.2)
co2wvl=wvl(CO2)

flux = jdat(0,CO2,0)/specmars(CO2)

;flux(38)=!VALUES.F_NAN
;flux(48)=!VALUES.F_NAN


plot, co2wvl, flux,color=1,psym=2,title='CO2_absorption',xtitle='wavelength',ytitle='reflectance'

;plot, co2wvl,flux, /nodata;,yrange=[0.2,1.2]
;oplot, co2wvl,flux, color=1,psym=2

stop


end

;;小暮さんへの質問
;;Q1.NANを代入するのはどうやる？⇒定義が必要。なんかいろんな宣言がいる。
;;   NANデータにするのなら、なにで補完？線形？
;;Q2.これはある1点のデータを持ってきているイメージがある。これをある高度(nadia ob)全部の足し合わせでCO2のラインを見たいときはどうやればよいでしょうか。。
;;　　地表から大気まで全部のデータのCO2吸収帯を見たい。高度指定のやり方を教えてほしい。。

