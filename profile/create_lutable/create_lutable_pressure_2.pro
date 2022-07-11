function leggi_file, file, colonne, head_size

  if head_size gt 0 then header=strarr(head_size)
  nlines=file_lines(file)-head_size
  data=dblarr(colonne, nlines)
  openr, lun, file, /get_lun
  if head_size gt 0 then readf, lun, header
  readf, lun, data
  free_lun, lun

  return, data
end

function phi_computing, I, E, G

theta=E
theta0=I
phase=G

p = !pi/180.
if (theta0 ne 0 and theta ne 0) then begin
        phi=(cos(phase*p)-cos(theta0*p)*cos(theta*p))/(sin(theta0*p)*sin(theta*p))
   ;print, phi
        phi=min([1d0,phi])
;        ;print, phi
        phi=max([-1d0,phi])
;    print, phi
        phi=acos(phi)/p
      endif else begin
        phi=0
      endelse

if (finite(PHI) eq 0) then PHI=0.0

return, phi
end

;------------------------------------------------------------------------------------------------------------------------
Pro Create_lutable_pressure_2
;------------------------------------------------------------------------------------------------------------------------

cd, '/home/nyonnkazama/' ;working dir
path_LUT_Input = '/work1/LUT/SP/table/input/'
path_LUT_Output = '/work1/LUT/SP/table/output/'

ISP = 2 ;1-15 for different SP

;parameters
IF ISP eq 1  then SP = 50.0d+0
IF ISP eq 2  then SP = 150.0d+0
IF ISP eq 3  then SP = 180.0d+0
IF ISP eq 4  then SP = 215.0d+0
IF ISP eq 5  then SP = 257.0d+0
IF ISP eq 6  then SP = 308.0d+0
IF ISP eq 7  then SP = 369.0d+0
IF ISP eq 8  then SP = 442.0d+0
IF ISP eq 9  then SP = 529.0d+0
IF ISP eq 10 then SP = 633.0d+0
IF ISP eq 11 then SP = 758.0d+0
IF ISP eq 12 then SP = 907.0d+0
IF ISP eq 13 then SP = 1096.0d+0
IF ISP eq 14 then SP = 1300.0d+0
IF ISP eq 15 then SP = 1500.0d+0

for IT1 = 1, 5 do begin    ;2) for T1
    for IT2 = 1, 3 do begin  ;3) for T2
      
      IF(IT1 EQ 1) then T1 = 135.0d+00
      IF(IT1 EQ 2) then T1 = 160.0d+00
      IF(IT1 EQ 3) then T1 = 213.0d+00
      IF(IT1 EQ 4) then T1 = 260.0d+00
      IF(IT1 EQ 5) then T1 = 285.0d+00
      
      R = 192.0d+00
      g = 3.72d+00
      H = R*T1/g ;scale height

      IF(IT2 EQ 1) then T2 = 80.0d+00
      IF(IT2 EQ 2) then T2 = 146.0d+00
      IF(IT2 EQ 3) then T2 = 200.0d+00
      
      kfike = path_LUT_Input + 'kfile/' + 'LUTable_T1_' + STRCOMPRESS(fix(T1), /REMOVE_AL) $
        + '_T2_'  + STRCOMPRESS(fix(T2), /REMOVE_AL) + '_PRS'  + STRCOMPRESS(fix(SP), /REMOVE_AL) + '.k'
      atmfile = path_LUT_Input + 'atmosfile/' + 'LUTable_T1_' + STRCOMPRESS(fix(T1), /REMOVE_AL) $
        + '_T2_'  + STRCOMPRESS(fix(T2), /REMOVE_AL) + '_PRS'  + STRCOMPRESS(fix(SP), /REMOVE_AL) + '.atmos'

      for ISZA = 1, 6 do begin ;4) SZA
        if ISZA eq 1 then SZA = 00.0
        if ISZA eq 2 then SZA = 15.0
        if ISZA eq 3 then SZA = 30.0
        if ISZA eq 4 then SZA = 45.0
        if ISZA eq 5 then SZA = 60.0
        if ISZA eq 6 then SZA = 75.0

        for IEA = 1, 3 do begin ;5) EA
          if IEA eq 1 then EA = 00.0
          if IEA eq 2 then EA = 05.0
          if IEA eq 3 then EA = 10.0

          for IPA = 1, 5 do begin ;6) Phase angle
            if IPA eq 1 then PA = 00.0
            if IPA eq 2 then PA = 45.0
            if IPA eq 3 then PA = 90.0
            if IPA eq 4 then PA = 135.0
            if IPA eq 5 then PA = 180.0

            for ID = 1, 6 do begin ;7) Dust
              if ID eq 1 then dust = 0.0
              if ID eq 2 then dust = 0.3
              if ID eq 3 then dust = 0.6
              if ID eq 4 then dust = 0.9
              if ID eq 5 then dust = 1.2
              if ID eq 6 then dust = 1.5

              for IWI = 1, 3 do begin ;8) Water ice 
                if IWI eq 1 then waterice = 0.0
                if IWI eq 2 then waterice = 0.5
                if IWI eq 3 then waterice = 1.0
                
                Azi_angle = phi_computing(SZA,EA,PA)
                
                for IAB = 1, 6 do begin ;9) Surface albedo
                  if IAB eq 1 then Surface_abledo = 0.05
                  if IAB eq 2 then Surface_abledo = 0.1
                  if IAB eq 3 then Surface_abledo = 0.2
                  if IAB eq 4 then Surface_abledo = 0.3
                  if IAB eq 5 then Surface_abledo = 0.4
                  if IAB eq 6 then Surface_abledo = 0.5

                  output_name_rad = path_LUT_Output + 'SP' + STRCOMPRESS(fix(ISP),/REMOVE_AL) $
                                                            + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                                            + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                                            + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                                            + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                                            + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                                            + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                                            + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                                            + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad.dat'
                                   
;                  output_name_rad_mono = path_LUT_Output + 'SP' + STRCOMPRESS(fix(ISP),/REMOVE_AL) $
;                                                            + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
;                                                            + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
;                                                            + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
;                                                            + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
;                                                            + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
;                                                            + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
;                                                            + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
;                                                            + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_rad_mono.dat'

                  output_name_trans = path_LUT_Output + 'SP' + STRCOMPRESS(fix(ISP),/REMOVE_AL) $
                                                            + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
                                                            + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
                                                            + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
                                                            + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
                                                            + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
                                                            + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
                                                            + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
                                                            + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_trans.dat'

;                  output_name_trans_mono = path_LUT_Output + 'SP' + STRCOMPRESS(fix(ISP),/REMOVE_AL) $
;                                                            + '_TA' + STRCOMPRESS(fix(IT1),/REMOVE_AL) $
;                                                            + '_TB' + STRCOMPRESS(fix(IT2),/REMOVE_AL) $
;                                                            + '_SZA' + STRCOMPRESS(fix(ISZA),/REMOVE_AL) $
;                                                            + '_EA' + STRCOMPRESS(fix(IEA),/REMOVE_AL) $
;                                                            + '_PA' + STRCOMPRESS(fix(IPA),/REMOVE_AL) $
;                                                            + '_Dust' + STRCOMPRESS(fix(ID),/REMOVE_AL) $
;                                                            + '_WaterI' + STRCOMPRESS(fix(IWI),/REMOVE_AL) $
;                                                            + '_SurfaceA' + STRCOMPRESS(fix(IAB),/REMOVE_AL) + '_trans_mono.dat'
;                                                            
                    ;------------------------------------------------------------------------------------------------------------------------
                    ;make ini file
                    ;------------------------------------------------------------------------------------------------------------------------
                    file = 'arsm_tmp_2.ini'
                    openw,lun,file,/get_lun
                    
                    printf,lun," 2 = .Info_Level."                                                                                                             
                    printf,lun,"' ' = .File_Info."
                    printf,lun,"'' = .File_Status. ! debug and multiple processes run control purpose"
                    printf,lun,"0 = .Intel_VML. ! 0 or 1"
                    printf,lun,""                                                                                                                              
                    printf,lun,"===   Output Files   ==="                                                                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun,"'"+output_name_trans+"'  = .File_Transmittance."                                                                                                  
                    printf,lun,"''  = .File_Monochromatic_Transmittance. "                                                                                    
                    printf,lun,"'"+output_name_rad+"'  = .File_Radiance."
                    printf,lun,"'' = .File_Monochromatic_Radiance."                                                                                          
                    printf,lun,""                                                                                                                              
                    printf,lun,"===   Wavenumber range and radiative transfer option   ================"                                                       
                    printf,lun,""                                                                                                                              
                    printf,lun," 4545 = .Vmin."                                                                                                                 
                    printf,lun," 5556  = .Vmax."           
                    printf,lun," 1     = .Geometry. ! 1 only"
                    printf,lun," 1     = .Planck_Source."                                                                                                       
                    printf,lun," 1     = .Solar_Source."                                                                                                        
                    printf,lun," 10000 = .V_Planck_Max."                                                                                                         
                    printf,lun," 1000 = .V_Solar_Min."                                                                                                          
                    printf,lun," 2 4 2 = .Radiative_Transfer_Options. ! see the manual, e.g.:"                                                                  
                    printf,lun," 16384 = .Computational_Subinterval_Size."                                                                                      
                    printf,lun," 1.1   = .Gas_Subinterval_Factor."                                                                                              
                    printf,lun,""                                                                                                                              
                    printf,lun,"===   Corr-I Correction   =============================================="                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun," 0     = .Cor_I_Partition.                  ! 10 is ok "                                                                        
                    printf,lun," 1000   = .Cor_I_Max_Subinterval_Size.       ! cm-1    "                                                                        
                    printf,lun," 0      = .Interpolate_Between_Subintervals.           "                                                                        
                    printf,lun,""                                                                                                                              
                    printf,lun,"=====   Model condition parameters   =================================="                                                       
                    printf,lun," 300.00   = .Observer_Altitude.             "
                    printf,lun," "+string(SZA,'(f5.2)')+" = .Theta_Sun."                                                                                                             
                    printf,lun," "+string(EA, '(f5.2)')+" = .Theta."                                                                                                             
                    printf,lun," 0.00  = .Phi_Sun."                                                                                                               
                    printf,lun," "+string(Azi_angle,'(f6.2)')+" = .Phi."                                                                                                             
                    printf,lun," 1.5 = .Heliocentric_Distance."                                                                                                             
                    printf,lun," 0.0 = .Surface_Temperature."                                                                                                 
                    printf,lun,""                                                                                                                              
                    printf,lun,""                                                                                                                              
                    printf,lun,"'/work1/data/ARS_example/pfsolspec_hr.dat' = .File_Solar_Radiance."     
                    printf,lun,"'"+atmfile+"' = .File_HTP."
                    printf,lun,"'HTP' = .HTP_Level_Selection."
                    printf,lun,""                                                                                                                              
                    printf,lun," "+string(Surface_abledo, '(f5.2)')+" = .Surface_Albedo."                                                                                                             
                    printf,lun,"' ' = .File_Surface_Albedo."                    
                    printf,lun,""                                                                                                                              
                    printf,lun,"=====   Monochromatic wavenumber grid   ================================"                                                     
                    printf,lun,""                                                                                                                              
                    printf,lun,"'/work1/LUT/SP/table/input/LUtable.v' = .File_Wavenumber_Grid."                    
                    printf,lun,""                                                                                                                              
                    printf,lun,"=====   Gases   ========================================================"                                                     
                    printf,lun,""                                                                                                                              
                    printf,lun,""                                                                                                                              
                    printf,lun," 1 = .Number_of_Species. "                                                                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun,""                                                                                                                              
                    printf,lun," 2 0 1                                  ! molecule and isotope codes, add.broad."
                    printf,lun,"' '          ! spectral lines file name"
                    printf,lun,"'/work1/LUT/SP/table/input/LUtable.v'                ! specific wavenumber grid file name"
                    printf,lun,"'"+kfike+"'          ! absorption coefficients file name"
                    printf,lun," 0.9532 ' ' 1. "
                    printf,lun," 1 0 0 0 0                                        ! line shape"                                                                 
                    printf,lun,""                                                                                                                              
                    printf,lun,""                                                                                                                              
                    printf,lun,"=====   Aerosols   ====================================================="                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun," 2 = .Number_of_Aerosols."          
                    printf,lun,""                                                                                            
                    printf,lun," 3 '/work1/LUT/Common/dust_2500_8500.aero'"                                                                                                                             
                    printf,lun," 2 5000.0 "+string(dust,format='(f6.4)')+" "+string(H*1e-3,format='(f5.2)')+" 0 60" 
                    printf,lun,""                                                                                            
                    printf,lun," 3 '/work1/LUT/Common/waterice_2500_8500.aero'"                                                                                                                             
                    printf,lun," 2 5000.0 "+string(waterice,format='(f6.4)')+" "+string(H*1e-3,format='(f5.2)')+" 0 60" 
                    printf,lun,""
                    printf,lun,"=====   Other settings   ==============================================="                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun," 1     = .Observation_Geometry. "                                                                                               
                    printf,lun," 1     = .Intel_VML."                                                                                                           
                    printf,lun,"' '   = .File_Info."                                                                                                           
                    printf,lun,"======   V-computation   ==============================================="                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun," 0.5d+0 = .Alpha."                                                                                                              
                    printf,lun," 1.5d+0 = .Beta."                                                                                                               
                    printf,lun," 100.d+0 = .Delta."                                                                                                             
                    printf,lun," 0.01d+0 = .Step_in_Wing."                                                                                                      
                    printf,lun," 1.e-5 = .Special_Grid_Taumin."                                                                                                 
                    printf,lun,""                                                                                                                              
                    printf,lun,"======   K-computation   ==============================================="                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun," 1 = .Line_Contour."                                                                                                            
                    printf,lun," 0 = .Normalization."                                                                                                           
                    printf,lun," 50 = .Line_Cutoff."                                                                                                            
                    printf,lun," 0 = .Taumin."                                                                                                                  
                    printf,lun," 1 = .Inform_Interval."                                                                                                         
                    printf,lun," 0 = .Grid_Centered."                                                                                                           
                    printf,lun," 0 = .HTP_Grid_Reduction. "                                                                                                     
                    printf,lun,""                                                                                                                              
                    printf,lun,"===   Limb sounding parameters (Unused)  ==============================="                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun," 0     = .Tangent_Sounding_Altitude.  [km]    ! for limb sounding only"                                                        
                    printf,lun," 43.4   = .Mean_Molecular_Weight.      [g//mole]! ---//---"                                                                      
                    printf,lun," 372.1  = .Gravitational_Acceleration. [cm//s2] ! ---//---"                                                                      
                    printf,lun," 3387.1 = .Radius_of_Planet.           [km]    ! ---//---"                                                                      
                    printf,lun," 0      = .Exponential_Tail.                   ! ---//---"                                                                      
                    printf,lun,""
                    printf,lun,"1 = .FBJ1J2_Warnings."
                    printf,lun,"0 = .Continuum_Absorption."
                    printf,lun,"' ' = .File_Continuum_Absorption."
                    printf,lun,""
                    printf,lun,"========================================================================"                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun," 3000 = .Sounding_Altitude."                                                                                                     
                    printf,lun," 0    = .Rayleigh_Scattering."                                                                                                  
                    printf,lun," 1    = .Delta_M."                                                                                                              
                    printf,lun," 32    = .Number_of_Streams.      "                                                                                              
                    printf,lun," 32   =  .Number_of_Moments.      "                                                                                              
                    printf,lun," 0    = .Set_Zero_Scattering."                                                                                                  
                    printf,lun,""                                                                                                                              
                    printf,lun,"===   Instrument specification   ======================================="                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun,"! 1-Rect, 2-Trian, 3-Gauss, 4-Sinc, 5-Sinc2, 6-Hamming, 7-Cos2, 8-User,File"                                                   
                    printf,lun,""                                                                                                                              
                    printf,lun," 3 = .Instrumental_Function."                                                                                                   
                    printf,lun,"'' = .File_Instrumental_Function."
                    printf,lun," 32.5D+0 = .Resolution.    [Ideal PFS: Sinc:1.211376319487D+0, Hamming:1.822245687D+0]"                                         
                    printf,lun,"                    [IRIS: Hamming => 2.1254448]"                                                                              
                    printf,lun,""                                                                                                                              
                    printf,lun," -1 = .Spectral_Step."                                                                                                          
                    printf,lun,"'/work1/LUT/SP/table/input/OMEGA_channel_cm.dat' = .File_Spectral_Channels.  ! if Spectral_Step < 0"                            
                    printf,lun,""                                                                                                                              
                    printf,lun,"===   Output format   =================================================="                                                      
                    printf,lun,""
                    printf,lun,"'e' = .Units_Output.                          ! for thermal region only"                                                       
                    printf,lun," 1 = .Two_Pass_Transmittance_Flag."                                                                                             
                    printf,lun,"'f12.7' 'f9.4' 'f10.6' '1p,e15.7,0p' 'f8.3' = .Output_Formats. [V,C,T,I,K]"                                                    
                    printf,lun,"   2   = .Transmittance_Path.                 "
                    printf,lun,""                                                                                                                              
                    printf,lun,"=== Unused parameters (from other programs) ============================"                                                      
                    printf,lun,""                                                                                                                              
                    printf,lun," 1048576 = .Buffer_Size."                                                                                                       
                    printf,lun," 1 = .Curtis_Godson_Modification."                                                                                             
                    printf,lun,""                                                                                                                              
                    printf,lun,"===  Other output  ====================================================="                                                      
                    printf,lun,""                                                                                                                             
                    printf,lun,""                                                                                                                              
                    printf,lun,""                                                                                                                              
                    printf,lun,"' '    = .File_Monochromatic_Reflectance. "                                                                                    
                    printf,lun,"' '    = .File_Reflectance. "
                    printf,lun,"' '    = .File_Monochromatic_Optical_Depth."
                    printf,lun,"' '    = .File_Monochromatic_Flux_Up."
                    printf,lun,"' '    = .File_Monochromatic_Flux_Down."
                    printf,lun,"' '    = .File_LibRadtran."
                    printf,lun,"' '    = .File_Medium_Properties."
                    printf,lun,""                                                                                                   
                    printf,lun," 0 = .Monochromatic_Dumps."
                    free_lun,lun
    
                    ;------------------------------------------------------------------------------------------------------------------------
                    ;run ARS
                    ;------------------------------------------------------------------------------------------------------------------------
                    spawn,'/work1/bin/arsm arsm_tmp_2.ini'


                endfor
              endfor
            endfor
          endfor
        endfor
      endfor
    endfor
endfor

stop
END