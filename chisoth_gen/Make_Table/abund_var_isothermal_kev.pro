
;+
; NAME
;
;       ABUND_VAR_ISOTHERMAL_KEV
;
; PURPOSE:
;       Computes CHIANTI isothermal spectra 
;       Note that this routine has a number of unique features that 
;       distinguish it from the other CHIANTI synthetic spectra routines. 
;       See the Programming Notes section.
;
; INPUTS:
;
;       MIN      Minimum of desired energy range in keV
;
;       MAX      Maximum of desired energy  range in keV
;
;       STEP      Bin size of spectrum (in keV)
;
;       TEMP      Electron temperature (or array) in K.
;
; OPTIONAL INPUTS
;
;       PRESSURE  Electron pressure in units of cm^-3 K.
;
;       EDENSITY  Electron density in units of cm^-3.
;
;       EM        Emission measure. The units of EM govern the intensity 
;                 units of the emission lines (i.e., column or volume 
;                 emission measure units). If EM is not specified, then the 
;                 emission measure is set to (N_e * N_H) where N_e is 
;                 derived from the user-specified PRESSURE or EDENSITY, 
;                 and N_H is derived from the routine PROTON_DENS.PRO.
;
;       SNGL_ION  Rather than include the entire list of CHIANTI ions in 
;                 the calculation, this input can be used to select a 
;                 single ion, or a number of different ions. E.g., 
;                 SNGL_ION='s_2' or SNGL_ION=['s_2','s_3','s_4'].
;
;       RADTEMP   The blackbody radiation field temperature (default 6000 K).
;
;       RPHOT    Distance from the centre of the star in stellar radius units.
;                I.e., RPHOT=1 corresponds to the star's surface. (Default is
;                infinity, i.e., no photoexcitation.)
;
;       MASTERLIST  The list of ions that will be considered for the 
;                   spectrum is contained in the masterlist file in the 
;                   CHIANTI directories. The user can specify his own file 
;                   through this keyword. E.g., 
;                   masterlist='/user/home/masterlist.ions'
;
;
;	ABUND_FILES:  Name of the abundance file/files(directory of a file 
;		      where .abund files exist) to use.  If not passed, then
;		     the user is prompted for it.
;
;	IONEQ_NAME:  Name of the ionization equilization name to use.  If not
;		     passed, then the user is prompted for it.
;
;
;  KEYWORDS
;
;       NOPROT     Switch off the inclusion of proton rates in the level 
;                  balance.
;
;       ERGS       The units of the output spectrum are by default in photons. 
;                  Setting /ERGS switches to erg units.
;
;       CONT       Adds continuum (free-free, free-bound, two-photon) to 
;                  spectrum.
;  
;       ALL        Include all lines, i.e. also those for which wavelengths 
;                  are only theoretical and not observed. 
;
;  OUTPUTS:
;
;        array_energies   array of energies for the calculated synthetic spectrum.
;
;        SPECTRUM Intensity array. The units depend on the user inputs to 
;                 ISOTHERMAL -- see note below. 
;
;        LIST_ENERGIES  A list of ENERGIES of the spectral lines for
;        use with synthetic_plot.pro 
;
;        LIST_IDENT A list of line identifications for use with 
;                   synthetic_plot.pro
;
; PROGRAMMING NOTES
;
;        Intensity Units
;        ---------------
;        The units by default are of the form photons cm^3 s^-1 sr^-1 * (units of EM), 
;        changing to ergs if the /ergs keyword is set.
;
;        The volume emission measure (N_e*N_H*V) has units cm^-3.        
;
;
;        The energy array is a linear one with STEP binwidth.
;
;        A thermal width is added to each spectral line. If that is
;        larger than the bin size the line intensity is spread around
;        pixels with a Gaussian.
;
;        The spectrum, unlike the result of ch_ss, is not per Angstrom
;        or per keV.
; 
;        to get an observed flux in  ergs (or photons) cm-2 s-1 you need to multiply the
;        spectrum for the volume EM (say 4d47), divide by the solar
;        distance^2, and integrate
;
;	 To get the spectra for more than one elemental abundance file, put the directory  
;	 of all .abund file as a input (like, abund_files='~/directory/*.abund' )
;
;	 If more than one temperature and abundance file is given, then the output 
;	 allspectrum will be a 3D array of [Tempeeraure_index,intensity,abundance_file_index]
;
; EXAMPLE
;
;  ; GOES 1-8 Angstroms: 
;  abund_var_isothermal_kev,1.5,12.4, 0.017, [1e7,3e7],array_energies,spectrum,list_en,list_ident,edens=1e10,/cont,/ergs
; 
;  sngl_ion='fe_23'
;
;
;  Sun distance suared:  (14959787070000.d*14959787070000.d) = 2.2379523e+26
;   a typical GOES flare of 10 MK has a volume EM of 4d47
;   
;   from ergs to Watts need to multiply by 1.e-3.
;
;   print, 1.e-3* 4.d21/2.2379523 *    total(spectrum[*,0]) = 8e-07
;   Watts/cm-2 which is an X-ray flux of a B8-class flare.
;
;
; COMMON BLOCKS
;
;        ELEMENTS
;
; CALLS
;
;        CH_SYNTHETIC, READ_ABUND, CH_GET_FILE, CONCAT_DIR, FREEFREE, 
;        FREEBOUND, TWO_PHOTON
;
; MODIFICATION HISTORY
;
;
;       v.1, Giulio Del Zanna (GDZ),  5 Jul 2016 
;             copied from isothermal and make_chianti_spec but works in kev
;	
;	V.2.0, Biswajit, 22 Oct 2019
;		Make isothermal_kev to abund_var_isothermal_kev for the calculation of
;		isothermal spectra for different elemental abundance file at a same time
;		by calling CH_ISOTHERMAL only one time in the calculation.
;		Updated for new version- of Chianti-9.0.1
;
;       V.2.1, Biswajit, 12 Feb 2020   
;		include 'allspectrum=0' in line-238 to set allspectrum=0 for no existing 
;		line in the  ch_synthetic calculation 
;
; VERSION     : 
;       Version 1. TEST version
;   	Version 2.0
;	Version 2.1
;  
;-


PRO abund_var_isothermal_kev, min, max, step,temp, array_energies,allspectrum,list_energies,list_ident,$
      pressure=pressure,edensity=edensity,photons=photons, ergs=ergs, $
      sngl_ion=sngl_ion, abund_files= abund_files , ioneq_name=ioneq_name, $
      noverbose=noverbose, min_abund=min_abund,  $
      masterlist=masterlist, noprot=noprot, radtemp=radtemp, $
      rphot=rphot, em=em,all=all, cont=cont
TIC
COMMON elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

IF n_params() LT 6 THEN BEGIN
  print,'Use: IDL> isothermal_kev, min, max, step,temp,lambda,allspectrum,list_wvl,'
  print,'                 list_ident, pressure= , edensity= , /photons,'
  print,'                 /noverbose, sngl_ion= , min_abund= , /cont,'
  print,'                 masterlist= ,abund_files= , ioneq_name= ,/noprot]'
  print,'                 radtemp= ,rphot= ,em= , /all '
  print,''
  print,'E.g., IDL> isothermal,170,180,0.1,1e6,lambda,allspectrum,list_wvl,'+$
       'list_ident,'
  print,'                  edensity=1e10'
  return
ENDIF

IF keyword_set(ergs) THEN photons=0 ELSE photons=1

;IF n_elements(spectrum) NE 0 THEN junk=temporary(spectrum)

IF n_elements(allspectrum) NE 0 THEN junk=temporary(allspectrum)

IF n_elements(min_abund) EQ 0 THEN min_abund=0.

IF (n_elements(pressure) NE 0 AND n_elements(edensity) NE 0) OR $
   (n_elements(pressure) EQ 0 AND n_elements(edensity) EQ 0) THEN BEGIN
  print,'%ISOTHERMAL: Please specify either EDENSITY or PRESSURE'
  return
ENDIF

IF keyword_set(noverbose) THEN verbose=0 ELSE verbose=1

n=n_elements(edensity)
nt=n_elements(temp)

IF n_elements(em) NE 0 THEN BEGIN
  IF n_elements(em) EQ nt THEN BEGIN
    logem_isothermal=alog10(em)
  ENDIF ELSE BEGIN
    logem_isothermal=alog10(em[0])
  ENDELSE
ENDIF ELSE BEGIN
  logem_isothermal=dblarr(nt)
ENDELSE

IF n_elements(ioneq_name) EQ 0 THEN BEGIN 
  dir=concat_dir(!xuvtop,'ioneq')
  ioneq_name=ch_get_file(path=dir,filter='*.ioneq', $
                         title='Select Ionization Equilibrium File')
IF NOT file_exist(ioneq_name) THEN message, 'Error, no file selected -- EXIT'
END

IF n_elements(abund_files) EQ 0 THEN BEGIN 
  dir=concat_dir(!xuvtop,'abundance')
  abund_files=ch_get_file(path=dir,filter='*.abund', $
                         title='Select Abundance File')
IF NOT file_exist(abund_files) THEN message, 'Error, no file selected -- EXIT'
END

IF nt NE 1 THEN no_sum_int=1

spectrum=-1
allspectrum=0
ang2kev=12.39854
wmin=ang2kev/max
wmax=ang2kev/min



ch_synthetic,wmin,wmax,out=transitions, press=pressure, err_msg=err_msg, $
     sngl_ion=sngl_ion,ioneq_name=ioneq_name,dem_name=dem_name, $
     photons=photons, masterlist=masterlist, noprot=noprot, $
     radtemp=radtemp, rphot=rphot, verbose=verbose, progress=progress, $
     density=edensity, no_sum_int=no_sum_int, logt_isothermal=alog10(temp), $
     logem_isothermal=logem_isothermal,all=all

IF err_msg [0]  NE  '' THEN BEGIN 
  print, '% ISOTHERMAL: Error - EXIT'
  return
END 

nl=n_elements(transitions.lines)

format='('+strtrim(string(nt),2)+'e10.2)'

; include abundances

abund_files_name=file_search(abund_files)

;nw=fix((max-min)/step)+1
nw=LONG64((max-min)/step)+1
;;array_energies=findgen(nw)*step+min
na=n_elements(abund_files_name)
spectrum=dblarr(nw,nt)
dim_spectrum=size(spectrum)
if dim_spectrum[0] eq 1 then  allspectrum=dblarr(nw,na)
if dim_spectrum[0] gt 1 then  allspectrum=dblarr(nw,nt,na)
;;binsize=findgen(nw)+step
;;format='('+strtrim(string(nt),2)+'e10.2)'
;;str1=transitions.lines.snote+' '+transitions.lines.ident
;;cc = SQRT(2.*!pi)/2./SQRT(2.*ALOG(2))
;;wvl = transitions.lines[*].wvl
;;line_energies = ang2kev/wvl


for abund_ind=0,n_elements(abund_files_name)-1 do begin

abund_name=abund_files_name[abund_ind]

read_abund,abund_name,abund,ref
;read_abund,abund_files_name[abund_ind],abund,ref


line_abunds = abund[transitions.lines.iz-1]
wvl = transitions.lines[*].wvl
; intensities = transitions.lines[*].int * line_abunds[*] 
line_energies = ang2kev/wvl


; setup a linear array 
;nw=fix((max-min)/step)+1
nw=LONG64((max-min)/step)+1
binsize=findgen(nw)+step

array_energies=findgen(nw)*step+min
spectrum=dblarr(nw,nt)

;;format='('+strtrim(string(nt),2)+'e10.2)'

str1=transitions.lines.snote+' '+transitions.lines.ident


 cc = SQRT(2.*!pi)/2./SQRT(2.*ALOG(2))

FOR i=0L, long(nl-1) DO BEGIN
   
; fwhm in Angstroms:
   fwhm = SQRT( 5.093d-13 * 10.^(transitions.lines[i].tmax) * $ 
                transitions.lines[i].wvl^2 / get_atomic_weight(transitions.lines[i].iz) )

; ?? FWHM in keV 
   fwhm_energy=  ang2kev/(transitions.lines[i].wvl-(fwhm/2.))-ang2kev/(transitions.lines[i].wvl+(fwhm/2.))
   
   int=transitions.lines[i].int * line_abunds[i]
   str1[i]=strpad(str1[i],50,/after)+'  Int='+string(format=format,int)
   
;; The if statements below check if the line extends over more than one 
;; pixel. If it doesn't, then the way the line is added to the spectrum is 
;; different.

  ind = WHERE( ABS(array_energies - line_energies[i]) LT 2.5*fwhm_energy )
   
        IF (ind[0] EQ -1) OR (N_ELEMENTS(ind) EQ 1) THEN BEGIN

            getmin = MIN(ABS(array_energies - line_energies[i]),ind2)
            spectrum[ind2,*] = spectrum[ind2,*] +  $
                              (transitions.lines[i].int * line_abunds[i])
     
            
                     ENDIF ELSE BEGIN
                        
                        
            int_lambda = [array_energies[ind[0]] - (binsize[ind[0]])/2., array_energies[ind]+binsize[ind]/2.]
            int_lambda = (int_lambda- line_energies[i])* SQRT(8.*ALOG(2)) / fwhm
            
            nind = N_ELEMENTS(ind)+1
            ints = DBLARR(nind, nt)
            amplitude = (transitions.lines[i].int * line_abunds[i])/fwhm/cc
            
            FOR j=0,nind-1 DO  $
               ints[j, *] = amplitude*cc*fwhm*GAUSSINT(int_lambda[j])
            
            line_profile = ints[1:nind-1, *]-ints[0:nind-2, *]
            spectrum[ind, *] = spectrum[ind, *] + line_profile 
            
         ENDELSE  
                                          
ENDFOR
         
list_ident=str1

i=sort(line_energies) 
list_energies=line_energies[i]
list_ident=list_ident[i]



if keyword_set(cont) then begin 
   
   print, '% ISOTHERMAL: calculating continuum..'
   
 ;
  em_int=double(10.)^logem_isothermal
 ;
  IF n_elements(edensity) EQ 0 THEN edensity=pressure/temp
  freebound, temp, array_energies,fb,/no_setup,min_abund=min_abund, $
             photons=photons, em_int=em_int,/kev
  
  ;freefree,temp, array_energies,ff,/no_setup,min_abund=min_abund, $
  ;     photons=photons, em_int=em_int,/kev
  
  freefree,temp, array_energies,ff,abund_file=abund_name,min_abund=min_abund, $
        photons=photons, em_int=em_int,/kev

  ;two_photon, temp, array_energies, two_phot,/no_setup,min_abund=min_abund, $
  ;     edensity=edensity, photons=photons, em_int=em_int,/kev

  two_photon, temp, array_energies, two_phot,abund_file=abund_name,min_abund=min_abund, $
       edensity=edensity, photons=photons, em_int=em_int,/kev
  
; the continuum is in  units /keV so to get it in units per bin we
; have to multiply for the step. 
  
  totcont=(fb+ff+two_phot)/1d40*step
  spectrum=spectrum+totcont
  
endif


print, '% ISOTHERMAL: calculation finished'


if dim_spectrum[0] eq 1 then  allspectrum[*,abund_ind]=spectrum
if dim_spectrum[0] gt 1 then  allspectrum[*,*,abund_ind]=spectrum


endfor
TOC
END
