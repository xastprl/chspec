pro writespecfits

fileline='20210508_line_spec_LogTbin0.004_E0.1toE40.1keV_eBin0.001_SingIn_CHIANTIv10.0.1.idl'
;'20210131_line_spec_300Tgrid_SingIn_UnityAbund_with_CHIANTIv10.idl';'20200929_line_spec_300Tgrid_SingIn_UnityAbund_with_mod_ch_synthec.idl';'202005_line_spec_350Tgrid_eDens1e10_SingIn_UnityAbund.idl';'20200406_line_spec_300Tgrid_eDens1e11_SingIn_UnityAbund.idl';'20200217_line_spec_300Tgrid_SingIn_UnityAbund.idl'
filecont='20210618_cont_spec_LogTbin0.004_E0.1toE30.1keV_eBin0.001_SingIn_CHIANTIv10.0.1.idl'
;'20210131_cont_spec_300Tgrid_eDens1e10_UnityAbund_with_CHIANTIv10.idl';'20200214_cont_spec_300Tgrid_UnityAbund.idl';'20200512__cont_spec_350Tgrid_eDens1e10_UnityAbund.idl';'20200407_cont_spec_300Tgrid_eDens1e11_UnityAbund.idl';'20200214_cont_spec_300Tgrid_UnityAbund.idl'
restore,fileline
    ;Present line spectrum are for E = 0.1 to 40.0 with dE = 1eV, where as Continuum spectrum are for E = 0.1 to 30.1 with a dE = 1eV.
    ALLENE_ALL = ALLENE_ALL[0:29999]
    ALLSPEC_ALL = ALLSPEC_ALL[0:29999, *, *]
restore,filecont
ELEMENTS=ALLELEMENTS
ENERGIES=ALLENE_CONT
TEMPERATURE=TEMPERATURE_CONT[0:550];TEMPERATURE_ALL ;[0:550] -> 0.3 to 50 MK
SPECTRUM=allspec_all[*,0:550,*]+ALLSPEC_CONT[*,0:550,*]

;specfile='20200311_ContPlusLineSpec_300Tgrid_UnityAbund.idl'

;restore,specfile

outfile='IsotSpec_0p1to30p1KeV_0p001Ebin_0p32MKto50MK_eDens1e10_CHIANTIv10.fits';'IsotSpec_0p1to30p1KeV_0p001Ebin_0p32MKto99MK_eDens1e10_CHIANTIv10.fits'
;'IsotSpec_0p1to30p1KeV_0p01Ebin_1MKto20MK_eDens1e10_UsedModified_chsythetic.fits';'IsotSpec_0p1to30p1KeV_0p01Ebin_1MKto50MK_eDens1e10.fits';'IsotSpec_0p1to30p1KeV_0p01Ebin_1MKto20MK_eDens1e11.fits';'IsotSpec_0p1to30p1KeV_0.01Ebin_1MKto20MK.fits';'test.fits'

elememt=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']

;allspec_all[*,0,25]

nt=n_elements(temperature)
nelem=n_elements(spectrum[0,0,*])

a=0 
mwrfits,a,outfile

c={ELEMENTS_Z:(indgen(n_elements(ELEMENTS))+1),TEMP:TEMPERATURE,ENE:ENERGIES}
hdr=(["N_ELEMENTS = "+strtrim(string(n_elements(ELEMENTS)),1),"N_ENERGY = "+strtrim(string(n_elements(ENERGIES)),1),"N_TEMPERATURE = "+strtrim(string(n_elements(TEMPERATURE)),1)])
mwrfits,c,outfile,hdr

for el=0,nelem-1 do begin
    ;b={ELEMENT:elememt[el], SPEC:spectrum[*,50,25], TEMP:TEMPERATURE[50]}
    b={SPEC:spectrum[*,50,25]}
    spec=REPLICATE(b, nt)
    for i=0,nt-1 do begin
	;spec[i]={ELEMENT:elememt[el], SPEC:spectrum[*,i,el], TEMP:TEMPERATURE[i]}
	spec[i]={SPEC:spectrum[*,i,el]}
    endfor
    mwrfits,spec,outfile
endfor

stop
end

