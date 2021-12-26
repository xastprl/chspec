;; This routine computes synthetic solar flare spectrum for different classes based 
;; on temperature and emission measure from ***REF****
;;
;; It calls abund_var_isothermal_kev module from Chianti. 
;; This is for creating table model for the xspec fitting of XSM data.
;;
;; Biswajit
;; 19/09/2019
;; Modified- 07/11/2019
;; 09/11/2019 : Corrected some bugs
;; 11/11/2019 : Make it use for parallel run in different idl terminal simultaneously for less computing and add a parameter sl_no (serial number)
;; 11/11/2019 : Same as var_abund_my_get_flarespec.pro, but use for parallel run in different idl terminal simultaneously for less computing and add a parameter sl_no (serial number) and new_abund_elem_name
;; 12/02/2020 : Coppied from var_abund_my_get_flarespec.pro. to calculate only line intensity of indivudul elements and save them in a array
;; 17/02/2020 : Coppied from isothermal_line_spec.pro. to give a single temperature inpu to ch_synthetic instead of a temperature array.
pro isothermal_line_spec_singleTin
TIC
	;Chianti init

        !PATH = expand_path('+/data/users/biswajit/software/chianti/CHIANTI_10.1.1_pro_standalone/')+':'+ !PATH
        use_chianti, '/data/users/biswajit/software/chianti/CHIANTI_10.0.1_database'	

        ;Energy in keV
        emin=0.1
        emax=40.1;20.495
        ebin=0.001;0.005
	
	nebins=long((emax-emin)/ebin)+1

	edensity=1d10 ;1d11;1d9;1d10	
	abund_files='/data/users/biswajit/software/chianti/CHIANTI_10.0.1_database/abundance/unity.abund'

	sundist=1.495978707 ; in 1e13 cm
	sundist2=sundist*sundist ; in 1e26 cm^2
	erg2watt=1.e-3
	keV2J = 1.6e-16

	na=n_elements(file_search(abund_files))
        allene_all=fltarr(nebins)

	vol_em=1.0e20;3.0e21
	Start=5.5;1.0d6 ;6.0 ;1.0d6
	Finish=8.0;2.0d7 ;7.75;2.0d7
        increment=0.004
        n_temperature=(Finish-Start)/increment
        alltemperature = 10^findgen(n_temperature,start=Start,increment=0.004)

	;alltemperature = cgLogGen(nfl_all, Start=Start, Finish=Finish)

        ;Selection of elements

        allelements=['h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe','co','ni','cu','zn']
	elements_Z=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
	
        masterlist_file='/data/users/biswajit/software/chianti/CHIANTI_10.0.1_database/masterlist/masterlist.ions';'/data/users/biswajit/software/chianti/CHIANTI_9.0.1_data/masterlist/masterlist.ions'
        read_masterlist,masterlist_file,all_ions

        n_elem=n_elements(allelements)
	allspec_all=dblarr(nebins,n_temperature,n_elem)

;for i1=0,(n_iter-1) do begin
for i1=0,(n_elements(alltemperature)-1) do begin
	temperature = alltemperature[i1]

	allene=fltarr(nebins)

    	for j=0,n_elem-1 do begin

	ion_list=all_ions[where(STRMATCH(all_ions, allelements[j]+'_'+'*') eq 1)]

		;; Run Chianti for each flare class parameters
		;Compute synthetic spectrum from Chianti

		;abund_var_isothermal_kev,emin,emax,ebin,temperature,allene,allspec,list_ene,list_iden,edensity=edensity,/cont,ioneq_name=!ioneq_file,abund_files= abund_files,/photons,/noverbose;,sngl_ion=['o_8'];ion_list;['o_8'];
		abund_var_isothermal_kev,emin,emax,ebin,temperature,allene,allspec,list_ene,list_iden,edensity=edensity,ioneq_name=!ioneq_file,abund_files= abund_files,/noverbose,/photons,sngl_ion=ion_list;['o_8'];
                ;abund_var_isothermal_kev_v10,emin,emax,ebin,temperature,allene,allspec,list_ene,list_iden,edensity=edensity,ioneq_name=!ioneq_file,abund_files= abund_files,/noverbose,/photons,sngl_ion=ion_list

		allspec_all[*,i1,j]=allspec* (vol_em/sundist2)
	endfor

endfor
	allene_all=allene       ; in units of KeV
	temperature_all=alltemperature
	save,abund_files,allelements,temperature_all,allene_all,allspec_all,filename='20210508_line_spec_LogTbin0.004_E0.1toE40.1keV_eBin0.001_SingIn_CHIANTIv10.0.1.idl'
        
print,"========================= END CALCULATION ============================"
TOC
stop

end
