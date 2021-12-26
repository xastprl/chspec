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
;; 14/02/2020 : Coppied from var_abund_my_get_flarespec.pro and abund_var_isothermal_ke.pro. to calculate only continuum intensity of indivudul elements and save them in a array
;; Example::
pro isothermal_continuum_spec
TIC
	;Chianti init
        ;!PATH = expand_path('+/data/users/biswajit/software/chianti/CHIANTI_9.2.1_pro_standalone/')+':'+ !PATH
        ;use_chianti, '/data/users/biswajit/software/chianti/CHIANTI_9.0.1_data'
        !PATH = expand_path('+/data/users/biswajit/software/chianti/CHIANTI_10.1.1_pro_standalone/')+':'+ !PATH
        use_chianti, '/data/users/biswajit/software/chianti/CHIANTI_10.0.1_database'

	;Energy in keV
        emin=0.1
        emax=30.1;40.1;20.495
        ebin=0.001;0.005
	
	nebins=long64((emax-emin)/ebin)+1
	;nw=fix((emax-emin)/ebin)+1
	;binsize=findgen(nebins)+ebin
	array_energies=findgen(nebins)*ebin+emin 

	edensity=1d10;1d11;1d9;1d10	
        ;abund_files='/data/users/biswajit/software/chianti/CHIANTI_9.0.1_data/abundance/biswajit_abund/'+new_abund_elem_name+'/*.abund';sun_coronal_2012_schmelz_ext.abund';biswajit_abund/*.abund'
	;abund_files='/data/users/biswajit/Desktop/projects/phd/xsm/xspec_tablemodel/'+new_abund_elem_name+'/*.abund'
	abund_files='/data/users/biswajit/software/chianti/CHIANTI_10.0.1_database/abundance/unity.abund'

	sundist=1.495978707 ; in 1e13 cm
	sundist2=sundist*sundist ; in 1e26 cm^2
	erg2watt=1.e-3
	keV2J = 1.6e-16

	na=n_elements(file_search(abund_files))
        allene_cont=fltarr(nebins)
	vol_em=1.0e20
	Start=5.5
	Finish=8.0
        increment=0.004
        n_temperature=(Finish-Start)/increment
        temperature_cont = 10^findgen(n_temperature,start=Start,increment=0.004)

        ;Selection of elements

        allelements=['h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe','co','ni','cu','zn']
        masterlist_file='/data/users/biswajit/software/chianti/CHIANTI_10.0.1_database/masterlist/masterlist.ions'
        read_masterlist,masterlist_file,all_ions

        n_elem=n_elements(allelements)
	allspec_cont=dblarr(nebins,n_elements(temperature_cont),n_elem)

;for i1=0,(n_temperature-1) do begin
;for i1=0,0 do begin

	temperature = temperature_cont


    	for j=0,n_elem-1 do begin

	ion_list=all_ions[where(STRMATCH(all_ions, allelements[j]+'_'+'*') eq 1)]

		;Compute Continum spectrum from Chianti for each elements
		;setup_elements,ioneq_file=!ioneq_file, abund_file=abund_files ; Requie to get the values(i.e, abund) in common block to call freebound.pro with /no_setup option for CHIANTIv9.01

		;freebound, temperature, array_energies,fb,/no_setup,min_abund=min_abund, /photons, em_int=em_int,/kev,iz=(j+1)
                freebound, temperature, array_energies,fb,ABUND_FILE=abund_files,min_abund=min_abund, /photons, em_int=em_int,/kev,Element=allelements[j],ioneq_file=!ioneq_file

		;two_photon, 1.0e7, array_energies, two_phot,abund_file=abund_file,min_abund=min_abund,edensity=1.0e10,sngl_ion=ion_list, photons=photons, em_int=em_int,/kev
		;two_photon, temperature, array_energies, two_phot,abund_file=abund_files,min_abund=min_abund,edensity=edensity,ioneq_file=!ioneq_file,sngl_ion=ion_list, /photons, em_int=em_int,/kev
                two_photon, temperature, array_energies, two_phot,abund_file=abund_files,min_abund=min_abund,edensity=edensity,ioneq_file=!ioneq_file,Element=allelements[j], /photons, em_int=em_int,/kev
		
		freefree,temperature, array_energies,ff,abund_file=abund_files,min_abund=min_abund,element=allelements[j],ioneq_file=!ioneq_file, /photons, em_int=em_int,/kev
		
		totalcount=(two_phot+ff+fb)/1d40*ebin
		allspec_cont[*,*,j]=totalcount* (vol_em/sundist2)
		allene_cont=array_energies	; in units of KeV

	endfor

;endfor


;spec=0
;for i=0,n_elements(allspec_cont[0,0,*])-1 do begin
;spec=spec+ALLSPEC_CONT[*,0,i]
;end


	save,abund_files,allelements_cont,temperature_cont,allene_cont,allspec_cont,filename='20210618_cont_spec_LogTbin0.004_E0.1toE30.1keV_eBin0.001_SingIn_CHIANTIv10.0.1.idl'
print,"========================= END CALCULATION ============================"
TOC
stop

end

;!PATH = expand_path('+/data/users/biswajit/software/chianti/CHIANTI_10.1.1_pro_standalone/')+':'+ !PATH
;use_chianti, '/data/users/biswajit/software/chianti/CHIANTI_10.0.1_database'
;abund_files='/data/users/biswajit/software/chianti/CHIANTI_10.0.1_database/abundance/sun_coronal_2012_schmelz_ext.abund'                                    
;ebin=0.01
;array_energies=findgen(40000)*ebin+0.1
;temperature=80.0e6
;freebound, temperature, array_energies,fb,ABUND_FILE=abund_files,min_abund=min_abund, /photons, em_int=em_int,/kev,ioneq_file=!ioneq_file
;two_photon, temperature, array_energies, two_phot,abund_file=abund_files,min_abund=min_abund,edensity=edensity,ioneq_file=!ioneq_file, /photons, em_int=em_int,/kev
;freefree,temperature, array_energies,ff,abund_file=abund_files,min_abund=min_abund,ioneq_file=!ioneq_file, /photons, em_int=em_int,/kev
;
;totalcount=(two_phot+ff+fb)/1d40*ebin
;plot,array_energies,totalcount,/yl,/xl
