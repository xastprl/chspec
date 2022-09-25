/*
CHIANTI Model for XSPEC to fit the X-ray Spectrum.
It Uses the CHIANTI database version-10.0.1.

Notes: 
      This version of the model can be used for an instrument 
      response of lower energy bound 0.1 keV or above and a 
      bin size of an integer multiple of 1eV.
      Calculated spectrum are valid in the energy rage of 0.1-30.0 keV.
      The present version of the model is not efficient for the model 
      calculation beyond the mentioned energy range. 
      Above 30 keV it sets all the counts to zero.

Last Modified:
      Oct,2020 :   Biswajit, Santosh, Mithun.
      June,2021:   Biswajit, (Generalized for the use of other instruments.)
*/

#include <functionMap.h>
#include <xsTypes.h>
#include <FunctionUtility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <iomanip>
#include <bits/stdc++.h>
#include "fitsio.h"
#include <cstdio>

#define I_n 30
#define T_n 551
#define E_n 30000

struct filehandle
{
    // Variable definitions

    double*** allspec;
    double temp[T_n];
    double ene[E_n];
    
    filehandle()
    {
        fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
        int hdunum,status, hdutype;
	    long frow, felem, nelem;
	    double spec[E_n];
        
        allspec=(double***)malloc(I_n * sizeof(double**));    
	    if (allspec == NULL){
                fprintf(stderr, "Out of memory");
                exit(0);
	    }
	    for (int i = 0; i < I_n; i++)
        {
            allspec[i] = (double**)malloc(T_n * sizeof(double*));
            for (int j = 0; j < T_n; j++)
                allspec[i][j] = (double*)malloc(E_n * sizeof(double));
		}

	    status = 0; 
	    char filename[]  = "";
	    fits_open_file(&fptr, filename, READONLY, &status); 
            frow      = 1;
            felem     = 1;
            nelem     = T_n;
	    fits_movabs_hdu(fptr, 2, &hdutype, &status);
	    fits_read_col(fptr, TDOUBLE, 2, frow, felem, nelem, NULL, temp,NULL, &status);
	    fits_read_col(fptr, TDOUBLE, 3, frow, felem, E_n, NULL, ene,NULL, &status);

	    for (hdunum = 3; hdunum < 33; hdunum++)
        {
            /* move to the HDU */
            fits_movabs_hdu(fptr, hdunum, &hdutype, &status);
            for (int j = 0; j < T_n; j++)
            {
                fits_read_col(fptr, TDOUBLE, 1, j+1, felem, E_n, NULL,  spec, NULL, &status);
                for (int k = 0; k < E_n; k++) { allspec[hdunum-3][j][k] = spec[k];}
            } 
        }

        fits_close_file(fptr, &status);

	    //printf("ChIsoth: Data Read!\n");

    }

    ~filehandle() 
    {
        for (int i = 0; i < I_n; i++)
        {
            for (int j = 0; j < T_n; j++) free(allspec[i][j]);
            free(allspec[i]);
	    }
        free(allspec);
    }

};


filehandle chfh;

extern "C"
void chisoth(const RealArray& energyArray, const RealArray& parms, 
               int spectrum, RealArray& flux, RealArray& fluxErr, 
               const string& init)

{
    double T,T1, mf, ad, tL, tR, eL, eR, dydx,dTRL,dTL ;
    float abund[I_n] ;
    int nE, en, inL, inR;
    
    T1=parms[0] ;// Log(T)
    T=pow(10,T1); // in K
    nE = energyArray.size()-1;
    flux.resize(nE);
    for (int p=1; p<I_n; p++) {abund[p]= parms[p];}
    abund[0]=12.0;
 
    inL = 0;
    double dlt=0.004;
	double logTMin=6.0;

    if (T >= chfh.temp[T_n - 2]) inL = T_n - 2;
	else if (T<=chfh.temp[1]) inL=0;
    else{
        inL=(int)((T1-logTMin)/dlt)-1;
        while ( T > chfh.temp[inL + 1] ) inL++;
    }

	inR=inL+1;

    tL = chfh.temp[inL];
    tR = chfh.temp[inR];
    
    // For interpolation in linear T
    //dTRL=tR - tL;
    //dTL=T - tL;

    // For interpolation in log T
    dTRL=log10(tR/tL);
    dTL=log10(T/tL);

    
    int LowerBound_ind = 0;
    int Ebound_error = 0;
    for (LowerBound_ind = 0; LowerBound_ind < E_n;LowerBound_ind++)
    {
             if (abs(float(chfh.ene[LowerBound_ind]) - float(energyArray[0])) < 1.0e-6) 
	     {   
		 Ebound_error = 0;
		 break;
	     }
	     else Ebound_error = 1;
    }
    if (int(Ebound_error) == 1) 
    {
	printf("ChIsoth Error: Lower bound of the energy response should be >= 100 eV or should be incremented as an integer multiple of 1 eV. ...Terminating the program.\n");
	exit(0);
    } 
    for (en = 0; en < nE; en++) flux[en]=0.0;

    for (int i = 0; i < I_n; i++)
	{
	    ad = abund[i] - abund[0];
	    mf = pow(10, ad);
	int ene_ind = LowerBound_ind;
	for (en = 0; en < nE; en++)
        {       eR = 0.0;
		eL = 0.0;
		while(chfh.ene[ene_ind] >= energyArray[en] and chfh.ene[ene_ind] < energyArray[en+1]) 
		{ 
                    eR = eR + chfh.allspec[i][inR][ene_ind];
                    eL = eL + chfh.allspec[i][inL][ene_ind];
		    ene_ind = ene_ind+1;
		}
            dydx= (eR-eL)/dTRL;
            flux[en]+= mf*(eL + dydx*dTL);		
        }
    }

    //printf("ChIsoth: Finished iteration\n");
}


/*
ln -sf lmodel_isothermal.dat lmodel.dat
echo "initpackage isoth lmodel_isothermal.dat `pwd` \nexit" | xspec
echo "lmod isoth . \nexit" | xspec
energies 0.1 30.1 3001
ulimit -s 500000
*/
