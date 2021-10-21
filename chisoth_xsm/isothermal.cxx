/*
CHIANTI Model for XSPEC to fit the XSM Spectrum.
It Uses the CHIANTI database version-9.01.

Biswajit, Santosh, Mithun .
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
#define T_n 300
#define E_n 3001

struct filehandle
{

    double*** allspec;
    double temp[T_n];
    double ene[E_n];
    
    filehandle()
    {
        fitsfile *fptr;      
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
	    char filename[]  = "/home/biswajit/Desktop/projects/phd/xsm/IsothModel/IsotSpec_0p1to30p1KeV_0p01Ebin_1MKto20MK_eDens1e10_UsedModified_chsythetic.fits";
	    fits_open_file(&fptr, filename, READONLY, &status); 
            frow      = 1;
            felem     = 1;
            nelem     = T_n;
	    fits_movabs_hdu(fptr, 2, &hdutype, &status);
	    fits_read_col(fptr, TDOUBLE, 2, frow, felem, nelem, NULL, temp,NULL, &status);
	    fits_read_col(fptr, TDOUBLE, 3, frow, felem, E_n, NULL, ene,NULL, &status);

	    for (hdunum = 3; hdunum < 33; hdunum++)
        {
            
            fits_movabs_hdu(fptr, hdunum, &hdutype, &status);
            for (int j = 0; j < T_n; j++)
            {
                fits_read_col(fptr, TDOUBLE, 1, j+1, felem, E_n, NULL,  spec, NULL, &status);
                for (int k = 0; k < E_n; k++) { allspec[hdunum-3][j][k] = spec[k];}
            } 
        }

        fits_close_file(fptr, &status);


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
    
    T1=parms[0] ;
    T=pow(10,T1); 
    nE = energyArray.size()-1;
    flux.resize(nE);

    for (int p=1; p<I_n; p++) {abund[p]= parms[p];}
    abund[0]=12.0;
 
    inL = 0;
    double dlt=0.004336767364;
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
    
    //dTRL=tR - tL;
    //dTL=T - tL;

    dTRL=log10(tR/tL);
    dTL=log10(T/tL);

    for (en = 0; en < nE; en++) flux[en]=0.0;

    for (int i = 0; i < I_n; i++)
	{
	    ad = abund[i] - abund[0];
	    mf = pow(10, ad);
		for (en = 0; en < nE; en++)
        {
            eR=chfh.allspec[i][inR][en];
            eL=chfh.allspec[i][inL][en];
            dydx= (eR-eL)/dTRL;
            flux[en]+= mf*(eL + dydx*dTL);		
        }
    }

}
