
#ifndef READTABLE_H_
#define READTABLE_H_

#include <xsTypes.h>
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
//#include <bits/stdc++.h>
#include "fitsio.h"
#include <cstdio>

// Model table details
#define I_n 30
#define T_n 551
#define E_n 30000
#define dlt 0.004 
#define logTMin 5.5
#define maxstoredBins 5

#define CHSPEC_TABLEPATH "/home/mithun/work/ch2/xsm/science/chspec/tables/"
#define CHSPEC_TABLEFILE "chspec_CHIANTIv10.fits.gz" 

int TABLE_READ_STATUS=0;

struct filehandle
{
    // Variable definitions

    double*** allspec;
    double*** allspecBinned;
    int nBinSpec=0;
    int nEbin[maxstoredBins];

    double temp[T_n];
    double ene[E_n];

    double deltaTemp[T_n];
    float logTemp[T_n];

    filehandle()
    {

        for (int i=0; i<maxstoredBins; i++) nEbin[i]=0;

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

        allspecBinned=(double***)malloc(I_n * sizeof(double**));
        if (allspecBinned == NULL){
                fprintf(stderr, "Out of memory");
                exit(0);
        }
        for (int i = 0; i < I_n; i++)
        {
            allspecBinned[i] = (double**)malloc(T_n * sizeof(double*));
            for (int j = 0; j < T_n; j++)
                allspecBinned[i][j] = (double*)malloc(E_n * sizeof(double));
        }

        status = 0;
        char filename[8096];
        sprintf(filename,"%s/%s",CHSPEC_TABLEPATH,CHSPEC_TABLEFILE);
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
        //printf("chspec: Table Data Read!\n");


        // For Gaussian DEM model
        for (int i = 0; i < T_n; i++) {
            logTemp[i]=log10(temp[i]);
            deltaTemp[i]=((pow(10,logTemp[i]+dlt*0.5))-(pow(10,logTemp[i]-dlt*0.5)))/1e6;
        }

    }

    ~filehandle()
    {
        for (int i = 0; i < I_n; i++)
        {
            for (int j = 0; j < T_n; j++) {
                free(allspec[i][j]);
                free(allspecBinned[i][j]);
            }
            free(allspec[i]);
            free(allspecBinned[i]);
        }
        free(allspec);
        free(allspecBinned);
    }

};

    filehandle chfh;

#endif
