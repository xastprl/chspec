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
      Aug,2020 :   Biswajit 
      Oct,2020 :   Santosh, Mithun (modified)
      Mar,2021 :   Biswajit (Updated for CHIANTI v10)
      June,2021:   Biswajit, (Generalized for the use of other instruments.)
      Nov 2021 :   Mithun: Fixed bug in generalization, rebinned table stored in memory to make the calculation 
                           faster
      Nov 2021 :   Mithun: Spectra from Single Gaussian DEM model added
      Sep 2022 :   Mithun: Restructured the codes to have combined isothermal and gaussian DEM models into a 
                           single chspec package

*/

#include "readTable.h"

//if(TABLE_READ_STATUS==0)
//{
//    filehandle chfh;
//    TABLE_READ_STATUS=1;
//}

extern "C"
void chisoth(const RealArray& energyArray, const RealArray& parms, 
               int spectrum, RealArray& flux, RealArray& fluxErr, 
               const string& init)

{
    double T,T1, mf, ad, tL, tR, eL, eR, dydx,dTRL,dTL ;
    float abund[I_n] ;
    int nE, en, inL, inR,ene_ind;
    int LowerBound_ind, Ebound_error;
    
    T1=parms[0] ;// Log(T)
    T=pow(10,T1); // in K

    nE = energyArray.size()-1;
    flux.resize(nE);
    for (en = 0; en < nE; en++) flux[en]=0.0;

    for (int p=1; p<I_n; p++) {abund[p]= parms[p];}
    abund[0]=12.0;
 
    inL = 0;

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

    // Check if energy bins match previous stored ones
    int matchBinSpec=-1;
    for (int i=0;i<chfh.nBinSpec;i++)
    {
        if(nE==chfh.nEbin[i])
            matchBinSpec=i;
    }
   
    if(matchBinSpec==-1) // Response energy bins are not same as previously stored ones
    {

        printf("Rebinning table to response energy bins: NEW NBINS %d OLD NBINS %d\n",nE,chfh.nEbin[chfh.nBinSpec]);

        LowerBound_ind = 0;
        Ebound_error = 0;

        for (LowerBound_ind = 0; LowerBound_ind < E_n;LowerBound_ind++)
        {
            if (fabs(chfh.ene[LowerBound_ind] - energyArray[0]) < 1.0e-6)
           {
            Ebound_error = 0;
            break;
           }
           else
                Ebound_error = 1;
        }

        if (int(Ebound_error) == 1)
        {
           printf("chspec Error: Lower bound of the energy response should be >= 100 eV or should be \
                    incremented as an integer multiple of 1 eV. ...Terminating the program.\n");
            exit(0);
        }

        // Rebin spectrum to required energy bins

        //printf("Lower bound %d\n",LowerBound_ind);

        int enBinOffset=0;
        for (int i=0;i<chfh.nBinSpec;i++) enBinOffset+=chfh.nEbin[i];

        for (int i = 0; i < I_n; i++)
        {
            for (int iT=0;iT<T_n;iT++)
            {
                ene_ind = LowerBound_ind;

                for (en = 0; en < nE; en++)
                {
                    chfh.allspecBinned[i][iT][en+enBinOffset]=0;
                    while(chfh.ene[ene_ind] >= energyArray[en] and chfh.ene[ene_ind] < energyArray[en+1])
                    {
                        chfh.allspecBinned[i][iT][en+enBinOffset]+=chfh.allspec[i][iT][ene_ind];
                        ene_ind = ene_ind+1;
                    }
                }
            }

        }
        chfh.nEbin[chfh.nBinSpec]=nE;
        matchBinSpec=chfh.nBinSpec;
        chfh.nBinSpec+=1;

    }

    int enBinOffset=0;
    for (int i=0;i<matchBinSpec;i++) enBinOffset+=chfh.nEbin[i];

    //printf("%ld\t%d\t%d\t%d\n",nE,chfh.nBinSpec,matchBinSpec,enBinOffset,chfh.nEbin[matchBinSpec]);

    for (int i = 0; i < I_n; i++)
    {
        ad = abund[i] - abund[0];
        mf = pow(10, ad);
        for (en = 0; en < nE; en++)
        {
            eR=chfh.allspecBinned[i][inR][en+enBinOffset];
            eL=chfh.allspecBinned[i][inL][en+enBinOffset];
            dydx= (eR-eL)/dTRL;
            flux[en]+= mf*(eL + dydx*dTL);
        }
    }

}

extern "C"
void chgausdem(const RealArray& energyArray, const RealArray& parms, 
               int spectrum, RealArray& flux, RealArray& fluxErr, 
               const string& init)

{

    int nE,en,ene_ind;
    float abund[I_n];
    double ad,mfArray[I_n];
    int LowerBound_ind, Ebound_error;

    double totEM[T_n],maxEM=0;
    float logT,sigma;

    for (int i=0;i<T_n;i++) totEM[i]=0.0;

    nE = energyArray.size()-1;
    flux.resize(nE);
    for (en = 0; en < nE; en++) flux[en]=0.0;
    
    for (int p=0; p<I_n; p++) {abund[p]= parms[p];}

//    for (int p=0; p<nGaus; p++) gausDEM[p]=parms[p+I_n];
    logT=parms[I_n];
    sigma=parms[I_n+1];    
   
    for (int i=0;i<T_n;i++)
    {
        totEM[i]=exp(-((chfh.logTemp[i]-logT)*(chfh.logTemp[i]-logT))/(2.0*sigma*sigma))*chfh.deltaTemp[i];
        if(totEM[i]>maxEM) maxEM=totEM[i];
    }
    
    //for (int i=0;i<T_n;i++) printf("%lf\n",totEM[i]);

    // Check if energy bins match previous stored ones
    int matchBinSpec=-1;
    for (int i=0;i<chfh.nBinSpec;i++)
    {
        if(nE==chfh.nEbin[i])
            matchBinSpec=i;
    }


    if(matchBinSpec==-1) // Response energy bins are not same as previously stored ones
    {

        printf("Rebinning table to response energy bins\n");

        LowerBound_ind = 0;
        Ebound_error = 0;
 
        for (LowerBound_ind = 0; LowerBound_ind < E_n;LowerBound_ind++)
        {
            if (fabs(chfh.ene[LowerBound_ind] - energyArray[0]) < 1.0e-6) 
 	       {   
 	   	    Ebound_error = 0;
 	   	    break;
 	       }
 	       else 
                Ebound_error = 1;
        }
 
        if (int(Ebound_error) == 1) 
        {
 	       printf("chspec Error: Lower bound of the energy response should be >= 100 eV or should be \
                    incremented as an integer multiple of 1 eV. ...Terminating the program.\n");
            exit(0);
        } 

        // Rebin spectrum to required energy bins

        //printf("Lower bound %d\n",LowerBound_ind);

        int enBinOffset=0;
        for (int i=0;i<chfh.nBinSpec;i++) enBinOffset+=chfh.nEbin[i];

        for (int i = 0; i < I_n; i++)
        {
            for (int iT=0;iT<T_n;iT++)
            {
                ene_ind = LowerBound_ind;

                for (en = 0; en < nE; en++)
                {
                    chfh.allspecBinned[i][iT][en+enBinOffset]=0;
                    while(chfh.ene[ene_ind] >= energyArray[en] and chfh.ene[ene_ind] < energyArray[en+1])
                    {
                        chfh.allspecBinned[i][iT][en+enBinOffset]+=chfh.allspec[i][iT][ene_ind];
                        ene_ind = ene_ind+1;
                    }
                }
            }

        }
        chfh.nEbin[chfh.nBinSpec]=nE;
        matchBinSpec=chfh.nBinSpec;
        chfh.nBinSpec+=1;

    }

    int enBinOffset=0;
    for (int i=0;i<matchBinSpec;i++) enBinOffset+=chfh.nEbin[i];

    //for (en = 0; en < nE; en++) printf("%lf\n",chfh.allspecBinned[25][400][en]);

    for (int i = 0; i < I_n; i++)
	{
	    ad = abund[i] - abund[0];
	    mfArray[i] = pow(10, ad);
    }
	    
    for (int iT=0;iT<T_n;iT++)
    {
        if((totEM[iT]/maxEM) > 1e-10)
        {
            for (int i = 0; i < I_n; i++)
            {
	            for (en = 0; en < nE; en++)
                    flux[en]+=mfArray[i]*chfh.allspecBinned[i][iT][en+enBinOffset]*totEM[iT];
            }
        }
    }


}

