/*
 * NAME:
 * probest
 *
 * PURPOSE:
 * To estimate the probability of sea ice/snow, clouds or open
 * water/clear land given a set of AVHRR observations, with the
 * possible extension of NWP data using a Bayes approach.
 *
 * REQUIREMENTS:
 * NA
 *
 * INPUT:
 * o Full set of AVHRR observations
 * o NWP data
 *
 * OUTPUT:
 * o Probabilities of sea ice/snow, open water/land and clouds
 *
 * NOTES:
 * Definitions of either Gamma or Normal distributions should be put in a
 * specific ASCII file and stored as an etc file of the software. then it
 * is easy to modify coefficients without recompiling the software.
 *
 * BUGS:
 * NA
 *
 * AUTHOR:
 * Øystein Godøy, met.no/FOU, 28.09.2004
 *
 * MODIFIED:
 * Øystein Godøy, met.no/FOU, 18.11.2004: Adapted for Ch3a and normal
 * distribution...
 * Øystein Godøy, METNO/FOU, 30.10.2006: Added the function Hanne Heiberg
 * have been using as probest_hanneh, changed function interface.
 * Øystein Godøy, METNO/FOU, 11.01.2007: Changed some names.
 * Øystein Godøy, METNO/FOU, 02.04.2007: Added 3B support.
 * Mari Anne Killie, METNO/FOU, 31.01.2008: coefficients are removed
 * and placed in ASCII file. Struct of type statcoeffstr now contains
 * the coeffs., and is passed down from avhrrice_pap. Function
 * findprob is added.
 * Mari Anne Killie, METNO/FOU, 13.10.2008: adopting for the
 * reflective part of 3B to be used (introducing r3b1 to replace d34 +
 * renaming existing r31 to r3a1)
 * Mari Anne Killie, METNO/FOU, 08.05.2009: d34 removed, snow introduced.
 * Mari Anne Killie, METNO/FOU, 08.09.2009: using 3 classes over land,
 * 3 over water and 4 over coast.
 * MAK, METNO/FOU, 19.12.2011: suddenly switching off class "snow"
 * over land (near lakes) creates noise -> Reintroducing the option to
 * use snow as 5th class in coastal zone, even though 5 classes did
 * not seem optimal to prevent a row of misclassified pixels along the
 * coast on the first try. Perhaps tuning of FMSNOWSEA and FMSNOWLAND will
 * help. Adding SNOWSWITCH to easily test the effect of the 5th class.
 * MAK, METNO/FOU, 19.12.2011: DTLIM added.
 * 
 * CVS_ID:
 * $Id: probest.c,v 1.11 2013-02-01 10:37:06 mariak Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include <fmsnowcover.h>

#define SNOWSWITCH 0 /*0: no snow in "coast", 1: snow + ice in coast.  */

#define DTLIM 0 /*Should perhaps be moved, but this decides wether
		    dT-signature is used or not! 277 K, approx
		    4celsius. DTLIM 273 used for OSI SAF. To easily
		    remove this test, set DTLIM to 0!*/

/* #undef FMSNOWCOVER_HAVE_LIBUSENWP */
int probest(pinpstr cpa, probstr *p, statcoeffstr cof) {

    double r21, r3a1, r3b1;
    double pa1gi, pr21gi, pr3a1gi, pdtgi, pr3b1gi;
    double pa1gc, pr21gc, pr3a1gc, pdtgc, pr3b1gc;
    double pa1gs, pr21gs, pr3a1gs, pdtgs, pr3b1gs;
    double pa1gw, pr21gw, pr3a1gw, pdtgw, pr3b1gw;
    double pa1gl, pr21gl, pr3a1gl, pdtgl, pr3b1gl;
    double denomsum;
    double pice=0.5, psnow=0.5, pcloud=0.5, pwater=0.5, pland=0.5;
    
    
    /*
     * Specify conditional probabilities according to statistical results.
     */
    r21 = cpa.A2/cpa.A1;
    if (cpa.daytime3b) {
      r3b1 = cpa.A3b/(cpa.A1/cos(fmdeg2rad(cpa.soz)));
    } else {
      r3a1 = cpa.A3/cpa.A1; 
    }

    /*Later, if needing to cut proc.time: sort the calculations below
      under the following categories, i.e.: only calculate those that
      will actually be used*/
    /* if (cpa.lmask <= FMSNOWSEA) { /\*use classes sea, sea ice, cloud*\/ */
    /* } else if (cpa.lmask >= FMSNOWLAND) { /\*classes land, snow, cloud*\/ */
    /* } else { /\*use classes land, sea, snow, sea ice, cloud*\/ */
    /* } */


    /*Ice and snow*/
    pa1gi = findprob( cof.ice.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"ice a1");
    pa1gs = findprob( cof.snow.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"snow a1");
    pr21gi = findprob( cof.ice.r21, r21, "ice r21" );
    pr21gs = findprob( cof.snow.r21, r21, "snow r21" );
    if (cpa.daytime3b) {
      pr3b1gi = findprob( cof.ice.r3b1, r3b1, "ice r3b1" );
      pr3b1gs = findprob( cof.snow.r3b1, r3b1, "snow r3b1" );
    } else {
      pr3a1gi = findprob( cof.ice.r3a1, r3a1, "ice r3a1" );
      pr3a1gs = findprob( cof.snow.r3a1, r3a1, "snow r3a1" );
    }
    pdtgi = findprob( cof.ice.dt, cpa.tdiff, "ice dt" );
    pdtgs = findprob( cof.snow.dt, cpa.tdiff, "snow dt" );
    
    /*Clouds*/
    pa1gc = findprob( cof.cloud.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"cloud a1");
    pr21gc = findprob( cof.cloud.r21, r21, "cloud r21" );
    if (cpa.daytime3b) {
	pr3b1gc = findprob( cof.cloud.r3b1, r3b1, "cloud r3b1" );
    } else {
	pr3a1gc = findprob( cof.cloud.r3a1, r3a1, "cloud r3a1" );
    }
    pdtgc = findprob( cof.cloud.dt, cpa.tdiff, "cloud dt");

    /*Land and water*/
    pa1gl=findprob( cof.land.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"land a1");
    pa1gw=findprob( cof.water.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"water a1");
    pr21gl = findprob( cof.land.r21, r21, "land r21" );
    pr21gw = findprob( cof.water.r21, r21, "water r21" );
    if (cpa.daytime3b) {
      pr3b1gl = findprob( cof.land.r3b1, r3b1, "land r3b1");
      pr3b1gw = findprob( cof.water.r3b1, r3b1, "water r3b1");
    } else {
      pr3a1gl = findprob( cof.land.r3a1, r3a1, "land r3a1" );
      pr3a1gw = findprob( cof.water.r3a1, r3a1, "water r3a1" );
    }
    pdtgl = findprob( cof.land.dt, cpa.tdiff, "land dt" );
    pdtgw = findprob( cof.water.dt, cpa.tdiff,"water dt" );
    

    /*
     * Use Bayes theorem and estimate probability for ice.
     */
    /*
    picegobs = (pr21gi*pr31gi*pa1gi*pdtgi*pice)/
	((pr21gi*pr31gi*pa1gi*pdtgi*pice)
	 +(pr21gw*pr31gw*pa1gw*pdtgw*pwater)
	 +(pr21gc*pr31gc*pa1gc*pdtgc*pcloud));
    */

    /*if (cpa.tdiff == NULL) {*/
    if (cpa.tdiff == 0) {
      pdtgi = pdtgs = pdtgc = pdtgl = pdtgw = 1.;
    }

    /*the pdtgX-test can fail over ice/snow and should only? be used when a positive model temperature. Fails over Greenland, but is needed over Norway..  Tdiff = T_model - T_4*/
    if(cpa.tdiff + cpa.T4 < DTLIM) {
      pdtgi  = pdtgs  = pdtgc  = pdtgl  = pdtgw  = 1.;
    }


    /*Used to easily remove signatures when testing*/
    /*pa1gi  = pa1gs  = pa1gc  = pa1gl  = pa1gw  = 1.;*/
    /*pr21gi = pr21gs = pr21gc = pr21gl = pr21gw = 1.;*/
    /*pr3b1gi= pr3b1gs= pr3b1gc= pr3b1gl= pr3b1gw= 1.;*/
    /*pr3a1gi= pr3a1gs= pr3a1gc= pr3a1gl= pr3a1gw= 1.;*/
    /*pdtgi  = pdtgs  = pdtgc  = pdtgl  = pdtgw  = 1.;*/

    if (cpa.lmask <= FMSNOWSEA) { /*Water: use classes sea ice/water/cloud*/
      if (cpa.daytime3b) {
	denomsum = (pr21gi*pr3b1gi*pa1gi*pdtgi*pice)
	          +(pr21gw*pr3b1gw*pa1gw*pdtgw*pwater)
	          +(pr21gc*pr3b1gc*pa1gc*pdtgc*pcloud);
	p->pice = (pr21gi*pr3b1gi*pa1gi*pdtgi*pice)/denomsum;
	p->pfree =(pr21gw*pr3b1gw*pa1gw*pdtgw*pwater)/denomsum;
	p->pcloud=(pr21gc*pr3b1gc*pa1gc*pdtgc*pcloud)/denomsum;
      } else {
	denomsum = (pr21gi*pr3a1gi*pa1gi*pdtgi*pice)
	  	  +(pr21gw*pr3a1gw*pa1gw*pdtgw*pwater)
	          +(pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud);
	p->pice = (pr21gi*pr3a1gi*pa1gi*pdtgi*pice)/denomsum;
	p->pfree =(pr21gw*pr3a1gw*pa1gw*pdtgw*pwater)/denomsum;
	p->pcloud=(pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud)/denomsum;
      }
    } else if (cpa.lmask >= FMSNOWLAND) { /*Land: use classes snow/land/cloud*/
      if (cpa.daytime3b) {
	denomsum = (pr21gs*pr3b1gs*pa1gs*pdtgs*psnow)
	          +(pr21gl*pr3b1gl*pa1gl*pdtgl*pland)
	    	  +(pr21gc*pr3b1gc*pa1gc*pdtgc*pcloud);
	p->pice = (pr21gs*pr3b1gs*pa1gs*pdtgs*psnow)/denomsum;
	p->pfree =(pr21gl*pr3b1gl*pa1gl*pdtgl*pland)/denomsum;
	p->pcloud=(pr21gc*pr3b1gc*pa1gc*pdtgc*pcloud)/denomsum;
      } else {
	denomsum = (pr21gs*pr3a1gs*pa1gs*pdtgs*psnow)
	          +(pr21gl*pr3a1gl*pa1gl*pdtgl*pland)
	          +(pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud);
	p->pice = (pr21gs*pr3a1gs*pa1gs*pdtgs*psnow)/denomsum;
	p->pfree =(pr21gl*pr3a1gl*pa1gl*pdtgl*pland)/denomsum;
	p->pcloud=(pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud)/denomsum;
      }
    } else { /*Coast: use classes ice/land/snow/cloud/water*/
 
     if (cpa.daytime3b) {
	denomsum = (pr21gi*pr3b1gi*pa1gi*pdtgi*pice)
	  +(pr21gs*pr3b1gs*pa1gs*pdtgs*psnow)*SNOWSWITCH
	  +(pr21gl*pr3b1gl*pa1gl*pdtgl*pland)
	  +(pr21gw*pr3b1gw*pa1gw*pdtgw*pwater)
	  +(pr21gc*pr3b1gc*pa1gc*pdtgc*pcloud);
	p->pice = ((pr21gi*pr3b1gi*pa1gi*pdtgi*pice)
		   +(pr21gs*pr3b1gs*pa1gs*pdtgs*psnow)*SNOWSWITCH)/denomsum;
	p->pfree =((pr21gl*pr3b1gl*pa1gl*pdtgl*pland)
		   +(pr21gw*pr3b1gw*pa1gw*pdtgw*pwater))/denomsum;
	p->pcloud = (pr21gc*pr3b1gc*pa1gc*pdtgc*pcloud)/denomsum;
      } else {
	denomsum = (pr21gi*pr3a1gi*pa1gi*pdtgi*pice)
	  +(pr21gs*pr3a1gs*pa1gs*pdtgs*psnow)*SNOWSWITCH
	  +(pr21gl*pr3a1gl*pa1gl*pdtgl*pland)
	  +(pr21gw*pr3a1gw*pa1gw*pdtgw*pwater)
	  +(pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud);
	p->pice = ((pr21gi*pr3a1gi*pa1gi*pdtgi*pice)
		   +(pr21gs*pr3a1gs*pa1gs*pdtgs*psnow)*SNOWSWITCH)/denomsum;
	p->pfree =((pr21gl*pr3a1gl*pa1gl*pdtgl*pland)
		   +(pr21gw*pr3a1gw*pa1gw*pdtgw*pwater))/denomsum;
	p->pcloud = (pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud)/denomsum;
     }

     /*  if (cpa.daytime3b) { /\*tester igjen 5 klasser i kystsoner*\/ */
/*        denomsum = (pr21gs*pr3b1gs*pa1gs*pdtgs*psnow) */
/* 	  +(pr21gi*pr3b1gi*pa1gi*pdtgi*pice) */
/* 	  +(pr21gl*pr3b1gl*pa1gl*pdtgl*pland) */
/* 	  +(pr21gw*pr3b1gw*pa1gw*pdtgw*pwater) */
/* 	  +(pr21gc*pr3b1gc*pa1gc*pdtgc*pcloud); */
/*         p->pice =((pr21gi*pr3b1gi*pa1gi*pdtgi*pice) */
/* 		  +(pr21gs*pr3b1gs*pa1gs*pdtgs*psnow))/denomsum; */
/* 	p->pfree =((pr21gl*pr3b1gl*pa1gl*pdtgl*pland) */
/* 		   +(pr21gw*pr3b1gw*pa1gw*pdtgw*pwater))/denomsum; */
/* 	p->pcloud = (pr21gc*pr3b1gc*pa1gc*pdtgc*pcloud)/denomsum; */
/*       } else { */
/* 	denomsum = (pr21gs*pr3a1gs*pa1gs*pdtgs*psnow) */
/* 	  +(pr21gi*pr3a1gi*pa1gi*pdtgi*pice) */
/* 	  +(pr21gl*pr3a1gl*pa1gl*pdtgl*pland) */
/* 	  +(pr21gw*pr3a1gw*pa1gw*pdtgw*pwater) */
/* 	  +(pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud); */
/* 	p->pice =((pr21gi*pr3a1gi*pa1gi*pdtgi*pice) */
/* 		  +(pr21gs*pr3a1gs*pa1gs*pdtgs*psnow))/denomsum; */
/* 	p->pfree =((pr21gl*pr3a1gl*pa1gl*pdtgl*pland) */
/* 		   +(pr21gw*pr3a1gw*pa1gw*pdtgw*pwater))/denomsum; */
/* 	p->pcloud = (pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud)/denomsum; */
/*      } */

    }

    /*
     * Remove ice from areas obviously open water. This could generate
     * noise over melt ponds and during spring season, but not
     * necessarily...
     */
    /*
    if ((cpa.T4 > 276.15) || (cpa.T4-cpa.T5 > 2.0) || (cpa.tdiff > 5.)) {
    picegobs = 0.;
    }
    */
    /*
    if ((cpa.T4-cpa.T5 > 2.0)) picegobs = 0.;
    */
    
    /*
    pigr21 = (pr21gi*pice)/
	((pr21gi*pice)
	 +(pr21gw*pwater)
	 +(pr21gc*pcloud));

    pigr31 = (pr31gi*pice)/
	((pr31gi*pice)
	 +(pr31gw*pwater)
	 +(pr31gc*pcloud));

    picegobs = pigr21*pigr31;
    picegobs = pigr31;
    */
	
    return(FM_OK);
}

/*
 * NAME:
 * findprob
 *
 * PURPOSE:
 * calculates the probability of <feature> given <surface> for the
 * observed value x, using either normalpdf, gammapdf og gammapdf3par
 *
 * INPUT:
 * 1: featstr struct containing the statistical coefficients for the
 * particular surface and feature + a keyword to determine which pdf
 * to use
 * 2: the observed value x
 * 3: a string with name of surface and feature in case of error 
 *
 * OUTPUT:
 * o the probability
 */
double findprob(featstr feat, double x, char *whereami) {
  double pdf;
  char *where="findpdf";
  char what[FMSNOWCOVER_MSGLENGTH];
  int errflg;
  errflg = 0;
  
  if (!feat.count) {
    sprintf(what,"Stat. coefficients have not been read for %s",whereami);
    errflg++;
  }

  else if (feat.count == 1){
    if (feat.key == 'n') pdf = normalpdf(feat.par1, feat.par2, x);
    else if (feat.key == 'g') pdf = gammapdf(feat.par1, feat.par2, x);
    /*else if (feat.key == 't') pdf = gammapdf3par(feat.par1,feat.par2,feat.par3,x);*/
    else {
      sprintf(what,"Could not recognize pdf routine key for %s",whereami);
      errflg++;
    }
  }

  else if (feat.count > 1) {
      sprintf(what,"Stat. coefficients for %s are read more than once\n",
	      whereami);
      errflg++;
    }

  if (errflg) {
    fmerrmsg(where,what);
    exit(FM_IO_ERR);
  }
  
  return(pdf);
}
