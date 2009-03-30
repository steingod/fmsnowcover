/*
 * NAME:
 * probest
 *
 * PURPOSE:
 * To estimate the probability of sea ice, open water and clouds given a
 * set of AVHRR observations, with the possible extension of NWP data
 * using a Bayes approach.
 *
 * REQUIREMENTS:
 * NA
 *
 * INPUT:
 * o Full set of AVHRR observations
 *
 * OUTPUT:
 * o Probabilities of sea ice, open water and clouds
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
 *
 * CVS_ID:
 * $Id: probest.c,v 1.4 2009-03-30 13:42:53 steingod Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include <fmsnowcover.h>

int probest(pinpstr cpa, probstr *p, statcoeffstr cof) {

    double r21, r3a1, d34, r3b1;
    double pa1gi, pa1gc, pa1gf;
    double pr21gi, pr3a1gi, pdtgi, pd34gi, pr3b1gi;
    double pr21gf, pr3a1gf, pdtgf, pd34gf, pr3b1gf;
    double pr21gc, pr3a1gc, pdtgc, pd34gc, pr3b1gc;
    /*double pigr21, pigr31;*/
    double denomsum;
    /*double pice=0.3333, pfree=0.3333, pcloud=0.3333;*/
    double pice=0.5, pfree=0.5, pcloud=0.5;
    /*double picegobs;*/
    
    /*
     * Specify conditional probabilities according to statistical results.
     */
    r21 = cpa.A2/cpa.A1;
    if (cpa.daytime3b) {
      d34 = cpa.T3-cpa.T4;
      r3b1 = cpa.A3b/cpa.A1;
    } else {
      r3a1 = cpa.A3/cpa.A1; 
    }
   
    pa1gi = findprob( cof.ice.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"ice a1");
   
    pr21gi = findprob( cof.ice.r21, r21,"ice r21" );
    
    if (cpa.daytime3b) {
     	pd34gi = findprob( cof.ice.d34, d34,"ice d34" );
	if (pd34gi < 0) pd34gi = 0;
	pr3b1gi = findprob( cof.ice.r3b1, r3b1, "ice r3b1" );
    } else {
   	pr3a1gi = findprob( cof.ice.r3a1, r3a1, "ice r3a1" );
    }
    
    pdtgi = findprob( cof.ice.dt, cpa.tdiff,"ice dt" );

    pa1gc = findprob( cof.cloud.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"cloud a1");

    pr21gc = findprob( cof.cloud.r21, r21,"cloud r21" );
   
    if (cpa.daytime3b) {
	pd34gc = findprob( cof.cloud.d34, d34, "cloud d34" );
	if (pd34gc < 0) pd34gc = 0;
	pr3b1gc = findprob( cof.cloud.r3b1, r3b1, "cloud r3b1" );
    } else {
	pr3a1gc = findprob( cof.cloud.r3a1, r3a1, "cloud r3a1" );
    }

    pdtgc = findprob( cof.cloud.dt, cpa.tdiff,"cloud dt");

    if (cpa.lmask > 1) { /*Landmask recognized land*/
      pa1gf = findprob( cof.land.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"land a1");
      
      pr21gf = findprob( cof.land.r21, r21,"land r21" );

      if (cpa.daytime3b) {
	pd34gf = findprob( cof.land.d34, d34,"land d34" );
	if (pd34gf < 0) pd34gf = 0;
	pr3b1gf = findprob( cof.land.r3b1, r3b1, "land r3b1");
      } else {
 	pr3a1gf = findprob( cof.land.r3a1, r3a1,"land r3a1" );
      }
      
      pdtgf = findprob( cof.land.dt, cpa.tdiff,"land dt" );

    } else {/*Landmask recognized water, or no landmask in use*/
    
      pa1gf=findprob( cof.water.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"water a1");

      pr21gf = findprob( cof.water.r21, r21,"water r21" );
   
      if (cpa.daytime3b) {
	pd34gf = findprob( cof.water.d34, d34,"water d34" );
	if (pd34gf < 0) pd34gf = 0;
	pr3b1gf = findprob( cof.water.r3b1, r3b1, "water r3b1");
      } else {
 	pr3a1gf = findprob( cof.water.r3a1, r3a1,"water r3a1" );
      }
      
      pdtgf = findprob( cof.water.dt, cpa.tdiff,"water dt" );
    }

    /*
     * Use Bayes theorem and estimate probability for ice.
     */
    /*
    picegobs = (pr21gi*pr31gi*pa1gi*pdtgi*pice)/
	((pr21gi*pr31gi*pa1gi*pdtgi*pice)
	 +(pr21gw*pr31gw*pa1gw*pdtgw*pwater)
	 +(pr21gc*pr31gc*pa1gc*pdtgc*pcloud));
    */

    #ifndef FMSNOWCOVER_HAVE_LIBUSENWP
    pdtgi = pdtgc = pdtgf = 1.;
    #endif

    /*Used to easily remove signatures when testing*/
    /*pa1gi = pa1gf = pa1gc = 1.; */
    /*pr21gi= pr21gf= pr21gc= 1.; */
    /*pr3b1gi=pr3b1gf=pr3b1gc=1.; */


    /*This one is to be properly removed when finished testing*/
    pd34gi= pd34gc= pd34gf= 1.;

    if (cpa.daytime3b) {
	denomsum = (pr21gi*pd34gi*pr3b1gi*pa1gi*pdtgi*pice)
	    +(pr21gf*pd34gf*pr3b1gf*pa1gf*pdtgf*pfree)
	    +(pr21gc*pd34gc*pr3b1gc*pa1gc*pdtgc*pcloud);
	p->pice = (pr21gi*pd34gi*pr3b1gi*pa1gi*pdtgi*pice)/denomsum;
	p->pfree = (pr21gf*pd34gf*pr3b1gf*pa1gf*pdtgf*pfree)/denomsum;
	p->pcloud = (pr21gc*pd34gc*pr3b1gc*pa1gc*pdtgc*pcloud)/denomsum;
    } else {
	denomsum = (pr21gi*pr3a1gi*pa1gi*pdtgi*pice)
	    +(pr21gf*pr3a1gf*pa1gf*pdtgf*pfree)
	    +(pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud);
	p->pice = (pr21gi*pr3a1gi*pa1gi*pdtgi*pice)/denomsum;
	p->pfree = (pr21gf*pr3a1gf*pa1gf*pdtgf*pfree)/denomsum;
	p->pcloud = (pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud)/denomsum;
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
