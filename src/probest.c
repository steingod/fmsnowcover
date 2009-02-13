/* *****************************************************************
 * COPYRIGHT:
 * EUMETSAT
 *
 * PRODUCED BY:
 * Norwegian Meteorological Institute (met.no)
 * Research and Development Department
 * P.O.BOX 43 - Blindern, N-0313 OSLO, NORWAY
 *
 * This SW was developed by met.no and the Danish Meteorological
 * Institute (DMI) within the context of the Co-operation Agreement
 * for the development of a pilot SAF on Ocean and Sea Ice.
 * *****************************************************************/
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
 * $Id: probest.c,v 1.1 2009-02-13 23:23:14 steingod Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include <avhrrice_pap.h>

int probest(pinpstr cpa, probstr *p, statcoeffstr cof) {

    double r21, r3a1, d34, r3b1;
    double pa1gi, pa1gc, pa1gw;
    double pr21gi, pr3a1gi, pdtgi, pd34gi, pr3b1gi;
    double pr21gw, pr3a1gw, pdtgw, pd34gw, pr3b1gw;
    double pr21gc, pr3a1gc, pdtgc, pd34gc, pr3b1gc;
    /*double pigr21, pigr31;*/
    double denomsum;
    /*double pice=0.3333, pwater=0.3333, pcloud=0.3333;*/
    double pice=0.5, pwater=0.5, pcloud=0.5;
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


    pa1gw = findprob( cof.water.a1, cpa.A1/cos(fmdeg2rad(cpa.soz)),"water a1");

    pr21gw = findprob( cof.water.r21, r21,"water r21" );
   
    if (cpa.daytime3b) {
	pd34gw = findprob( cof.water.d34, d34,"water d34" );
	if (pd34gw < 0) pd34gw = 0;
	pr3b1gw = findprob( cof.water.r3b1, r3b1, "water r3b1");
    } else {
 	pr3a1gw = findprob( cof.water.r3a1, r3a1,"water r3a1" );
    }
    
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


    #ifndef AVHRRICE_HAVE_NWP
    pdtgi = pdtgc = pdtgw = 1.;
    #endif

    /*Used to easily remove signatures when testing*/
    /*pa1gi = pa1gw = pa1gc = 1.; */
    /*pr21gi= pr21gw= pr21gc= 1.; */
    /*pr3b1gi=pr3b1gc=pr3b1gw=1.; */


    /*This one is to be properly removed when finished testing*/
    pd34gi= pd34gc= pd34gw= 1.;


    if (cpa.daytime3b) {
	denomsum = (pr21gi*pd34gi*pr3b1gi*pa1gi*pdtgi*pice)
	    +(pr21gw*pd34gw*pr3b1gw*pa1gw*pdtgw*pwater)
	    +(pr21gc*pd34gc*pr3b1gc*pa1gc*pdtgc*pcloud);
	p->pice = (pr21gi*pd34gi*pr3b1gi*pa1gi*pdtgi*pice)/denomsum;
	p->pwater = (pr21gw*pd34gw*pr3b1gw*pa1gw*pdtgw*pwater)/denomsum;
	p->pcloud = (pr21gc*pd34gc*pr3b1gc*pa1gc*pdtgc*pcloud)/denomsum;
    } else {
	denomsum = (pr21gi*pr3a1gi*pa1gi*pdtgi*pice)
	    +(pr21gw*pr3a1gw*pa1gw*pdtgw*pwater)
	    +(pr21gc*pr3a1gc*pa1gc*pdtgc*pcloud);
	p->pice = (pr21gi*pr3a1gi*pa1gi*pdtgi*pice)/denomsum;
	p->pwater = (pr21gw*pr3a1gw*pa1gw*pdtgw*pwater)/denomsum;
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
	
    return(0);
}



int probest_hanneh(pinpstr cpa, probstr *p) {

    /*char *where="probest_hanneh";*/
    double a21,a31;
    double pa21gi, pa21gw, pa21gc, pa31gi, pa31gw, pa31gc;
    /* if pice = pwater = pcloud, these factors cancel in Bayes  */
    /* formula and are not needed: */
    /* double pice=1.0/3.0, pwater=1.0/3.0, pcloud=1.0/3.0; */
    double denomsum;
    int debug=0 /* 0 = false */;
    
    /*
     * Estimate A3a/A1 and A2/A1 for use as input to PDF functions.
     */
    if (cpa.A3 < 0.0 || cpa.A2 < 0.0 || cpa.A1 < 0.0){
      if (debug) {
	printf("%s %7.3f %7.3f %7.3f", "A1 A2 A3 ", cpa.A1,cpa.A2,cpa.A3);
      }
      return(1);
    }
    a21 = (double) (cpa.A2 / cpa.A1);
    a31 = (double) (cpa.A3 / cpa.A1);

    /*
     * Estimate the probability of getting A2/A3 and A3/A1 given ice, cloud or water
     */

    pa21gi = normalpdf(0.747654, 0.0616349,a21);
    /*    pa21gi = gammapdf3(16.38219,0.03767022,.1,a21); */
    /*    printf(" %s%f\n"," pa21gi =",pa21gi); */
    /*     if (pa21gi < 0) { */
    /* 	fprintf(stderr,"%s%s\n",progerr,"Error in distribution function a21i."); */
    /* 	return(1); */
    /*     } */

    pa21gw = normalpdf(0.404327,0.0572932,a21);
    /*     pa21gw = gammapdf3(13.77803,0.01353032,0.25,a21); */
    /*     printf(" %s%f\n"," pa21gw =",pa21gw); */
    /*     if (pa21gw < 0) { */
    /* 	fprintf(stderr,"%s%s\n",progerr,"Error in distribution function a21w."); */
    /* 	return(1); */
    /*     } */
   
    pa21gc = normalpdf(0.813544, 0.064344, a21);
    /*     pa21gc = gammapdf3(8.400453,0.01743947,0.4,a21);*/
    /*     printf(" %s%f\n"," pa21gc =",pa21gc); */
    /*     if (pa21gc < 0) { */
    /* 	fprintf(stderr,"%s%s\n",progerr,"Error in distribution function a21c."); */
    /* 	return(1); */
    /*     }  */

    pa31gi = normalpdf(0.095215, 0.03729273, a31);
    /*    pa31gi = gammapdf3(2.0891368,0.03008235, 0.0 ,a31);*/
    /*     printf(" %s%f\n"," pa31gi =",pa31gi);   */
    /*     if (pa31gi < 0) { */
    /* 	fprintf(stderr,"%s%s\n",progerr,"Error in distribution function a31i."); */
    /* 	return(1); */
    /*     } */
 
    pa31gw = normalpdf(0.07901805, 0.06690694, a31);
    /*    pa31gw = gammapdf3(3.56321, 0.03292881, -0.01,a31); */
    /*     printf(" %s%f\n"," pa31gw =",pa31gw);  */
    /*     if (pa31gw < 0) { */
    /*       fprintf(stderr,"%s%s\n",progerr,"Error in distribution function a31w."); */
    /*       return(1); */
    /*     } */
    
    pa31gc = normalpdf(0.6288227, 0.2006739, a31);
    /*    pa31gc = gammapdf3(44.04446,0.02609957,-0.5,a31);*/
    /*     printf(" %s%f\n"," pa31gc =",pa31gc);  */
    /*     if (pa31gc < 0) { */
    /*       fprintf(stderr,"%s%s\n",progerr,"Error in distribution function a31c."); */
    /*       return(1); */
    /*     } */

    /*
     * Use Bayes theorem and estimate probability for ice, water and cloud.
     */

    /* 
     * Use this formula if EQUAL pice, pwater and pcloud,
     * i.e. assume equal chance of ice, water, and cloud 
     */
    denomsum = (pa21gi*pa31gi)+(pa21gw*pa31gw)+(pa21gc*pa31gc);
    /*
    *picegobs = (float)((pa21gi*pa31gi) / denomsum);
    *pwgobs =   (float)((pa21gw*pa31gw) / denomsum);
    *pcgobs =   (float)((pa21gc*pa31gc) / denomsum);
    */
    p->pice = (float)((pa21gi*pa31gi) / denomsum);
    p->pwater =   (float)((pa21gw*pa31gw) / denomsum);
    p->pcloud =   (float)((pa21gc*pa31gc) / denomsum);

    /* 
     *  Use this formula if pice, pwater and pcloud ARE NOT EQUAL
     */
    /*     denomsum = (pa21gi*pa31gi*pice)+(pa21gw*pa31gw*pwater)+(pa21gc*pa31gc*pcloud); */
    /*     picegobs = (pa21gi*pa31gi*pice)  / denomsum; */
    /*     pwgobs =   (pa21gw*pa31gw*pwater)/ denomsum; */
    /*     pcgobs =   (pa21gc*pa31gc*pcloud)/ denomsum; */

/*     printf(" %s %6.3f %6.3f %6.3f %6.3f %s %6.3f %6.3f %6.3f %6.3f %s %6.3f %6.3f %6.3f\n",  */
/* 	   "a21 pgiwc ",  a21, pa21gi, pa21gw, pa21gc,  */
/* 	   " a31 pgiwc ", a31, pa31gi, pa31gw, pa31gc, */
/* 	   " piwc ",picegobs, pwgobs, pcgobs); */

    return(0);
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
  char what[OSI_MSGLENGTH];
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
    exit(0);
  }
  
  return(pdf);
}
