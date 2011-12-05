/*
 * NAME:
 * NA
 *
 * PURPOSE:
 * Pixelwise processing of snow parameters using AVHRR data as input.
 *
 * Index numbers correspond to specific surface classes. Se code for
 * further description of which index numbers that are used as these
 * may change frequently according to specific needs. For the time
 * the program write TIFF files for presentation program gl_image
 * which interpretes 10 different classes and uses an internal colortable
 * with indexes 0-9.
 *
 * REQUIREMENTS:
 * NA
 *
 * INPUT:
 * img - AVHRR image data and header
 * cmask - cloud mask
 * lmask - land/sea mask
 * algo - flag determining whether night time or day time data are used
 *
 * OUTPUT:
 * pice - Probability of ice given the AVHRR observations
 * pfree - Probability of open water or land given the AVHRR observations
 * pcloud - Probability of cloud given the AVHRR observations
 * classed - Classed ice probability
 *
 * NOTES:
 * Use of algo variable changed, code should be checked for unnecessary use of
 * this variable. This is to make it easier to implement a quality flag and to
 * ensure that pixels with enough sun is processed instead of discarding the
 * whole tile...
 *
 * A hack to handle saturation problems within 3A is imlemented. This
 * should be handled more properly by the preprocessing in time.
 *
 * BUGS:
 * NA
 *
 * AUTHOR:
 * Øystein Godøy, DNMI/FOU, 28/04/1999
 *
 * MODIFIED:
 * Øystein Godøy, met.no/FOU, 27.09.2004
 * Modification for full Bayes approach started...
 * Øystein Godøy, met.no/FOU, 03.11.2004
 * Testing glm fit, , testing ended but code not removed
 * Øystein Godøy, METNO/FOU, 30.10.2006: Adapting code to results achieved
 * by Vibeke W. Thyness and Hanne Heiberg.
 * Øystein Godøy, METNO/FOU, 11.01.2007: Minor changes, added missing
 * values due to 3A trouble.
 * Øystein Godøy, METNO/FOU, 02.04.2007: Adding 3B support based upon old
 * training data.
 * Mari Anne Killie, METNO/FOU, 31.01.2008: small modification to
 * allow statcoeffstr being passed from avhrrice_pap.c down to probest.c
 * Mari Anne Killie, METNO/FOU, 13.10.2008: satimg removed +
 * introducing fm_ch3brefl.
 * MAK, METNO/FOU, 22.09.2009: (temp.) adding "cat" to categorize each
 * pixel in class with highest probability.
 *
 * CVS_ID:
 * $Id: pix_proc.c,v 1.10 2011-12-05 09:58:47 mariak Exp $
 */ 

#include <fmsnowcover.h>
/*#undef FMSNOWCOVER_HAVE_LIBUSENWP*/

int process_pixels4ice(fmio_img img, unsigned char *cmask[], 
       unsigned char *lmask, nwpice nwp, datafield *probs, 
       unsigned char *class, unsigned char *cat, short algo, statcoeffstr cof) {
    
    char *where="process_pixels4ice";
    char what[FMSNOWCOVER_MSGLENGTH];
    int i, j, size;
    int xc, yc;
    /* double x; */
    pinpstr cpar;
    probstr p;
    fmucsref ucs0;
   /*  fmxy cart; */
    fmindex cart;
    fmucspos ucspos;
    fmgeopos geop;
    fmtime timeid;
    fmsec1970 timeidsec, tst;
    float zsun;
    int doy;
    fmscale calib; /*will contain gain and intercept values*/

    fmlogmsg(where,
	    "Now processing the individual pixels to gain ice probability...");
    /*
     * Convert structures to lesser units for later use...
     */
    size = img.iw*img.ih;

    ucs0.Ax = img.Ax;
    ucs0.Ay = img.Ay;
    ucs0.Bx = img.Bx;
    ucs0.By = img.By;
    ucs0.iw = img.iw;
    ucs0.ih = img.ih;
    
    timeid.fm_year = img.yy;
    timeid.fm_mon = img.mm;
    timeid.fm_mday = img.dd;
    timeid.fm_hour = img.ho;
    timeid.fm_min = img.mi;
    timeid.fm_sec = 0;
    timeidsec = tofmsec1970(timeid);
    
    cpar.A1 = -999.;
    cpar.A2 = -999.;
    cpar.A3 = -999.;
    cpar.A3b= -999.;
    cpar.T3 = -999.;
    cpar.T4 = -999.;
    cpar.T5 = -999.;
    cpar.soz= -99;
    cpar.saz= -99;
    cpar.tdiff=-99;
    cpar.algo = algo;

    fm_img2slopes(img,&calib); /*collects gain and intercept*/

    doy = fmdayofyear(timeid);

    /*
     * Start of nested loops that run through alle pixels.
     */
    for (yc=0; yc < img.ih; yc++) {
	for (xc=0; xc < img.iw; xc++) {

	    /*
	     * 2D -> 1D indexing...
	     */
	    i=fmivec(xc, yc, img.iw);
	    class[i] = 0;
	    cat[i] = 5; /*undef.*/

	    /*
	     * Quick exit when test processing large tile (not
	     * interested in all of the tile)
	     */
	    /*if (i > 10000000) {return(10);}*/

	    /*
	     * Better safe than sorry...
	     */
	    for (j=0; j<FMSNOWCOVER_OLEVELS; j++) {
		((float *) probs[j].data)[i] = FMSNOWCOVERMISVAL_NOCOV;
	    }

	    /*
	     * Estimate solar zenith angle for each pixel.
	     */
	    cart.row = yc;
	    cart.col = xc;
	    ucspos = fmind2ucs(ucs0, cart);
	    geop = fmucs2geo(ucspos,MI);
	    /*tst needed to compensate for changes in fmsolarzenith:*/
	    tst = fmutc2tst(timeidsec, geop.lon);
	    zsun = fmsolarzenith(tst, geop);


	    if (zsun < FMSNOWSUNZEN) {
		cpar.algo = 2;
	    } else {
		class[i] = 0;
		for (j=0; j<FMSNOWCOVER_OLEVELS; j++) {
		    ((float *) probs[j].data)[i] = FMSNOWCOVERMISVAL_NIGHT;
		}
		continue;
	    }
    
	    /*
	     * Classification is only performed if the infrared channels
	     * are available and only for the parts of the image were
	     * the satellite has passed. 
	     */
	    if ((img.image[3][i] == 0) && (img.image[4][i] == 0)) {
		for (j=0; j<FMSNOWCOVER_OLEVELS; j++) {
		    ((float *) probs[j].data)[i] = FMSNOWCOVERMISVAL_NOCOV;
		}
		continue;
	    }

	    /*
	     * Handle 3B if 3A is missing
	     */
	    cpar.daytime3b = 0;
	    if (cpar.algo == 2) {
		if (img.image[2][i] > 0 && img.image[5][i] == 0) {
		    cpar.daytime3b = 1;
		}
	    }

	    cpar.lmask = 0;
	    if (lmask == NULL) {
		if (i == 0) {
		    fmlogmsg(where,
			"Landmask not in use, using coefficients for sea/ice/cloud");
		}
	    }
	    else {
		cpar.lmask = (short) lmask[i];
	    }

	    /*
	     * Added hack on 3A due to saturation problems...
	     */
	    if (!cpar.daytime3b) {
		if ((img.image[5][i] == 0) && (img.image[3][i] > 50)) {
		    class[i] = 0;
		    for (j=0; j<FMSNOWCOVER_OLEVELS; j++) {
			((float *) probs[j].data)[i] = FMSNOWCOVERMISVAL_3A;
		    }
		    continue;
		}
	    }

	    /*
	     * Estimate geophysical parameters to be used in classification.
	     * At present only a nighttime algoritm utilizing infrared 
	     * satellite information which are converted to 
	     * brightnesstemperatures. T3, T4 and T5 is the brightness-
	     * temperatures of AVHRR channels 3, 4 and 5 respectively.
	     * Dimension of temperatures are Kelvin.
	     */

	    if (cpar.algo == 2 && img.z > 3) {
	        cpar.A1 = fm_byte2float(img.image[0][i], calib, "Reflectance");
		cpar.A2 = fm_byte2float(img.image[1][i], calib, "Reflectance");
		cpar.A3 = fm_byte2float(img.image[5][i], calib, "Reflectance");
	    }
	    cpar.T3 = fm_byte2float(img.image[2][i], calib, "Temperature");
	    cpar.T4 = fm_byte2float(img.image[3][i], calib, "Temperature");
	    cpar.T5 = fm_byte2float(img.image[4][i], calib, "Temperature");
	    cpar.soz = zsun;
	    cpar.saz = 0.;

	    cpar.tdiff = 0.0;	  
	    #ifdef FMSNOWCOVER_HAVE_LIBUSENWP
	    cpar.tdiff = nwp.t0m[i]-cpar.T4;
            #endif

	    /* Estimate the reflective part of daytime channel 3b */
	    if (cpar.daytime3b){
	      cpar.A3b = fm_ch3brefl(cpar.T3,cpar.T4,cpar.soz,img.sa,doy);
	    }
	    

	    if (i == 0) {
	     fmlogmsg(where,"Using probest to estimate pixel probabilities...");
	    } 
	    if (probest(cpar, &p, cof)) {
		sprintf(what,
			"Something went wrong in pixel processing of %d",i);
		fmerrmsg(where,what);
	    }
	    
	    /*
	     * Adding this to prevent classification when probabilities do
	     * not sum to 1.
	     */
	    if (p.pice+p.pfree+p.pcloud<0.95 || p.pice+p.pfree+p.pcloud>1.05){
		continue; 
	    }
 
	    /*
	     * Also, if r3b1 is too large the probabilities can end up as
	     * nan (not fixed by statement above). Trying this:
	     */
	    if (isnan(p.pice) || isnan(p.pfree) || isnan(p.pcloud)) {
		continue;
	    }

	    ((float *) probs[0].data)[i] = p.pice;
	    ((float *) probs[1].data)[i] = p.pfree;
	    ((float *) probs[2].data)[i] = p.pcloud;

	    /* Do not remove, I would like to test this further later...
	     * It did not converge at first attempt...
	     * Øystein Godøy, METNO/FOU, 12.04.2007 
	    p = -8.48841152
		+2.90687352*(cpar.A2/cpar.A1)
		-8.06381521*(cpar.A3/cpar.A1)
		+0.06169313*cpar.soz
		+0.01925053*cpar.saz;
	    x = -6.92968721
		+3.15424360*(cpar.A2/cpar.A1)
		-8.51300241*(cpar.A3/cpar.A1)
		+0.04520431*cpar.soz;
	    x = -7.147843141
		+0.006599748*(cpar.A1/cos(deg2rad(cpar.soz)))
		+2.962917376*(cpar.A2/cpar.A1)
		-8.593505367*(cpar.A3/cpar.A1)
		+0.047032443*(cpar.soz);
	    x = -5.012037873
		+0.008460015*(cpar.A1/cos(deg2rad(cpar.soz)))
		-6.927675661*(cpar.A3/cpar.A1)
		+0.044061292*(cpar.soz);
	    x = -6.57091311
		+0.00738986*(cpar.A1/cos(deg2rad(cpar.soz)))
		-5.72279895*(cpar.A3/cpar.A1)
		+0.05807489*(cpar.soz);

	    x = -4.761735
		+3.649293*(cpar.A2/cpar.A1)
		-7.423763*(cpar.A3/cpar.A1);
	    p = exp(x)/(1+exp(x));
	    */
	    
	    if (p.pice < 0.0) {
		class[i] = 0;
	    } else if (p.pice < 0.05) {
		class[i] = 1;
	    } else if (p.pice < 0.10) {
		class[i] = 2;
	    } else if (p.pice < 0.15) {
		class[i] = 3;
	    } else if (p.pice < 0.20) {
		class[i] = 4;
	    } else if (p.pice < 0.25) {
		class[i] = 5;
	    } else if (p.pice < 0.30) {
		class[i] = 6;
	    } else if (p.pice < 0.35) {
		class[i] = 7;
	    } else if (p.pice < 0.40) {
		class[i] = 8;
	    } else if (p.pice < 0.45) {
		class[i] = 9;
	    } else if (p.pice < 0.50) {
		class[i] = 10;
	    } else if (p.pice < 0.55) {
		class[i] = 11;
	    } else if (p.pice < 0.60) {
		class[i] = 12;
	    } else if (p.pice < 0.65) {
		class[i] = 13;
	    } else if (p.pice < 0.70) {
		class[i] = 14;
	    } else if (p.pice < 0.75) {
		class[i] = 15;
	    } else if (p.pice < 0.80) {
		class[i] = 16;
	    } else if (p.pice < 0.85) {
		class[i] = 17;
	    } else if (p.pice < 0.90) {
		class[i] = 18;
	    } else if (p.pice < 0.95) {
		class[i] = 19;
	    } else if (p.pice <= 1.0) {
		class[i] = 20;
	    } else {
		class[i] = 0;
	    }

	    if ((p.pice > p.pfree) && (p.pice > p.pcloud)) {
	      cat[i] = ICE;
	    } else if ((p.pfree > p.pice) && (p.pfree > p.pcloud)) {
	      cat[i] = CLEAR;
	    } else if ((p.pcloud > p.pice) && (p.pcloud > p.pfree)){
	      cat[i] = CLOUD;
	    } else { /*some probs. are equal*/
	      cat[i] = UNCL;
	    }

	}
    }
    fmlogmsg(where,"Now returning to main...");
   
    return(FM_OK);
}

