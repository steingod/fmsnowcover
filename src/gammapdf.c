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
 * gammapdf
 *
 * PURPOSE:
 * To estimate the PDF for satellite observations using a two parameter Gamma
 * distribution. 
 *
 * NOTES:
 * The Gamma distribution is only defined for positive values.
 * The Gamma function may create overflow for large values.
 * The Gamma function is not part of the ANSI C standard. In this function the
 * SGI mathematical library is used.
 *
 * RETURN VALUES:
 * The probability (positive value) is returned unless an error occured.
 * A negative value is returned if a negative argument or a too large
 * alpha value is submitted.
 *
 * AUTHOR:
 * �ystein God�y, DNMI/FOU, 10/05/1999
 *
 * MODIFIED:
 * �ystein God�y, met.no/FOU, 28.09.2004
 *
 * CVS_ID:
 * $Id: gammapdf.c,v 1.1 2009-02-13 23:23:13 steingod Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #include <satimg.h> */
#include <fmutil.h>

double gammapdf(double alpha, double beta, double x) {
    double g, lg;
    double gpdf;
    char *where="gammapdf";

    if (x <= 0) {
	fmerrmsg(where,
	    "The Gamma distribution is only defined for positive values.");
	return(-1);
    }
    if (alpha > 170) {
	fmerrmsg(where,
	    "The alpha parameter is too large.");
	return(-2);
    }

    lg = gamma(alpha);
    g = signgam*exp(lg);

    gpdf = (pow(x,(alpha-1))*exp(-x/beta))/(g*pow(beta,alpha));

    return(gpdf);
}

