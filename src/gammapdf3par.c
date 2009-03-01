/*
 * NAME:
 * gammapdf3par
 *
 * PURPOSE:
 * To estimate the PDF for satellite observations using a three parameter
 * Gamma distribution. 
 *
 * NOTES:
 * The Gamma distribution is only defined for positive values.
 * The Gamma function may create overflow for large values.
 *
 * RETURN VALUES:
 * The probability (positive value) is returned unless an error occured.
 * A negative value is returned if a negative argument or a too large
 * alpha value is submitted.
 *
 * REQUIRES:
 * tgamma function of the libc math implementation...
 * 
 * AUTHOR:
 * Øystein Godøy, met.no/FOU, 28.09.2004 
 *
 * MODIFIED:
 * NA
 *
 * CVS_ID:
 * $Id: gammapdf3par.c,v 1.2 2009-03-01 21:24:15 steingod Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double gammapdf3par(double alpha, double beta, double gamma, double x) {
    double g, lg;
    double gpdf;
    char *where="gammapdf3par";

    if (x <= 0) {
	errmsg(where,
	    "The Gamma distribution is only defined for positive values.");
	return(-1);
    }
    if (alpha > 170) {
	errmsg(where,
	    "The alpha parameter is too large.");
	return(-2);
    }

    /*
    lg = gamma(alpha);
    g = signgam*exp(lg);
    */

    gpdf = (pow(x-gamma,(alpha-1))*exp(-(x-gamma)/beta))/
	(pow(beta,alpha)*tgamma(alpha));

    return(gpdf);
}

