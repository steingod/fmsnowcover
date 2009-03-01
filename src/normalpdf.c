/*
 * NAME:
 * normalpdf
 *
 * PURPOSE:
 * To estimate probability using a normal distribution.
 *
 * REQUIREMENTS:
 *
 * INPUT:
 *
 * OUTPUT:
 *
 * NOTES:
 *
 * BUGS:
 *
 * AUTHOR:
 * Øystein Godøy, met.no/FOU, 18.11.2004 
 *
 * MODIFIED:
 * NA
 *
 * CVS_ID:
 * $Id: normalpdf.c,v 1.2 2009-03-01 21:24:15 steingod Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #include <satimg.h> */
#include <fmutil.h>

double normalpdf(double mean, double sdev, double x) {

    /*char *where="normalpdf";*/
    double npdf;

    npdf = (1./(sdev*sqrt(2.*fmPI)))*exp(-pow((x-mean),2.)/(2.*pow(sdev,2.)));

    return(npdf);
}
