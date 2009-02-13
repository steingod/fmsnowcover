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
 * $Id: normalpdf.c,v 1.1 2009-02-13 23:23:14 steingod Exp $
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
