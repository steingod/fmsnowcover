/*
 * NAME:
 * nwp_read.h
 *
 * PURPOSE:
 * 
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
 * Øystein Godøy, met.no/FOU, 18.10.2004 
 *
 * MODIFIED:
 * NA
 *
 * CVS_ID:
 * $Id: getnwp.h,v 1.3 2009-03-30 13:42:53 steingod Exp $
 */
#ifndef NWPICE_READ
#define NWPICE_READ

#include <stdio.h>
#include <stdlib.h>
#include <fmio.h>
#include <fmutil.h>

/* #include <satimg.h> */
#include <mifield.h>

#define FFNS 4
#define NOFIELDS1 1
#define NOFIELDS2 0
#define NOFIELDS NOFIELDS1+NOFIELDS2

typedef struct {
    fmsec1970 validtime;
    int leadtime;
    fmucsref refucs;
    float *t950hpa; /* temp at 950 hPa */
    float *t800hpa; /* temp at 800 hPa */
    float *t700hpa; /* temp at 700 hPa */
    float *t500hpa; /* temp at 500 hPa */
    float *t0m; /* surface temp */
    float *t2m; /* 2m temp */
    float *ps; /* mean sea level pressure */
    float *topo; /* model topography */
    float *pw; /* precipitable water */
    float *rh; /* relative humidity at surface */
} nwpice;

int nwpice_init(nwpice *nwp); 
int nwpice_read(char *fpath, char **filenames, int nrf, int nruns,
	fmtime reqtime, fmucsref refucs, nwpice *nwp);
int nwpice_free(nwpice *nwp);

#endif /* NWP_READ */
