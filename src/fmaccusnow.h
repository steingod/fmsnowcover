/*
 * NAME: 
 * fmaccusnow.h
 *
 * PURPOSE: 
 * Header file used within time integration of passage snow products.
 *
 * AUTHOR: 
 * Steinar Eastwood, DNMI/FOU, 21.08.2000
 *
 * MODIFIED: 
 * Øystein Godøy, METNO/FOU, 23.04.2009: Modified for use within the
 * fmsnowcover package.
 *
 * CVS_ID:
 * $Id: fmaccusnow.h,v 1.1 2009-04-23 10:43:43 steingod Exp $
 */

#ifndef _FMACCUSNOW_H
#define _FMACCUSNOW_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <safhdf.h>
#include <tiffio.h>
#include <fmutil.h>
#include <fmio.h>
#include <dirent.h>


/* File name specific parameters */
#define BASEFNAME  "fmsnow_"     /* Base name of AVHRR ice files */
#define MINLENFNAME   27      /* Minimum AVHRR ice file name length */

#define MINPROBAVHRR 0.
#define MAXPROBAVHRR 100.


/* Adopting these from avhrrice_pap */
#define FMACCUSNOWMISVAL_NOCOV -991 
#define FMACCUSNOWMISVAL_NIGHT -990 
#define FMACCUSNOWMISVAL_LAND -992
#define FMACCUSNOWMISVAL_3A -993


/* Parameters for HDF5 files */
#define ACCUSNOWH5P_PROJSTR "+proj=stere +a=6371000 +lon_0=0 +lat_ts=60 +b=6371000 +lat_0=90"

#define CLASS_DT OSI_INT
#define PROB_DT OSI_FLOAT
#define PROB_ICE_DESC "P(ice/snow)"
#define PROB_CLOUD_DESC "P(cloud)"
#define PROB_CLEAR_DESC "P(clear)"
#define C_ICE     1
#define C_CLEAR   2
#define C_CLOUDED 4
#define C_UNCLASS 3
#define C_UNDEF   5
#define PROB_MISVAL -199

#define PROBLIMITS 20
#define CLASSLIMITS 5
#define CLASSLIMITSSTR 66
#define DUMMYSTR 100
#define FILELEN 256        /* standard length of filenames including path */
#define FMACCUSNOWPROD_LEVELS 3
/*
 * Function prototypes.
 */


int average_merge_files(char **infSST, int nrInput, fmucsref safucs, 
                          unsigned char *class, unsigned char *probclass,
			  float *probice, float *probclear, float cloudlim);

int check_headers(int nrInput, PRODhead hrSSThead[]);

int check_sat_area(char **satlist, int numsat, char *filename);

int find_sat_area_index(char **satlist, int numsat, char *filename);

int read_sat_area_list(char *listfile, char **elemlist);

int store_snow(char *filename, PRODhead ph, unsigned char *im, 
			 int numcat, char *desc[], char *satstring);

void usage();

/*
 * End function prototypes.
 */

#endif /* _FMACCUSNOW_H */
