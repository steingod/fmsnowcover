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
 * $Id: fmaccusnow.h,v 1.5 2013-02-01 10:31:28 steingod Exp $
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

/*
 * Useful constants
 */
#define LISTLEN 100
#define DEFAULTCLOUD 0.2
#define TOTAREAS 8
#define MAXSAT 20
#define MAXAREA 20
#define DATESTRINGLENGTH 50 /* Holds even ISO strings */


/* 
 * File name specific parameters 
 */
#define BASEFNAME "fmsnow_" /* Base name of AVHRR passage ice files */
#define MINLENFNAME 27 /* Minimum AVHRR ice file name length */

#define MINPROBAVHRR 0.
#define MAXPROBAVHRR 100.

/* 
 * Codes to be used within the output file
 */
#define FMACCUSNOWMISVAL_NOCOV -991  /* No coverage */
#define FMACCUSNOWMISVAL_NIGHT -990  /* Night scene */
#define FMACCUSNOWMISVAL_LAND -992 /* Land pixel */
#define FMACCUSNOWMISVAL_3A -993 /* AVHRR 3A missing */

/* 
 * Parameters for HDF5 files
 */
#define ACCUSNOWH5P_PROJSTR "+proj=stere +a=6371000 +lon_0=0 +lat_ts=60 +b=6371000 +lat_0=90"

#define CLASS_DT OSI_INT
#define PROB_DT OSI_FLOAT
#define PROB_ICE_DESC "P(ice/snow)"
#define PROB_CLOUD_DESC "P(cloud)"
#define PROB_CLEAR_DESC "P(clear)"
#define C_ICE     1
#define C_CLEAR   2
#define C_CLOUDED 3
#define C_UNCLASS 4
#define C_UNDEF   5
#define PROB_MISVAL -199

#define PROBLIMITS 20
#define CATLIMITS 5
#define CLASSLIMITSSTR 66
#define DUMMYSTR 100
#define FILELEN 256 /* standard length of filenames including path */
#define FMACCUSNOWPROD_LEVELS 3

/*
 * Function prototypes.
 */
int average_merge_files(char **infSST, int nrInput, fmucsref safucs, 
			unsigned char *class, unsigned char *probclass,
			float *probice, float *probclear, float cloudlim,
			int *numCloudfree);

int check_headers(int nrInput, PRODhead hrSSThead[]);

int check_sat_area(char **satlist, int numsat, char *filename);

int find_sat_area_index(char **satlist, int numsat, char *filename);

int read_sat_area_list(char *listfile, char **elemlist);

int store_snow(char *fname,unsigned char *im,fmio_mihead clinfo,int image_type);

void usage();

/*
 * End function prototypes.
 */

#endif /* _FMACCUSNOW_H */
