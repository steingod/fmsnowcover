/******************************************************************
 * COPYRIGHT: EUMETSAT
 *
 * PRODUCED BY:
 * Norwegian Meteorological Institute (DNMI)
 * Research and Development Department
 * P.O.BOX 43 - Blindern, N-0313 OSLO, NORWAY
 *       
 * This SW was developed by DNMI and DMI within the context of the
 * Co-operation Agreement for the development of a pilot SAF on
 * Ocean and Sea Ice.
 *****************************************************************/
 
/*
 * NAME: accuice.h
 *
 * PURPOSE: Header file for High Latitude O&SI SAF SST
 *          software.
 *
 * AUTHOR: Steinar Eastwood, DNMI/FOU, 21.08.2000
 * MODIFIED: 
 * SE, DNMI, 23.01.2001, 11.06.2001
 * SE, DNMI, 11.09.2001    Changed names on BIAS variables
 * SE, DNMI, 17.12.2001    NUMSAT = 5
 * SE, met.no, 10.01.2003  New function interface.
 * SE, met.no, 09.07.2003  New Time SST HDF variables.
 * SE, met.no, 28.11.2003  New name parameters and function prototypes.
 * SE, met.no, 02.12.2003  New function prototypes, included time_conv.h.
 * SE, met.no, 18.02.2004  LENNMHRSST 34 -> MINLENNMHRSST 31
 * SE, met.no, 20.02.2004  Removed T4BIASN14 test for NOAA14. Using mean SST
 *                         clim in upper gross error check if indicated.
 * SE, met.no, 09.08.2004  New parameter in function calc_SST.
 * SE, met.no, 11.08.2004  New version of SST HDF5 file formats.
 * SE, met.no, 18.10.2004  NOSATPROD -> MAXSATPROD, now = 8.
 * SE, met.no, 23.02.2005  New funcion extrapol_SeaIceEdge.
 * SE, met.no, 20.04.2005  Introduced THRSZA, change in average_merge_filesQF.
 * SE, met.no, 29.04.2005  Introduced LML and NAR products on HDF format.
 * SE, met.no, 07.06.2005  Allowing both signed and unsigned short for time SST field.
 * SE, met.no, 26.02.2006  Added iceflag.
 * SE, met.no, 06.03.2006  Added cloudflag.
 *
 */

#ifndef _SAF_SST_H
#define _SAF_SST_H

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
#define BASENMAVHRRICE  "fmsnow_"     /* Base name of AVHRR ice files */
#define MINLENNMAVHRRICE   27      /* Minimum AVHRR ice file name length */

#define MINPROBAVHRR 0.
#define MAXPROBAVHRR 100.


/* Adopting these from avhrrice_pap */
#define OSIMISVAL_NOCOV -991 
#define OSIMISVAL_NIGHT -990 
#define OSIMISVAL_LAND -992
#define OSIMISVAL_3A -993


/* Parameters for HDF5 SST files */
#define HLSSTH5P_PROJSTR "+proj=stere +a=6371000 +lon_0=0 +lat_ts=60 +b=6371000 +lat_0=90"

#define CLASS_DT OSI_INT
#define PROB_DT OSI_FLOAT
#define PROB_ICE_DESC "P(ice)"
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
#define AVHRRICEPROD_LEVELS 3
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

#endif /* SAF_SST_H */
