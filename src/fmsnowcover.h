/*
 * NAME:
 * fmsnowcover.h
 *
 * PURPOSE:
 * NA
 *
 * REQUIREMENTS:
 * NA
 *
 * INPUT:
 * NA
 *
 * OUTPUT:
 * NA
 *
 * NOTES:
 * NOTES:
 * For further information, check the main program file pap_ice.c
 * A description of program features are given in the header.
 *
 * Further chacking of string lengths should be imposed, now possible
 * memory leaks can be caused.
 *
 * BUGS:
 * NA
 *
 * AUTHOR:
 * Øystein Godøy, DNMI/FoU, 26/11/2000.
 *
 * MODIFIED:
 * Øystein Godøy, met.no/FOU, 27.09.2004
 * Modification for full Bayes approach started...
 * Øystein Godøy, met.no/FOU, 03.11.2004: See pix_proc.c
 * Øystein Godøy, METNO/FOU, 30.10.2006: See pap_avhrrice.c
 * Øystein Godøy, METNO/FOU, 11.01.2007: Changed some names, added missing
 * 3A code
 * Øystein Godøy, METNO/FOU, 02.04.2007: Added 3B support.
 * Mari Anne Killie, METNO/FOU, 31.01.2008: See avhrrice_pap.c
 * Mari Anne Killie, METNO/FOU, 26.08.2008: Added A3b in struct
 * pinpstr and edited for r3a1/r3b1 in struct surfstr
 * Mari Anne Killie, METNO/FOU, 08.05.2009: snow added, d34 removed.
 *
 * CVS_ID:
 * $Id: fmsnowcover.h,v 1.13 2012-01-04 11:37:07 mariak Exp $
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <safhdf.h>
#include <projects.h>
#include <fmutil.h>
#include <fmio.h>
#include <getnwp.h>

#define FMSNOWCOVER_MSGLENGTH 255 /* String length for system messages */
#define MAXCHANNELS 6
#define MAXIMGSIZE 1440000
#define CLASSLIMITS 20	    /* Maximum number of classes in image */
#define CLASSLIMITSSTR 66   /* Length of string classlimit */
#define DUMMYSTR 100
#define NWP_NOFIELDS 5
#define FMSNOWCOVER_OLEVELS 3 /* Number of output levels */
#define FMSNOWCOVERMISVAL_NOCOV -991 
#define FMSNOWCOVERMISVAL_NIGHT -990 
#define FMSNOWCOVERMISVAL_LAND -992
#define FMSNOWCOVERMISVAL_3A -993
#define FMSNOWSUNZEN 85.
#define FMSNOWSEA 0 
#define FMSNOWLAND 191 /*works better than 255?!*/
/*The following 5 can be removed:*/
#define ICE 1
#define CLEAR 2
#define CLOUD 3
#define UNCL 4
#define UNDEF 5
#define CATLIMITS 5

/*
 * Some useful data constants to use in the software.
 */
/* PI should be defined by the standard C library on the UNIX platform,
 * define if not already existing... */
#ifndef PI
#define PI 3.141592654		/* Mathematical constant PI */
#endif
#define DEG2RAD PI/180.		/* Factor to multiply with to get radians */
#define RAD2DEG 180./PI		/* Factor to multiply with to get degrees */
#define FILENAME 50
#define FILELEN 256		/* Standard length of filenames incl path */
#define ANGBOX 10 	 	/* Box size in pixels for viewing geomtry */

/*
 * Data structure to hold configuration info etc.
 */
typedef struct {
    char imgpath[FILELEN];
    char nwppath[FILELEN];
    char cmpath[FILELEN];
    char lmpath[FILELEN];
    char productpath[FILELEN];
    char probtabname[FILELEN];
    char indexfile[FILELEN];
} cfgstruct;

/*
 * Data structure to hold time identification of satellite scene or equivalent.
 */
typedef struct {
    uint year;
    ushort month;
    ushort day;
    ushort hour;
    ushort minute;
} timestrct;

/*
 * Data structure to use when processing the individual pixels using the
 * Bayes approach...
 */
typedef struct {
    float A1;
    float A2;
    float A3;
    float A3b; /*The reflective part of channel 3b*/
    float T3;
    float T4;
    float T5;
    float soz;
    float saz;
    float tdiff;
    short lmask;
    short cmask;
    short algo;
    short daytime3b;
} pinpstr;

/*
 * Data structure to hold probaility estimates given miclpa
 */
typedef struct {
    double pice;
    double pfree;
    double pcloud;
} probstr;

/*
 * Data structure to hold probability coefficients read from file 
 */
typedef struct { /*struct for coefficients!*/
  char key;   /*'n' for normalpdf, 'g' for gamma, 't' for 3par-gamma*/
  double par1;   
  double par2;
  double par3;
  int count; /*counts the number of times coeffs are read!*/
} featstr;

typedef struct {
  featstr a1;
  featstr r21;
  featstr r3a1;
  featstr r3b1;
  featstr dt;
} surfstr;

typedef struct {
  surfstr ice;
  surfstr snow;
  surfstr cloud;
  surfstr water;
  surfstr land;
} statcoeffstr;

typedef struct {
  char feat[10];
  char surf[10];
  char key;
  double par1;
  double par2;
  double par3;
} dummystr;
  
/*
 * Prototypes
 */
void usage(void);

int decode_cfg(char cfgfile[],cfgstruct *cfg);

int process_pixels4ice(fmio_img img, 
    unsigned char *cmask[], unsigned char *lmask, nwpice nwp, 
    datafield *probs, unsigned char *class, unsigned char *cat,
    short algo, statcoeffstr cof);

void moment(float data[], int n, float *ave, float *adev, float *sdev,
    float *var, float *skew, float *curt);

int probest(pinpstr cpa, probstr *p, statcoeffstr cof);
double gammapdf(double alpha, double beta, double x);
double normalpdf(double mean, double sdev, double x);

/*void store_mitiff_result(char *outfile,unsigned char *icep,fmio_mihead img);*/
/*void store_mitiff_cat(char *outfile, unsigned char *cat, fmio_mihead img);*/
int rdstatcoeffs(char *coeffsfile, statcoeffstr *coeffs);
double findprob(featstr feat, double x, char *whereami);
int locstatcoeffs (dummystr dummies, statcoeffstr *cof);
int putcoeffs(featstr *feat, dummystr dummies);
float findcloudfree(datafield *d, int xsize, int ysize);
int updateindexfile(char *filename, char *avhrrfile, char *fmsnowfile,
    char *datetime, char *areaname, float validraw, float cloudfree); 
