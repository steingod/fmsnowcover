/*
 * PURPOSE:
 * Writes image data on TIFF formatted file, ready for visualization
 * on any standard image viewer. TIFF tag number 262 is photometric
 * interpretation. This tag is 1 for grayscale images and 3 for
 * palette images. In situations when grayscale images are wanted
 * TIFF tag number 320 can be commented out.
 * 
 * REQUIRES:
 * libfmio
 * libtiff
 * 
 * AUTHOR: 
 * Øystein Godøy, DNMI/FOU, 01/08/1995
 * MODIFIED:
 * Øystein Godøy, DNMI/FOU, 22/06/2000
 *
 * CVS_ID:
 * $Id: store_mitiff_product.c,v 1.2 2012-01-04 11:40:03 mariak Exp $
 */

#include <fmio.h>
#include <tiffio.h>
#include <fmsnowcover.h>

#define CLASSLIMITS 20
#define CLASSLIMITSSTR 66
#define DUMMYSTR 100

void store_mitiff_result(char *outfile, unsigned char *icep, fmio_mihead img) {
    
    int i;
    uint16 cm[3][256];
    char *info, *strclali;
    char *par="PROBABILITY OF SNOW\n";
    char satinfo[FMIO_TIFFHEAD], date[17];
    char *clali[CLASSLIMITS] = {
	"[0.00 -  0.05>", 
	"[0.05 -  0.10>", 
	"[0.10 -  0.15>", 
	"[0.15 -  0.20>", 
	"[0.20 -  0.25>", 
	"[0.25 -  0.30>", 
	"[0.30 -  0.35>", 
	"[0.35 -  0.40>", 
	"[0.40 -  0.45>", 
	"[0.45 -  0.50>", 
	"[0.50 -  0.55>", 
	"[0.55 -  0.60>", 
	"[0.60 -  0.65>", 
	"[0.65 -  0.70>", 
	"[0.70 -  0.75>", 
	"[0.75 -  0.80>", 
	"[0.80 -  0.85>", 
	"[0.85 -  0.90>", 
	"[0.90 -  0.95>", 
	"[0.95 -  1.00]"
    };

    /*
     * Create string containing AVHRR image characteristics.
     */
    sprintf(date, "%02d:%02d %02d/%02d-%4d", img.hour, 
	img.minute, img.day, img.month, img.year);
    fm_MITIFF_create_head(satinfo, img.satellite, date, 0,
	1, "C", img.xsize, img.ysize,
	"Polar Stereographic", "60 N", 0.,
	1000., 1000., 0., 0., 
	img.Ax, img.Ay, img.Bx, img.By, 
	"");     

    /*
     * Change satellite information to correspond with the present
     * classified image.
     */
    strclali = (char *) malloc(DUMMYSTR);
    if (!strclali) {
	fprintf(stderr," ERROR(store_mitiff_result): ");
	fprintf(stderr,"Could not allocate memory\n");
	exit(0);
    }
    info = (char *) malloc(DUMMYSTR+CLASSLIMITS*CLASSLIMITSSTR+strlen(par));
    if (!info) {
	fprintf(stderr," ERROR(store_mitiff_result): ");
	fprintf(stderr,"Could not allocate memory\n");
	exit(0);
    }
    sprintf(info, "\n COLOR INFO:\n %s", par);
    sprintf(strclali, " %d\n", CLASSLIMITS);
    strcat(info, strclali);
    for (i=0; i<CLASSLIMITS; i++) {
	sprintf(strclali, " %s\n", clali[i]);
	strcat(info, strclali);
    }
    strcat(satinfo, info);

    /*
     * Actual writing of image.
     */
    fmheatmap(CLASSLIMITS,cm[0],cm[1],cm[2]);
    fm_MITIFF_write_imagepal(outfile, icep, satinfo, img, cm);

    free(strclali);
    free(info);
}


void store_mitiff_cat(char *outfile, unsigned char *cat, fmio_mihead img) {
    
    int i;
    uint16 cm[3][256];
    char *info, *strclali;
    char *par="PIXEL CLASSIFICATION\n";
    char satinfo[FMIO_TIFFHEAD], date[17];
    char *clali[CATLIMITS] = {
	"1: Ice/snow", 
	"2: Clear", 
	"3: Cloud", 
	"4: Unclass",
	"5: Undef"
    };

    /*
     * Create string containing AVHRR image characteristics.
     */
    sprintf(date, "%02d:%02d %02d/%02d-%4d", img.hour, 
	img.minute, img.day, img.month, img.year);
    fm_MITIFF_create_head(satinfo, img.satellite, date, 0,
	1, "C", img.xsize, img.ysize,
	"Polar Stereographic", "60 N", 0.,
	1000., 1000., 0., 0., 
	img.Ax, img.Ay, img.Bx, img.By, 
	"");     

    /*
     * Change satellite information to correspond with the present
     * classified image.
     */
    strclali = (char *) malloc(DUMMYSTR);
    if (!strclali) {
	fprintf(stderr," ERROR(store_mitiff_result): ");
	fprintf(stderr,"Could not allocate memory\n");
	exit(0);
    }
    info = (char *) malloc(DUMMYSTR+CATLIMITS*CLASSLIMITSSTR+strlen(par));
    if (!info) {
	fprintf(stderr," ERROR(store_mitiff_result): ");
	fprintf(stderr,"Could not allocate memory\n");
	exit(0);
    }
    sprintf(info, "\n COLOR INFO:\n %s", par);
    sprintf(strclali, " %d\n", CATLIMITS);
    strcat(info, strclali);
    for (i=0; i<CATLIMITS; i++) {
	sprintf(strclali, " %s\n", clali[i]);
	strcat(info, strclali);
    }
    strcat(satinfo, info);

    /*
     * Actual writing of image.
     */
    /*fmheatmap(CATLIMITS,cm[0],cm[1],cm[2]);*/
    for (i=0; i<256; i++) {
	cm[0][i] = cm[1][i] = cm[2][i] = 0;
    }
    cm[0][ICE] = cm[1][ICE] = cm[2][ICE] = 255*255;/*snow/ice - white*/
    cm[0][CLEAR] = 0;/*clear - blue (in Diana: cut this colour + set backgr!)*/
    cm[1][CLEAR] = 151*255;
    cm[2][CLEAR] = 255*255;
    cm[0][CLOUD] = 99*255;/*clouds - grey*/
    cm[1][CLOUD] = 86*255;
    cm[2][CLOUD] = 82*255;
    cm[0][UNCL]  = cm[1][UNCL] = cm[2][UNCL] = 151*255;/*unclass. - grey*/
    cm[0][UNDEF] = cm[1][UNDEF] = cm[2][UNDEF] = 0;/*udefined - black*/
    
    fm_MITIFF_write_imagepal(outfile, cat, satinfo, img, cm);

    free(strclali);
    free(info);
}
