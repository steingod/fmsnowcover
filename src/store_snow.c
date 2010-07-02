/*
 * NAME:
 * store_snow
 *
 * PURPOSE:
 * Writes image data on TIFF formatted file, ready for visualization
 * on any standard image viewer. 
 * 
 * REQUIRES:
 * libfmio
 * libtiff
 * 
 * AUTHOR: 
 * Øystein Godøy, DNMI/FOU, 01/08/1995
 * MODIFIED:
 * Øystein Godøy, DNMI/FOU, 22/06/2000
 * Mari Anne Killie, DNMI/FOU, 03/04/2009
 * Mari Anne Killie, DNMI/FOU, 02/07/2010 Introduced image_type to
 * select between different types. 0: probability for snow/ice (20
 * classes), 1: categorized image (5 classes), 2: SAR testing
 *
 * CVS_ID:
 * $Id: store_snow.c,v 1.3 2010-07-02 15:10:27 mariak Exp $
 */

#include <fmio.h>
#include <tiffio.h>
#include <fmaccusnow.h>

int store_snow(char *fname,unsigned char *im,fmio_mihead clinfo,int image_type){

  int i, ret, numcat;
  uint16 cm[3][256];
  char *strclali, *info;   
  char *par="PIXEL CLASSIFICATION\n";
  char date[17], satinfo[FMIO_TIFFHEAD];
  
  if (image_type == 0) {
    numcat = PROBLIMITS;
  } else if (image_type == 1) {
    numcat = CATLIMITS;
  } else if (image_type == 2) {
    numcat = 2;
  } else { /*invalid value for image_type, exit!*/
    fprintf(stderr," ERROR(store_snow): ");
    fprintf(stderr,"invalid value for image_type (%d), exiting.\n",image_type);
    return(FM_IO_ERR);
  }
  
  char *desc[numcat];
  for(i=0; i<numcat; ++i ){
    desc[i] = (char*)malloc(20*sizeof(char)); 
  }
  
  if (image_type == 0) {
    strcpy(desc[0],"[0.00 -  0.05>");
    strcpy(desc[1],"[0.05 -  0.10>"); 
    strcpy(desc[2],"[0.10 -  0.15>");
    strcpy(desc[3],"[0.15 -  0.20>");
    strcpy(desc[4],"[0.20 -  0.25>");
    strcpy(desc[5],"[0.25 -  0.30>");
    strcpy(desc[6],"[0.30 -  0.35>");
    strcpy(desc[7],"[0.35 -  0.40>");
    strcpy(desc[8],"[0.40 -  0.45>");
    strcpy(desc[9],"[0.45 -  0.50>"); 
    strcpy(desc[10],"[0.50 -  0.55>");
    strcpy(desc[11],"[0.55 -  0.60>"); 
    strcpy(desc[12],"[0.60 -  0.65>");
    strcpy(desc[13],"[0.65 -  0.70>"); 
    strcpy(desc[14],"[0.70 -  0.75>");
    strcpy(desc[15],"[0.75 -  0.80>"); 
    strcpy(desc[16],"[0.80 -  0.85>");
    strcpy(desc[17],"[0.85 -  0.90>"); 
    strcpy(desc[18],"[0.90 -  0.95>");
    strcpy(desc[19],"[0.95 -  1.00]");
  } else if (image_type == 1) {
    strcpy(desc[0],"1: Ice/snow");
    strcpy(desc[1],"2: Clear");
    strcpy(desc[2],"3: Clouded");
    strcpy(desc[3],"4: Unclass");
    strcpy(desc[4],"5: Undef");
  } else if (image_type == 2) {
    strcpy(desc[0],"0: no SAR update");
    strcpy(desc[1],"1: SAR update");
  }
  
  /*Create string containing AVHRR image characteristics.*/
  sprintf(date, "%02d:%02d %02d/%02d-%4d", clinfo.hour, clinfo.minute, 
	  clinfo.day, clinfo.month, clinfo.year);
  fm_MITIFF_create_head(satinfo, clinfo.satellite, date, 0, 1, "C", 
			clinfo.xsize, clinfo.ysize,"Polar Stereographic", 
			"60 N", 0., 1000., 1000., 0., 0., 
			clinfo.Ax, clinfo.Ay, clinfo.Bx, clinfo.By, "");   
  
  strclali = (char *) malloc(DUMMYSTR);
  if (!strclali) {
    fprintf(stderr," ERROR(store_snow): ");
    fprintf(stderr,"Could not allocate memory\n");
    return(FM_MEMALL_ERR);
  }
 
  info = (char *) malloc(DUMMYSTR+numcat*CLASSLIMITSSTR+strlen(par));
  if (!info) {
    fprintf(stderr," ERROR(snow_mitiff_product): ");
    fprintf(stderr,"Could not allocate memory\n");
    return(FM_MEMALL_ERR);
  }

  sprintf(info, "\n COLOR INFO:\n %s",par);
  sprintf(strclali, " %d\n", numcat);
  strcat(info,strclali);
  for (i=0; i<numcat; i++) {
    sprintf(strclali, " %s\n", desc[i]);
    strcat(info, strclali);
  }
  strcat(satinfo, info);

  /*Probability image:*/
  if (image_type == 0) {
    fmheatmap(numcat,cm[0],cm[1],cm[2]);
  }
  else {
    /*Classified (category) image:*/
    for (i=0; i<256; i++) {
      cm[0][i] = cm[1][i] = cm[2][i] = 0;
    }
    cm[0][C_ICE]=cm[1][C_ICE]=cm[2][C_ICE]=255*255;/*snow/ice: white*/
    cm[0][C_CLEAR]  = 0;/*clear: blue (in Diana: cut this colour + set backgr!*/
    cm[1][C_CLEAR]  = 151*255;
    cm[2][C_CLEAR]  = 255*255;
    cm[0][C_CLOUDED]= 99*255;/*clouds: grey*/
    cm[1][C_CLOUDED]= 86*255;
    cm[2][C_CLOUDED]= 82*255;
    cm[0][C_UNCLASS]=cm[1][C_UNCLASS]=cm[2][C_UNCLASS]=151*255;/*uncl.: grey*/
    cm[0][C_UNDEF]=cm[1][C_UNDEF]=cm[2][C_UNDEF]= 0;/*udefined: black*/
  }
  /*Later: add option for when image_type == 1 and image_type == 2 separately*/

  ret = fm_MITIFF_write_imagepal(fname, im, satinfo, clinfo, cm); 
  
  free(strclali);
  free(info);
  
  return(FM_OK);
}

