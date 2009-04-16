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
 * Mari Anne Killie, DNMI/FOU, 03/04/2009
 *
 * CVS_ID:
 */

#include <fmio.h>
#include <tiffio.h>
#include <accusnow.h>

int store_snow(char *filename, PRODhead ph, unsigned char *im, int numcat, char *desc[], char *satstring) {

  int i, ret;
  char date[17], satinfo[FMIO_TIFFHEAD];
  char *strclali, *info;   
  char *par="PIXEL CLASSIFICATION\n";
  uint16 cm[3][256];
  /*unsigned short cmny[3][256];*/
  fmio_mihead clinfo = {
	"Not known",
	00, 00, 00, 00, 0000, -9, 
	{0, 0, 0, 0, 0, 0, 0, 0}, 
	0, 0, 0, 0., 0., -999., -999.
  };

  /*Create string containing AVHRR image characteristics.*/
  sprintf(date, "%02d:%02d %02d/%02d-%4d", ph.hour, ph.minute, ph.day, 
	  ph.month, ph.year);
  
  fm_MITIFF_create_head(satinfo, satstring, date, 0, 1, "C", ph.iw, ph.ih,
			"Polar Stereographic", "60 N", 0., 1000., 1000., 
			0., 0., ph.Ax, ph.Ay, ph.Bx, ph.By, "");   
  
  sprintf(clinfo.satellite,"%s",satstring);
  clinfo.hour = ph.hour;
  clinfo.minute = ph.minute;
  clinfo.day = ph.day;
  clinfo.month = ph.month;
  clinfo.year = ph.year;
  clinfo.zsize = 1;
  clinfo.xsize = ph.iw;
  clinfo.ysize = ph.ih;
  clinfo.Ax = ph.Ax;
  clinfo.Ay = ph.Ay;
  clinfo.Bx = ph.Bx;
  clinfo.By = ph.By;

  strclali = (char *) malloc(DUMMYSTR);
    if (!strclali) {
	fprintf(stderr," ERROR(snow_mitiff_product): ");
	fprintf(stderr,"Could not allocate memory\n");
	return(2);
    }
    
    info = (char *) malloc(DUMMYSTR+numcat*CLASSLIMITSSTR+strlen(par));
    if (!info) {
	fprintf(stderr," ERROR(snow_mitiff_product): ");
	fprintf(stderr,"Could not allocate memory\n");
	return(2);
    }

    sprintf(info, "\n COLOR INFO:\n %s",par);
    sprintf(strclali, " %d\n", numcat);
    strcat(info,strclali);
    for (i=0; i<numcat; i++) {
	sprintf(strclali, " %s\n", desc[i]);
	strcat(info, strclali);
    }
    strcat(satinfo, info);
     if (numcat == 20) {
      fmheatmap(numcat,cm[0],cm[1],cm[2]);
    }
    else {
      for (i=0; i<256; i++) {
	cm[0][i] = cm[1][i] = cm[2][i] = 0;
      }
      cm[0][1] = cm[1][1] = cm[2][1] = 255*255;/*snø/is - hvitt*/
      cm[0][3] = cm[1][3] = cm[2][3] = 151*255;/*uklass - grått*/
      cm[0][2] = 99*255;/*klart - grålilla*/
      cm[1][2] = 86*255;
      cm[2][2] = 82*255;
      cm[0][4] = 0;/*blått - skyer*/
      cm[1][4] = 151*255;
      cm[2][4] = 255*255;
      cm[0][5] = cm[1][5] = cm[2][5] = 0;/*udef - sort*/
    }

     ret = fm_MITIFF_write_imagepal(filename, im, satinfo, clinfo, cm); 
     
     free(strclali);
     free(info);
    
  return(0);
}

