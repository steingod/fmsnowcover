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
 * NAME: func_accusnow.c (working name)
 * 
 * PURPOSE: Functions used by accusnow.
 * 
 *
 * AUTHOR: Steinar Eastwood, DNMI, 16.11.2000
 * MODIFIED: 
 * SE, DNMI, 23.01.2001
 * SE, met.no, 10.01.2003  New version of average_merge_filesQF.
 * SE, met.no, 10.01.2003  Error in average_merge_filesQF.
 * SE, met.no, 28.10.2003  Improved error handling when reading ice file.
 * SE, met.no, 28.11.2003  New time functions.
 * SE, met.no, 02.12.2003  Moved time functions to func_time_conv.c.
 * SE, met.no, 09.12.2003  Skipping HDF files if reading fails.
 *                         New merging func -> standard.
 * SE, met.no, 13.01.2004  exit replaced by return.
 * SE, met.no, 10.08.2004  HR SST files now contain SST, QFLG and SZA.
 * SE, met.no, 23.02.2005  Extrapolating to unclassified values in ice field.
 *                         Check for land before ice in average_merge_filesQF.
 * SE, met.no, 20.04.2005  Optional use of satellite senith angle in selecting
 *                         obs for averaging.
 * SE, met.no, 01.11.2006  Use dummy values for quality flags and sat zenith if 
 *                         old hdf format.
 *
 *MAK, met.no, 08.01.2009 Using func_average_SST.c as starting point
 *for files that will compile a snowmap for the cryorisk project.
 */ 

#include <accusnow.h>



/* 
 *  Function to loop through all input files checking each
 *  pixel. Pixels with probability of cloud larger than given
 *  'cloudlim' are thrown away. The remaining are summed, and then
 *  averaged to find a probability for snow for cloudfree case.
 *
 */

int average_merge_files(char **infAVHRRICE, int nrInput, fmucsref safucs, 
                          unsigned char *class, unsigned char *probclass, 
			  float *probice, float *probclear, float cloudlim)
{

  char *errmsg="\n\tERROR(average_merge_files): ";
  int i, elem, pn, status, ret, size_n;
  unsigned int xc, yc;
  int *numPix, *numIce, *numLand, *numCloud, *numCloudfree, *numUndef;
  float Pice_val, Pclear_val, Pcloud_val, probsum, sumCloudfree;
  float *sumIce, *sumClear;
  osihdf ice_h5p;

  /* Allocate memory */
  size_n = safucs.iw*safucs.ih;

  sumIce   = (float *) malloc(size_n*sizeof(float));
  sumClear = (float *) malloc(size_n*sizeof(float));
  numIce   = (int *) malloc(size_n*sizeof(int));
  numCloudfree = (int *) malloc(size_n*sizeof(int*));
  numPix   = (int *) malloc(size_n*sizeof(int));
  numLand  = (int *) malloc(size_n*sizeof(int));
  numCloud = (int *) malloc(size_n*sizeof(int));
  numUndef = (int *) malloc(size_n*sizeof(int));  

  if (!sumIce || !sumClear || !numIce || !numPix || !numLand
      || !numCloud || !numUndef || !numCloudfree) {
     fprintf(stderr," Could not allocate memory for data field\n");
     return(3);
  }

  /* Initialize */
  for (i=0;i<size_n;i++)  {
    sumIce[i]   = 0.0;
    sumClear[i] = 0.0;
    numIce[i]   = 0;
    numCloudfree[i] = 0;
    numPix[i]   = 0;
    numLand[i]  = 0;
    numCloud[i] = 0;
    numUndef[i] = 0;
  }


  for (pn=0;pn<nrInput;pn++)  {      /* Loop through all sat.passes */   
   
    init_osihdf(&ice_h5p);

    ret = read_hdf5_product(infAVHRRICE[pn],&ice_h5p,0); /*0:reads everything*/
    if (ret) {
      fprintf(stderr,
	      "%s, Trouble encountered when reading data file %s (%d).\n", 
	      errmsg, infAVHRRICE[pn],status);
      fprintf(stderr,"\t Skipping file.\n");
      continue;
    }
      

    for (yc=0;yc<ice_h5p.h.ih;yc++) {
    for (xc=0;xc<ice_h5p.h.iw;xc++) {

      elem = fmivec(xc, yc, ice_h5p.h.iw);

      Pice_val   = ((float *) ice_h5p.d[0].data)[elem];
      Pclear_val = ((float *) ice_h5p.d[1].data)[elem];
      Pcloud_val = ((float *) ice_h5p.d[2].data)[elem];


      /*1) check that pixel has prob.value */
      if ( (Pcloud_val>=MINPROBAVHRR) && (Pcloud_val<=MAXPROBAVHRR) && (Pclear_val>=MINPROBAVHRR) && (Pclear_val<=MAXPROBAVHRR) && (Pice_val>=MINPROBAVHRR) && (Pice_val<=MAXPROBAVHRR) ){
	
	/*2) check that prob.values sum to ~1*/
	probsum = Pcloud_val + Pclear_val + Pice_val;
	if (probsum > 1.05 || probsum < 0.95) { 
	/*this should never be true due to similar check in avhrrice_pap!*/
	/* printf("P(ice): %f, P(clear): %f, P(cloud): %f\n",
	   Pice_val, Pclear_val, Pcloud_val); 
	  fprintf(stderr,"The probability does not add up to 1 (%f),",probsum);
	  fprintf(stderr," check input file %s\n", infAVHRRICE[pn]);*/
	  continue;
	}

	/*3) check cloud probability -> if too high, throw away pixel*/
	if (Pcloud_val >= cloudlim) {
	  numCloud[elem] ++;
	  numPix[elem] ++;
	  continue;
	}
	
	/*4) compute a prob based on the ratio between clear and ice/snow*/
	else { 
 	  sumCloudfree = Pclear_val + Pice_val;
	  if (sumCloudfree <= MINPROBAVHRR) {
	    /* will not happen unless cloudlim > 0.95 (still unlikely)*/
	    fprintf(stderr,"Not nice to divide by zero, check cloudlim!\n");
	    return(8); /*random return value used.. */
	  }
	  numCloudfree[elem] ++; 
	  numPix[elem] ++;
	  sumIce[elem] += Pice_val/sumCloudfree;
	  sumClear[elem] += Pclear_val/sumCloudfree;
	}
      }

      /* if NOT prob.value for this pixel: */
      else if (Pice_val == OSIMISVAL_NOCOV || Pice_val == OSIMISVAL_NIGHT || Pice_val == OSIMISVAL_3A){ /* Undefined*/
	if (Pcloud_val != Pice_val || Pclear_val != Pice_val) {
	  /*not supposed to happen, check avhrrice_pap routines!*/
	  fprintf(stderr,
		  "Strange values encountered for pixel %d in file %s\n",
		  elem,infAVHRRICE[pn]);
	  fprintf(stderr,"(P(ice) = %f, P(clear) = %f, P(cloud) = %f)\n",
		  Pice_val,Pclear_val,Pcloud_val);
	  continue;
	}
	numUndef[elem] ++;	 
	numPix[elem] ++; 
      }
      
      else { /*also not supposed to happen, check avhrrice_pap/input files*/
	/*fprintf(stderr,"Invalid pixel values encountered in file %s\n",
		infAVHRRICE[pn]);
	fprintf(stderr,"\tP(ice) = %f, P(clear) = %f, P(cloud) = %f\n",
	Pice_val,Pclear_val,Pcloud_val);*/
	continue;
      }
   
    }
    } /*finished looping through all pixels for current sat.pass*/
 
    
    if (free_osihdf(&ice_h5p) != 0) {
      fprintf(stderr,"%s Could not free ice_h5p properly.",errmsg);
      return(3);
    }

  } /*finished looping through all sat.passes*/


  /* Loop through grid and calculate average probabilities */

  for (elem=0;elem<size_n;elem++) {

    /* First control that things add up*/
    if (numCloudfree[elem] + numCloud[elem] + numUndef[elem] != numPix[elem]) {
      fprintf(stderr,"Something is wrong, check this!\n");
      printf("Element: %d\n",elem);
      printf("cloudfree: %d, cloud: %d, undef: %d, numpix: %d\n",numCloudfree[elem],numCloud[elem],numUndef[elem],numPix[elem]);
      return(8); /*again random return value chosen*/
    }
     
    /* If pixel is cloudfree for at least one sat.pass: */
    if (numCloudfree[elem] > 0) { 
      probice[elem] = sumIce[elem]/numCloudfree[elem];
      probclear[elem] = sumClear[elem]/numCloudfree[elem];
      if (probice[elem] > probclear[elem]) {  /*snow/ice*/
	class[elem] = C_ICE;
      }
      else if (probice[elem] < probclear[elem]) { /*clear*/
	class[elem] = C_CLEAR;
      } 
      else { /* Ice and clear equally likely */
	class[elem] = C_UNCLASS;
      }
    }
    /* Alternatively the pixel is clouded or undef. for all sat.passes */
    else if (numCloudfree[elem] == 0 && numCloud[elem] > 0) { 
      class[elem] = C_CLOUDED;
    }
    else { /* numPix == numUndef, No sat. data */
      if (numPix[elem] != numUndef[elem]) { /*unødv. kontroll*/
	fprintf(stderr,"Something wrong during pixel classification\n");
	return(8);
      }
      class[elem] = C_UNCLASS;
    }
    


    if (probice[elem] < 0.0) {
      probclass[elem] = 0;
    } else if (probice[elem] < 0.05) {
      probclass[elem] = 1;
    } else if (probice[elem] < 0.10) {
      probclass[elem] = 2;
    } else if (probice[elem] < 0.15) {
      probclass[elem] = 3;
    } else if (probice[elem] < 0.20) {
      probclass[elem] = 4;
    } else if (probice[elem] < 0.25) {
      probclass[elem] = 5;
    } else if (probice[elem] < 0.30) {
      probclass[elem] = 6;
    } else if (probice[elem] < 0.35) {
      probclass[elem] = 7;
    } else if (probice[elem] < 0.40) {
      probclass[elem] = 8;
    } else if (probice[elem] < 0.45) {
      probclass[elem] = 9;
    } else if (probice[elem] < 0.50) {
      probclass[elem] = 10;
    } else if (probice[elem] < 0.55) {
      probclass[elem] = 11;
    } else if (probice[elem] < 0.60) {
      probclass[elem] = 12;
    } else if (probice[elem] < 0.65) {
      probclass[elem] = 13;
    } else if (probice[elem] < 0.70) {
      probclass[elem] = 14;
    } else if (probice[elem] < 0.75) {
      probclass[elem] = 15;
    } else if (probice[elem] < 0.80) {
      probclass[elem] = 16;
    } else if (probice[elem] < 0.85) {
      probclass[elem] = 17;
    } else if (probice[elem] < 0.90) {
      probclass[elem] = 18;
    } else if (probice[elem] < 0.95) {
      probclass[elem] = 19;
    } else if (probice[elem] <= 1.0) {
      probclass[elem] = 20;
    } else {
      probclass[elem] = 0;
    }
    
  }


  free(sumIce);
  free(sumClear);
  free(numIce);
  free(numPix);
  free(numLand);
  free(numCloud);
  free(numUndef);
  free(numCloudfree);

  return(0);
}



/*
 *  Function to intercompare the headers of the SST products that are
 *  supposed to be from the same satellite pass.
 */

int check_headers(int nrInput, PRODhead hrSSThead[])
{

  int i, icheck;

  icheck = 0;

  for (i=0;i<nrInput;i++)  {
    if (hrSSThead[i].year != hrSSThead[0].year ||
        hrSSThead[i].month != hrSSThead[0].month ||
        hrSSThead[i].day != hrSSThead[0].day ||
        hrSSThead[i].hour != hrSSThead[0].hour ||
        hrSSThead[i].minute != hrSSThead[0].minute ||
        strcmp(hrSSThead[i].source,hrSSThead[0].source) != 0 ||
        strcmp(hrSSThead[i].product,hrSSThead[0].product) != 0)  {
      icheck = 10;
    }
  }
  return(icheck);
}




/*
 *  Function to read files with list of satellites or area tiles to 
 *  be used.
 *
 *  Return values:
 *  -2 : File contains no data (no lines).
 *  -1 : Can't find file.
 *   0 -> : Number of elements found on file.
 *
 *  Implemented 04.02.2005
 */


int read_sat_area_list(char *listfile, char **elemlist) 
{
  char *errmsg="\n\tERROR(read_sat_area_list): ";
  int i, numelem;
  char linein[100];
  FILE *fpin;

  numelem = i = 0;

  fpin = fopen(listfile,"r");
  if (!fpin) {
    fprintf(stderr,"%s Could not open file %s\n",errmsg,listfile);
    return(-1);
  }

  while (fgets(linein,99,fpin) != NULL) {
    i ++;
    if (linein[0] != '#' && strlen(linein) > 1) {
      elemlist[numelem] = (char *) malloc((strlen(linein)+1)*sizeof(char));
      sscanf(linein,"%s",elemlist[numelem]);
      numelem ++;
    }
  }
  
  if (i == 0) {
    fprintf(stderr,"%s File contains no data: %s\n",errmsg,listfile);
    return(-2);
  }
  
  return(numelem);
}



/*
 *  Function to check filename for matching satellite or area tile from
 *  list of approved values read from file.
 *
 *  The elements in elemlist are tested for occurence in filename on the
 *  form  "_string_", "_string." or "string_".
 *
 *  Return values:
 *  0 : Could not find any elements from list in file name.
 *  1 : Found one of the elements in file name.
 *
 *  Implemented 04.02.2005
 */

int check_sat_area(char **elemlist, int numelem, char *filename)
{

  int i;
  char tmp1[32], tmp2[32], tmp3[32];

  for (i=0;i<numelem;i++) {
    sprintf(tmp1,"_%s_",elemlist[i]);
    sprintf(tmp2,"_%s.",elemlist[i]);
    sprintf(tmp3,"%s_",elemlist[i]);

    if (strstr(filename,tmp1) != NULL || strstr(filename,tmp2) != NULL || 
	(strstr(filename,tmp3) != NULL && filename[0]==tmp3[0])) {
      return(1);
    }
  }

  return(0);
}

/*identical to above, but returning elemlist index!*/
int find_sat_area_index(char **elemlist, int numelem, char *filename)
{

  int i;
  char tmp1[32], tmp2[32], tmp3[32];

  for (i=0;i<numelem;i++) {
    sprintf(tmp1,"_%s_",elemlist[i]);
    sprintf(tmp2,"_%s.",elemlist[i]);
    sprintf(tmp3,"%s_",elemlist[i]);

    if (strstr(filename,tmp1) != NULL || strstr(filename,tmp2) != NULL || 
	(strstr(filename,tmp3) != NULL && filename[0]==tmp3[0])) {
      return(i);
    }
  }

  return(-1);
}

