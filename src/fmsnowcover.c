/*
 * PURPOSE:
 * To determine whether a pixel is covered by snow for AVHRR pixels. 
 * 
 * INPUT:
 * o Multichannel AVHRR satellite imagery in MEOS/HDF format, filename
 *   specified on commandline. Full path is required.
 * o Mask file in DNMI TIFF format specified by symbolic link. 
 *
 * OUTPUT:
 * o DNMI TIFF file containing colortable.
 *   Filename has to be specified at commandline. Full path is required.
 * 
 * BUGS:
 * Not thoroughly tested yet...
 *
 * Not all dynamically allocated memory is properly freed yet, code should be
 * checked... Especially code connected with use of strtok...
 *
 * AUTHOR: 
 * Øystein Godøy, DNMI/FOU, 13/12/2000
 *
 * MODIFICATIONS:
 * Øystein Godøy, DNMI/FOU, 02/03/2001
 * Datastream altered, removed some parts which were only used for testing.
 * Output format MITIFF was updated. Still a lot of cleaning of code is 
 * required.
 * Øystein Godøy, DNMI/FOU, 25/03/2001
 * Implemented use of configuration file, additional changes in user
 * interface and operation is required for easy porting to other sites.
 * Øystein Godøy, met.no/FOU, 27.09.2004: Modification for full Bayes
 * approach started...
 * Øystein Godøy, METNO/FOU, 30.10.2006: Adapting software to the results
 * achieved by Vibeke W. Thyness and Hanne Heiberg.  First the use of
 * NWCSAF PPS cloud mask etc is omitted, then add Bayes handling of
 * clouds, and increase number of output variables (add probability of
 * cloud and probability of open water), finally check if some other
 * statistics is better to use than the present one.
 * Øystein Godøy, METNO/FOU, 11.01.2007: Changed some names, and some
 * other modifications.
 * Mari Anne Killie, METNO/FOU, 31.01.2008: Some modifications in main
 * and decode_cfg + added functions rdstatcoeffs, locstatcoeffs and
 * putcoeffs, all so that the statistical coeffs needed in probest can
 * be read from text file rather than being hardcoded in probest.c.
 * Mari Anne Killie, METNO/FOU, 13.10.2008: introduced #ifdef
 * AVHRRICE_HAVE_NWP as quick way to comment out nwp.
 *
 * CVS_ID:
 * $Id: fmsnowcover.c,v 1.3 2009-03-10 13:22:23 mariak Exp $
 */
 
/* #include <satimg.h> */
#include <fmsnowcover.h>
#include <unistd.h>

int main(int argc, char *argv[]) {

    char *where="avhrrice_pap";
    char what[OSI_MSGLENGTH];
    extern char *optarg;
    int ret;
    short errflg = 0, iflg = 0, cflg = 0;
    short status;
    unsigned int size;
    char fname[FILENAME];
    /*char cfname[FILENAME];*/
    char pname[4];
    char *lmaskf, *opfn;
    char *infile, *cfgfile, *coffile;
    unsigned char *classed;
    cfgstruct cfg;
    FILE *lmask_located; /*Can be removed later*/
    fmio_mihead iinfo = {
	"Not known",
	00, 00, 00, 00, 0000, -9, 
	{0, 0, 0, 0, 0, 0, 0, 0}, 
	0, 0, 0, 0., 0., -999., -999.
    };
    /*
    struct mihead cminfo = {
	"Not known",
	00, 00, 00, 00, 0000, -9, 
	{0, 0, 0, 0, 0, 0, 0, 0}, 
	0, 0, 0, 0., 0., -999., -999.
    };
    */
    fmio_mihead clinfo = {
	"Not known",
	00, 00, 00, 00, 0000, -9, 
	{0, 0, 0, 0, 0, 0, 0, 0}, 
	0, 0, 0, 0., 0., -999., -999.
    };
    /*struct miheadpal cmpal;*/
    fmio_img img;
    fmucsref refucs;
    fmtime reftime;
    nwpice nwp;
    /*PRODhead header;*/
    osihdf lm;
    osihdf ice;
    osi_dtype ice_ft[OSI_OLEVELS]={OSI_FLOAT,OSI_FLOAT,OSI_FLOAT};
    char *ice_desc[OSI_OLEVELS]={"P(ice)","P(water)","P(cloud)"};

    statcoeffstr coeffs = {{{0}}};
    
    /*
     * Interprete commandline arguments.
     */
     while ((ret = getopt(argc, argv, "c:i:o:")) != EOF) {
	switch (ret) {
	    case 'c':
		cfgfile = (char *) malloc(FILELEN);
		if (!cfgfile) {
		    fmerrmsg(where,"Memory trouble.");
		    exit(FM_MEMALL_ERR);
		}
		if (!strcpy(cfgfile, optarg)) exit(0);
		cflg++;
                break;
	    case 'i':
		if (!strcpy(fname, optarg)) exit(0);
		iflg++;
                break;
	    default:
		usage();
	}
    }
    if (!iflg || !cflg) errflg++;
    if (errflg) usage();

    fprintf(stdout,"\n");
    fprintf(stdout," ================================================\n");
    fprintf(stdout," |                  AVHRRICE_PAP                |\n");
    fprintf(stdout," ================================================\n");
    fprintf(stdout,"\n");
 
    /*
     * Decode configuration file.
     */
    if (decode_cfg(cfgfile,&cfg) != 0) {
	fmerrmsg(where,"Could not decode configuration");
	exit(FM_IO_ERR);
    }

    /*
     * Set up datapaths etc.
     */
    infile = (char *) malloc(FILELEN);
    if (!infile) {
	fprintf(stderr,"%s\n"," Trouble processing:");
	fprintf(stderr,"%s\n",infile);
	fmerrmsg(where,"Could not allocate memory for infile");
	exit(FM_MEMALL_ERR);
    }
    sprintf(infile,"%s/%s",cfg.imgpath,fname);
    lmaskf = (char *) malloc(FILELEN);
    if (!lmaskf) {
	fprintf(stderr,"%s\n"," Trouble processing:");
	fprintf(stderr,"%s\n",infile);
	fmerrmsg(where,"Could not allocate memory for lmaskf");
	exit(FM_MEMALL_ERR);
    }
    if (strstr(fname,"ns") != NULL) {
	sprintf(pname,"%s","ns");
	sprintf(lmaskf,"%s/lmask_%s.hdf5",cfg.lmpath,"nse");
    } else if (strstr(fname,"at") != NULL) {
	sprintf(pname,"%s","at");
	sprintf(lmaskf,"%s/lmask_%s.hdf5",cfg.lmpath,"atl");
    } else if (strstr(fname,"nr") != NULL) {
	sprintf(pname,"%s","nr");
	sprintf(lmaskf,"%s/lmask_%s.hdf5",cfg.lmpath,"nrw");
    } else if (strstr(fname,"gr") != NULL) {
	sprintf(pname,"%s","gr");
	sprintf(lmaskf,"%s/lmask_%s.hdf5",cfg.lmpath,"grl");
    } else {
	fprintf(stderr,"%s\n"," Trouble processing:");
	fprintf(stderr,"%s\n",infile);
	fprintf(stderr," ERROR(main):  area not recognised\n");
	exit(FM_VAROUTOFSCOPE_ERR);
    }
    /*setting path to file containing probability coeffs*/
    coffile = (char *) malloc(FILELEN);
    if (!coffile) {
	fprintf(stderr,"%s\n"," Trouble processing:");
	fprintf(stderr,"%s\n",coffile);
	fmerrmsg(where,"Could not allocate memory for coffile");
	exit(FM_MEMALL_ERR);
    }
    sprintf(coffile,"%s",cfg.probtabname);
   
    /*
     * Open file with AVHRR information and read image
     * data and information
     */
    fprintf(stdout," Reading input AVHRR data...\n");
    fprintf(stdout," %s\n", fname);
    fm_init_fmio_img(&img);
    if (fm_readdata(infile, &img)) {
	fmerrmsg(where,"Could not open file...\n");
	exit(FM_IO_ERR);
    }

    if (strstr(img.sa,"NOAA-15") ||  strstr(img.sa,"NOAA-16")) {
	printf(" File discarded due to operational constraints on channel 3A\n");
	fm_clear_fmio_img(&img);
	exit(FM_VAROUTOFSCOPE_ERR);
	
    }
    printf(" Satellite: %s\n", img.sa);
    printf(" Time: %02d/%02d/%4d %02d:%02d\n", img.dd, img.mm, img.yy,
    img.ho, img.mi);
    sprintf(iinfo.satellite,"%s",img.sa);
    iinfo.hour = img.ho;
    iinfo.minute = img.mi;
    iinfo.day = img.dd;
    iinfo.month = img.mm;
    iinfo.year = img.yy;
    iinfo.zsize = img.z;
    iinfo.xsize = img.iw;
    iinfo.ysize = img.ih;
    iinfo.Ax = img.Ax;
    iinfo.Ay = img.Ay;
    iinfo.Bx = img.Bx;
    iinfo.By = img.By;
    size = img.iw*img.ih;

    fm_img2fmtime(img,&reftime);
    fm_img2fmucsref(img,&refucs);


    /*
     * Get NWP data...
     * This is probably not necessary in the future, but is kept until it
     * is clear whether T4 will be used to avoid cloud contamination or not.
     */

    nwpice_init(&nwp);

    #ifdef AVHRRICE_HAVE_NWP
    if (nwpice_read("xxx",reftime,refucs,&nwp)) {
    	fmerrmsg(where,"No NWP data available.");
    	fm_clear_fmio_img(&img);
    	nwpice_free(&nwp);
    	exit(FM_IO_ERR);
    }
    #endif

    /*
     * Get land/sea mask, accepted if within 0.5 km of the image.
     * This should probably be added agin in the future, but is left out
     * until further evaluation...
     */
    
    lm.d = NULL;
    if (lmask_located = fopen(lmaskf,"r")) {
      fprintf(stdout," Reading land/sea mask (GTOPO30 based):\n %s\n", lmaskf);
      status = read_hdf5_product(lmaskf, &lm, 0);
      if (status != 0) {
	fprintf(stderr,"%s\n"," Trouble processing:");
	fprintf(stderr,"%s\n",infile);
	fprintf(stderr,"%s%s\n", fmerrmsg,"Could not read land/sea mask");
	return(1);
      }
      fprintf(stdout," Checking for area consistency with land/sea mask...\n");
      if (((int) floorf(lm.h.Bx*10.)) != ((int) floorf(img.Bx*10.)) || 
	  ((int) floorf(lm.h.By*10.)) != ((int) floorf(img.By*10.)) ||
	  ((int) floorf(lm.h.Ax*10.)) != ((int) floorf(img.Ax*10.)) || 
	  ((int) floorf(lm.h.Ay*10.)) != ((int) floorf(img.Ay*10.)) ||
	  lm.h.iw != img.iw || lm.h.ih != img.ih) {
	fprintf(stderr,"%s\n"," Trouble processing:");
	fprintf(stderr,"%s\n",infile);
	fprintf(stderr,"%s%s",fmerrmsg,"Inconsistency between land/sea mask ");
	fprintf(stderr,"%s\n","and data input.");
	fprintf(stderr," Ax: %.2f %.2f\n", lm.h.Ax, img.Ax);
	fprintf(stderr," Ay: %.2f %.2f\n", lm.h.Ay, img.Ay);
	fprintf(stderr," Bx: %.2f %.2f\n", lm.h.Bx, img.Bx);
	fprintf(stderr," By: %f %f\n", lm.h.By, img.By);
	fprintf(stderr," iw: %d %d\n", lm.h.iw, img.iw);
	fprintf(stderr," ih: %d %d\n", lm.h.ih, img.ih);
	return(1);
      }
    }
    else {
      fprintf(stdout," Landmask unavailable, coefficients for ");
      fprintf(stdout,"sea/ice/cloud will be used throughout\n");
    }
    
    
 

    /*
     * Loading the statistical coeffs into statcoeffs struct           
     */
    fprintf(stdout," Reading file containing statistical coefficients:\n %s\n",
	    coffile);
      
    if (rdstatcoeffs(coffile,&coeffs) != 0) {
      sprintf(what," Trouble when reading coefficients, exiting");
      fmerrmsg(where,what);
      exit(FM_IO_ERR);
    }


    /*
     * Function "process_pixels4ice" is called to perform the objective
     * classification of the present satellite scene. Further description
     * of the function is given in the code.
     */
    init_osihdf(&ice);
    sprintf(ice.h.source, "%s", img.sa);
    sprintf(ice.h.product, "%s", where);
    ice.h.iw = img.iw;
    ice.h.ih = img.ih;
    ice.h.z = OSI_OLEVELS;
    ice.h.Ax = img.Ax;
    ice.h.Ay = img.Ay;
    ice.h.Bx = img.Bx;
    ice.h.By = img.By;
    ice.h.year = img.yy;
    ice.h.month = img.mm;
    ice.h.day = img.dd;
    ice.h.hour = img.ho;
    ice.h.minute = img.mi;
    status = malloc_osihdf(&ice,ice_ft,ice_desc);

    classed = (unsigned char *) malloc(size*sizeof(char));
    if (!classed) {
	sprintf(what,
	"Could not allocate memory for classed array while processing : %s\n",infile);
	fmerrmsg(where,what);
	exit(FM_MEMALL_ERR);
    }
    fmlogmsg(where,"Estimating ice probability");

    if (lm.d == NULL) {
      status = process_pixels4ice(img, NULL, NULL, nwp,
				ice.d, classed, 2, coeffs);
    }
    else {
      status = process_pixels4ice(img, NULL, (unsigned char *)(lm.d->data), 
				  nwp, ice.d, classed, 2, coeffs);
    }
    
    if (status) {
	sprintf(what,"Something failed while processing pixels of %s",infile);
	fmerrmsg(where,what);
    } else {
	fmlogmsg(where,"Finished estimating ice probability");
    }
    
    /*
     * Clean up memory used.
     *
     * Should add freeing of lmask here if needed in future...
     */
    fmlogmsg(where,"Cleaning memory");
    fm_clear_fmio_img(&img);

    /*
     * Write results to files, HDF5 file for internal use and TIFF 6.0 
     * (MITIFF) file for visual presentation on Internet/DIANA etc.
     *
     * MITIFF generation will be moved to a separate application in
     * time...
     */
    sprintf(clinfo.satellite,"%s",img.sa);
    clinfo.hour = img.ho;
    clinfo.minute = img.mi;
    clinfo.day = img.dd;
    clinfo.month = img.mm;
    clinfo.year = img.yy;
    clinfo.zsize = 1;
    clinfo.xsize = img.iw;
    clinfo.ysize = img.ih;
    clinfo.Ax = img.Ax;
    clinfo.Ay = img.Ay;
    clinfo.Bx = img.Bx;
    clinfo.By = img.By;

    opfn = (char *) malloc(FILELEN+5);
    if (!opfn) exit(3);
    sprintf(opfn,"%s/ice_%s_%4d%02d%02d%02d%02d.hdf5", 
	cfg.productpath,pname,
	img.yy, img.mm, img.dd, img.ho, img.mi);
    sprintf(what,"Creating output file: %s", opfn);
    fmlogmsg(where,what);
    status = store_hdf5_product(opfn,ice);
    free(opfn);
    if (status != 0) {
	sprintf(what,"Trouble processing: %s",infile);
	fmerrmsg(where,what);
    }
    /*
    opfn = (char *) malloc(FILELEN+5);
    if (!opfn) exit(3);
    sprintf(opfn,"%s/ice_%s_%4d%02d%02d%02d%02d.mitiff", 
	cfg.productpath,pname,
	img.yy, img.mm, img.dd, img.ho, img.mi);
    sprintf(what,"Creating output file: %s", opfn);
    fmlogmsg(where,what);
    store_mitiff_result(opfn,classed,clinfo);
    free(opfn);
    */

    fprintf(stdout," ================================================\n");

    exit(FM_OK);
}

/*
 * NAME:
 * usage
 *
 * PURPOSE:
 * To give info to the user on syntax etc.
 */
void usage() {
    fprintf(stdout,"\n");
    fprintf(stdout," SYNTAX:\n");
    fprintf(stdout,
	    " ice_avhrr -c <cfgfile> -i <infile>\n\n");
    fprintf(stdout,
	    " <cfgfile>: Configuration file containing data paths etc.\n");
    fprintf(stdout,
	    " <infile>: Input METSAT file, path is taken from cfgfile.\n");
    fprintf(stdout,"\n");
    fprintf(stdout," The configuration file contains all necessary data\n");
    fprintf(stdout," paths for production of ice tiles. Output names are\n");
    fprintf(stdout," automatically created.\n");
    fprintf(stdout,"\n");
    exit(FM_OK);
}

/*
 * NAME:
 * decode_cfg
 *
 * PURPOSE:
 * To decode configuration file.
 *
 * NOTES:
 * Standard return values 0, 2 and 3, OK, I/O and Memory trouble.
 */

int decode_cfg(char cfgfile[],cfgstruct *cfg) {
    FILE *fp;
    char *errmsg=" ERROR(decode_cfg): ";
    char *dummy,*pt;
    char *token=" ";

    dummy = (char *) malloc(FILELEN*sizeof(char));
    if (!dummy) {
	fprintf(stderr,"%s%s\n",errmsg,"Could not allocate memory");
	return(FM_MEMALL_ERR);
    }

    fp = fopen(cfgfile,"r");
    if (!fp) {
	fprintf(stderr,"%s%s\n",errmsg,"Could not open config file.");
	return(FM_IO_ERR);
    }

    while (fgets(dummy,FILELEN,fp) != NULL) {
	if (strncmp(dummy,"#",1) == 0) continue;
	if (strlen(dummy) > (FILELEN-50)) {
	    fprintf(stderr,"%s%s\n",errmsg,
		    "Input string larger than FILELEN");
	    free(dummy);
	    return(FM_IO_ERR);
	}

	pt = strtok(dummy,token);

	if (!pt) {
	    fprintf(stderr,"%s%s\n",errmsg,"strtok trouble.");
	    free(dummy);
	    return(FM_IO_ERR);
	}
	if (strncmp(pt,"IMGPATH",7) == 0) {
	    pt = strtok(NULL,token);
	    if (!pt) {
		fprintf(stderr,"%s%s\n",errmsg,"strtok trouble.");
		free(dummy);
		return(FM_IO_ERR);
	    }
	    fmremovenewline(pt);
	    sprintf(cfg->imgpath,"%s",pt);
	} else if (strncmp(pt,"CMPATH",6) == 0) {
	    pt = strtok(NULL,token);
	    if (!pt) {
		fprintf(stderr,"%s%s\n",errmsg,"strtok trouble.");
		free(dummy);
		return(FM_IO_ERR);
	    }
	    fmremovenewline(pt);
	    sprintf(cfg->cmpath,"%s",pt);
	} else if (strncmp(pt,"LMPATH",6) == 0) {
	    pt = strtok(NULL,token);
	    if (!pt) {
		fprintf(stderr,"%s%s\n",errmsg,"strtok trouble.");
		free(dummy);
		return(FM_IO_ERR);
	    }
	    fmremovenewline(pt);
	    sprintf(cfg->lmpath,"%s",pt);
	} else if (strncmp(pt,"PRODUCTPATH",11) == 0) {
	    pt = strtok(NULL,token);
	    if (!pt) {
		fprintf(stderr,"%s%s\n",errmsg,"strtok trouble.");
		free(dummy);
		return(FM_IO_ERR);
	    }
	    fmremovenewline(pt);
	    sprintf(cfg->productpath,"%s",pt);
	}else if (strncmp(pt,"PROBTABNAME",11) == 0) {
	    pt = strtok(NULL,token);
	    if (!pt) {
		fprintf(stderr,"%s%s\n",errmsg,"strtok trouble.");
		free(dummy);
		return(FM_IO_ERR);
	    }
	    fmremovenewline(pt);
	    sprintf(cfg->probtabname,"%s",pt);
	}
    }

    fclose(fp);

    free(dummy);

    return(FM_OK);
}


/*
 * NAME:
 * rdstatcoeffs
 *
 * PURPOSE:
 * To read the statistical coefficients needed in probest from file
 * with name&path given in config-file.
 *
 */

int rdstatcoeffs (char *coeffsfile, statcoeffstr *cof){
  char *where="rdstatcoeffs";
  char what[OSI_MSGLENGTH];
  FILE *fpi;
  char *line = NULL;
  ssize_t read;
  size_t len = 0;
  dummystr dummies = {" "," ",' ',0.,0.,0.};
  int j;
  
  fpi = fopen(coeffsfile,"r");
  if (!fpi) { 
    sprintf(what,"Unable to open file %s",coeffsfile);
    fmerrmsg(where,what);
    return(FM_IO_ERR);
  }

  while ( (read = getline(&line,&len,fpi)) != -1 ) {
    j = 0;
    while ((line[j] == ' ') && j < OSI_MSGLENGTH) j++;
    if (j >= OSI_MSGLENGTH) {
      sprintf(what,"Line length exceeds maximum length");
      fmerrmsg(where,what);
      return(FM_IO_ERR);
    }
    if ((line[j] == '#') || line[j] == '\n') continue;
    if (sscanf(line,"%s%s %c%lf%lf%lf",dummies.surf,dummies.feat,
	       &dummies.key,&dummies.par1,&dummies.par2,&dummies.par3) != 6) {
      sprintf(what,"Wrong format on line '%s %s..'\n",dummies.surf,
	      dummies.feat);
      fmerrmsg(where,what);
      return(FM_IO_ERR);
    }
    else {
      locstatcoeffs(dummies,cof);
      sprintf(dummies.surf," ");
      sprintf(dummies.feat," ");
      dummies.key = ' ';
      dummies.par1 = dummies.par2 = dummies.par3 = 0.;
    }
  }

  if (line) free(line);
  fclose(fpi);
  
  return(FM_OK);
}


/*search through statcoeffstr to find right location for parameters*/
int locstatcoeffs (dummystr dummies, statcoeffstr *cof){
  char *where="locstatcoeffs";
  char what[OSI_MSGLENGTH];
  int featflg,surfflg;

  featflg = surfflg = 0;
 
  if (!strcmp(dummies.surf,"ice")) { /*strcmp returns 0 if match!*/
    if (!strcmp(dummies.feat,"a1")) putcoeffs(&cof->ice.a1,dummies);
    else if (!strcmp(dummies.feat,"r21")) putcoeffs(&cof->ice.r21,dummies);
    else if (!strcmp(dummies.feat,"d34")) putcoeffs(&cof->ice.d34,dummies);
    else if (!strcmp(dummies.feat,"r3a1")) putcoeffs(&cof->ice.r3a1,dummies);
    else if (!strcmp(dummies.feat,"r3b1")) putcoeffs(&cof->ice.r3b1,dummies);
    else if (!strcmp(dummies.feat,"dt"))  putcoeffs(&cof->ice.dt,dummies);
    else featflg++;
  } 
  else if (!strcmp(dummies.surf,"cloud")) {
    if (!strcmp(dummies.feat,"a1")) putcoeffs(&cof->cloud.a1,dummies);
    else if (!strcmp(dummies.feat,"r21")) putcoeffs(&cof->cloud.r21,dummies);
    else if (!strcmp(dummies.feat,"d34")) putcoeffs(&cof->cloud.d34,dummies);
    else if (!strcmp(dummies.feat,"r3a1")) putcoeffs(&cof->cloud.r3a1,dummies);
    else if (!strcmp(dummies.feat,"r3b1")) putcoeffs(&cof->cloud.r3b1,dummies);
    else if (!strcmp(dummies.feat,"dt"))  putcoeffs(&cof->cloud.dt,dummies);
    else featflg++;
  }
  else if (!strcmp(dummies.surf,"water")) {
    if (!strcmp(dummies.feat,"a1")) putcoeffs(&cof->water.a1,dummies);
    else if (!strcmp(dummies.feat,"r21")) putcoeffs(&cof->water.r21,dummies);
    else if (!strcmp(dummies.feat,"d34")) putcoeffs(&cof->water.d34,dummies);
    else if (!strcmp(dummies.feat,"r3a1")) putcoeffs(&cof->water.r3a1,dummies);
    else if (!strcmp(dummies.feat,"r3b1")) putcoeffs(&cof->water.r3b1,dummies);
    else if (!strcmp(dummies.feat,"dt"))  putcoeffs(&cof->water.dt,dummies);
    else featflg++;
  }
  else if (!strcmp(dummies.surf,"land")) {
    if (!strcmp(dummies.feat,"a1")) putcoeffs(&cof->land.a1,dummies);
    else if (!strcmp(dummies.feat,"r21")) putcoeffs(&cof->land.r21,dummies);
    else if (!strcmp(dummies.feat,"d34")) putcoeffs(&cof->land.d34,dummies);
    else if (!strcmp(dummies.feat,"r3a1")) putcoeffs(&cof->land.r3a1,dummies);
    else if (!strcmp(dummies.feat,"r3b1")) putcoeffs(&cof->land.r3b1,dummies);
    else if (!strcmp(dummies.feat,"dt"))  putcoeffs(&cof->land.dt,dummies);
    else featflg++;
  }
  else surfflg++;

  if (featflg || surfflg) {
    if (featflg) sprintf(what,"Feature '%s' not recognised for surface %s\n",
	      dummies.feat,dummies.surf);
    if (surfflg) sprintf(what,"Surface '%s' not recognised\n",dummies.surf);
    fmerrmsg(where,what);
    exit(0);
  }

  return(FM_OK);
}

/*place stat. parameters at right location in statcoeffstr*/
int putcoeffs(featstr *feat, dummystr dummies) {
  feat->key  = dummies.key;
  feat->par1 = dummies.par1;
  feat->par2 = dummies.par2;
  feat->par3 = dummies.par3;
  feat->count++; /*will equal the number of times a specific set of
		    parameters are read. Should ideally not differ from one..*/
  return(FM_OK);
}


