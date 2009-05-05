/*
 * NAME: 
 * fmaccusnow
 * 
 * PURPOSE: 
 * Average (AVHRR) fmsnowcover product files to compound a cloudfree
 * snowmap. The resulting product is tagged with the end time of the
 * integration period.
 * 
 * SYNTAX: accusnow -s <dir_fmsnow> -d <date_end> 
 *         -p <period> -a <pref_outf> -o <path_outf>
 *         (-t <satellite> -l <satlist> -m <arealist> -z)
 *
 *    <dir_fmsnow>  : Directory with hdf5 files with fmsnow data.
 *    <date_end>     : End date of merging period.
 *    <period>       : Integration period in hours.
 *    <pref_outf>    : Prefix for output filename.
 *    <path_outf>    : Path for output file.
 *    <satellite>    : Name of one satellite if processing only one (optional).
 *    <satlist>      : File with satellites to use (optional).
 *    <arealist>     : File with tile areas to use (optional).
 *    -z             : Use threshold on satellite zenith angle (value from header file).
 *
 * NOTE:
 * NA
 * 
 * AUTHOR: 
 * Steinar Eastwood, DNMI, 21.08.2000
 *
 * MODIFIED: 
 * Mari Anne Killie, METNO/FOU, 08.01.2009: Original software created by
 * Steinar Eastwoord, modified for use within the fmsnowcover package.
 * Øystein Godøy, METNO/FOU, 23.04.2009: More cleaning of software.
 *
 * CVS_ID:
 * $Id: fmaccusnow.c,v 1.6 2009-05-05 12:34:23 steingod Exp $
 */ 

 
#include <fmaccusnow.h>
#include <unistd.h>

int main(int argc, char *argv[]) {    
    char *where="fmaccusnow";
    extern char *optarg;
    char *dir_avhrrice, *date_start, *date_prod, *date_end;
    int sflg, dflg, pflg, aflg, oflg, tflg, lflg, mflg, zflg, cflg;
    int period, i, j, f, t, tile, nrInput, ret, ind, numf;
    int numsat, numarea;
    fmsec1970 stime, ftime, etime, prodtime;
    float cloudlim;
    char **infile_avhrrice, *outfHDF, *outfMITIFF_class, *outfMITIFF_psnow;
    char **infile_sorted,**infile_currenttile;
    char *pref_outf, *path_outf, *checkfile, *sret, datestr[13], *procsat; 
    char datestr_ymdhms[15];
    char *satlistfile, *arealistfile, **satlist, **arealist;
    fmtime timedate;
    fmucsref safucs;
    struct dirent *dirl_avhrrice;
    DIR *dirp_avhrrice;
    char *defpref = "accusnow";
    char *pref_ps = "sp";
    char *pref_cl = "cl";
    int satfound;
    char *satstring; /*To be used in output filenames*/
    osihdf snowprod, checkfileheader; /*To be renamed*/
    osihdf *inputhdf;
    char *prod_desc[FMACCUSNOWPROD_LEVELS] = {"class","P(snow)","P(clear)"};
    osi_dtype prod_ft[FMACCUSNOWPROD_LEVELS] = {CLASS_DT,PROB_DT,PROB_DT};
    unsigned char *class, *snowclass;
    float *probsnow, *probclear;
    int sorted_index, index_offset;
    int *num_files_area, *num_files_area_counter;
    char *default_arealist[TOTAREAS] = {"ns","nr","at","gr","gn","gf","gm","gs"};
    char *class_desc[CLASSLIMITS] = {
	"1: Ice/snow",
	"2: Clear",
	"3: Unclass",
	"4: Clouded",
	"5: Undef"
    };
    char checkstr[7];
    char *snow_desc[PROBLIMITS] = {"[0.00 -  0.05>", "[0.05 -  0.10>", 
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

    if (!(argc >= 9 && argc <= 18)) usage();

    fprintf(stdout,"\n");
    fprintf(stdout,"\t=================================================\n");
    fprintf(stdout,"\t|                  FMACCUSNOW                     |\n");
    fprintf(stdout,"\t=================================================\n");
    fprintf(stdout,"\n");

    /* Interprete commandline arguments */
    sflg=dflg=pflg=aflg=oflg=tflg=lflg=mflg=zflg=cflg=0;
    while ((ret = getopt(argc, argv, "s:d:p:a:o:t:l:m:c:z")) != EOF) {
	switch (ret) {
	    case 's':
		dir_avhrrice = (char *) malloc(strlen(optarg)+1);
		if (! dir_avhrrice) {
		    fmerrmsg(where,"Could not allocate dir_avhrrice");
		    exit(FM_MEMALL_ERR);
		}
		if (!strcpy(dir_avhrrice, optarg)) exit(1);
		sflg++;
		break;
	    case 'd':
		date_end = (char *) malloc(DATESTRINGLENGTH*sizeof(char));
		if (! date_end) {
		    fmerrmsg(where,"Could not allocate date_end");
		    exit(FM_MEMALL_ERR);
		}
		if (!strcpy(date_end, optarg)) exit(FM_IO_ERR);
		dflg++;
		break;
	    case 'p':
		period = atoi(optarg);
		pflg++;
		break;
	    case 'a':
		pref_outf = (char *) malloc(strlen(optarg)+1);
		if (! pref_outf) {
		    fmerrmsg(where,"Could not allocate pref_outf");
		    exit(FM_MEMALL_ERR);
		}
		if (!strcpy(pref_outf, optarg)) exit(FM_IO_ERR); 
		aflg++; 
		break;
	    case 'o':
		path_outf = (char *) malloc(strlen(optarg)+1);
		if (! path_outf) {
		    fmerrmsg(where,"Could not allocate path_outf");
		    exit(FM_MEMALL_ERR);
		}
		if (!strcpy(path_outf, optarg)) exit(FM_IO_ERR);
		oflg++;
		break;
	    case 't':
		procsat = (char *) malloc(strlen(optarg)+1);
		if (! procsat) {
		    fmerrmsg(where,"Could not allocate procsat");
		    exit(FM_MEMALL_ERR);
		}
		if (!strcpy(procsat, optarg)) exit(FM_IO_ERR);
		tflg++;
		break;
	    case 'l':
		satlistfile = (char *) malloc(strlen(optarg)+1);
		if (! satlistfile) {
		    fmerrmsg(where,"Could not allocate satlistfile");
		    exit(FM_MEMALL_ERR);
		}
		if (!strcpy(satlistfile, optarg)) exit(FM_IO_ERR);
		lflg++;
		break;
	    case 'm':
		arealistfile = (char *) malloc(strlen(optarg)+1);
		if (! arealistfile) {
		    fmerrmsg(where,"Could not allocate arealistfile");
		    exit(FM_MEMALL_ERR);
		}
		if (!strcpy(arealistfile, optarg)) exit(FM_IO_ERR);
		mflg++;
		break;
	    case 'c':
		cloudlim = atof(optarg);
		cflg++;
		break;
	    case 'z':
		zflg++;
		break;
	    default:
		usage();
	}
    }

    if (!sflg || !dflg || !pflg || !oflg) usage();
    if (lflg && tflg) {
	fprintf(stdout,
		"\n ERROR: do not give arguments l and t simultaneously\n\n");
	exit(FM_OK);
    }

    if (strlen(date_end) == 8 ) strcat(date_end,"00");
    if (strlen(date_end) != 10) {
	fmerrmsg(where,"Incorrect length of date_end, should be yyyymmddhh.");
	exit(FM_IO_ERR);
    }

    etime = ymdh2fmsec1970(date_end,0);
    stime = etime-period*3600;
    date_start = (char *) malloc(DATESTRINGLENGTH*sizeof(char));
    if (! date_start) {
	fmerrmsg(where,"Could not allocate date_start.");
	exit(FM_MEMALL_ERR);
    }
    fmsec19702isodatetime(stime, date_start);
    fmsec19702isodatetime(etime, date_end);
    prodtime = etime;
    date_prod= (char *) malloc(DATESTRINGLENGTH*sizeof(char));
    if (! date_prod) {
	fmerrmsg(where,"Could not allocate date_prod.");
	exit(FM_MEMALL_ERR);
    }
    fmsec19702isodatetime(prodtime, date_prod);
    if (tofmtime(prodtime,&timedate)) {
	fmerrmsg(where,"tofmtime failed on prodtime.");
	exit(FM_IO_ERR);
    }

    fprintf(stdout,"\tAVHRR ice files dir:   %s\n", dir_avhrrice);
    fprintf(stdout,"\tPeriod:                %d hours\n", period);
    fprintf(stdout,"\tDate start:            %s \n", date_start);
    fprintf(stdout,"\tDate end:              %s \n", date_end);
    fprintf(stdout,"\tDate product:          %s \n", date_prod);
    fprintf(stdout,"\tPeriod start:          %d \n", (int)stime);
    fprintf(stdout,"\tPeriod end:            %d \n", (int)etime);
    fprintf(stdout,"\tPath outfile:          %s \n", path_outf);
    if (aflg) {
	fprintf(stdout,"\tPrefix outfile:        %s \n", pref_outf);
    } else {
	pref_outf = defpref;
    }

    if (tflg) {
	fprintf(stdout,"\tProcess only satellite: %s \n", procsat);
    }
    if (lflg) {
	fprintf(stdout,"\tUsing list of satellites from file: %s \n", satlistfile);
    }
    if (!cflg) {
	cloudlim = DEFAULTCLOUD;
    }
    fprintf(stdout,"\tUsing cloud probability limit: %.1f \n", cloudlim);
    if (mflg) {
	fprintf(stdout,"\tUsing area tiles from file: %s \n", arealistfile);
    }

    satlist = (char **) malloc(MAXSAT*sizeof(char *)); /*flytte til etter if lfgl?*/
    if (! satlist) {
	fmerrmsg(where,"Could not allocate satlist");
	exit(FM_MEMALL_ERR);
    }

    /* 
     * Read satellite name list 
     */
    numsat = 0;
    if (lflg) {
	numsat = read_sat_area_list(satlistfile,satlist);
	if (numsat == 0) {
	    fprintf(stderr,"\n\tFound NO satellites on file %s\n\texiting..\n\n",
		    satlistfile);
	} else if (numsat < 0) {
	    fmerrmsg(where,"Could not read %s, processing all passages.",satlistfile);
	    lflg = 0;
	} else {
	    fprintf(stdout,"\n\tProcessing satellites (%d): \n\t  ",numsat);
	    satstring = (char *) malloc(LISTLEN*sizeof(char));
	    if (! satstring) {
		fmerrmsg(where,"Could not allocate satstring");
		exit(FM_MEMALL_ERR);
	    }
	    for (i=0;i<numsat;i++) { 
		fprintf(stdout,"%s ",satlist[i]);
		if (i != 0) strcat(satstring,"_");
		strcat(satstring,satlist[i]);
	    }
	    fprintf(stdout,"\n");
	}
    } else if (tflg) {
	satstring = (char *) malloc(strlen(procsat)*sizeof(char));
	if (! satstring) {
	    fmerrmsg(where,"Could not allocate satstring");
	    exit(FM_MEMALL_ERR);
	}
	sprintf(satstring,"%s",procsat);
    } else { /* No sat.name list given, need componet for filename */
	satstring = (char *) malloc(LISTLEN*sizeof(char));
	if (! satstring) {
	    fmerrmsg(where,"Could not allocate satstring");
	    exit(FM_MEMALL_ERR);
	}
	sprintf(satstring,"allsats");
	fprintf(stdout,"\tProcessing all available satnames.\n");
    }

    numarea = 0;
    arealist = (char **) malloc(MAXAREA*sizeof(char *));
    if (! arealist) {
	fmerrmsg(where,"Could not allocate arealist");
	exit(FM_MEMALL_ERR);
    }
    if (mflg) {/* Read tile list */
	numarea = read_sat_area_list(arealistfile,arealist);
	if (numarea == 0) {
	    fmerrmsg(where,"Found no area tiles on %s",arealistfile);
	    exit(FM_IO_ERR);
	}
	else if (numarea < 0) {
	    fmerrmsg(where,"Could not read %s, processing all area tiles.", arealistfile);
	    mflg = 0;
	}
	else {
	    fprintf(stdout,"\tProcessing area tiles (%d):\n",numarea);
	    for (i=0;i<numarea;i++) { 
		fprintf(stdout,"\t%s\n",arealist[i]);
	    }
	}
    }
    if (mflg == 0) { /*hard code a list of all known areas*/
	numarea = TOTAREAS;
	arealist = default_arealist;
	fprintf(stdout,
		"\tProcessing files for tiles on default tilelist (%d):\n\t ",numarea);
	for (i=0;i<numarea;i++) {
	    printf("%s ",arealist[i]);
	}
	fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");

    /* 
     * Allocate and initialise variables to handle sorting of files 
     * based on tile
     */
    num_files_area = (int *)malloc(numarea*sizeof(int));
    if (! num_files_area) {
	fmerrmsg(where,"Could not allocate num_files_area");
	exit(FM_MEMALL_ERR);
    }
    num_files_area_counter = (int *)malloc(numarea*sizeof(int));
    if (! num_files_area_counter) {
	fmerrmsg(where,"Could not allocate num_files_area_counter");
	exit(FM_MEMALL_ERR);
    }
    for (i=0;i<numarea;i++){
	num_files_area[i] = num_files_area_counter[i] = 0;
    }

    /* 
     * Reading directory, find all avhrrice files to process. 
     */
    dirp_avhrrice = opendir(dir_avhrrice);
    if (!dirp_avhrrice) {
	fmerrmsg(where,"Could not open %s",dir_avhrrice);
	return(FM_IO_ERR);
    }

    numf = 0;
    while ((dirl_avhrrice = readdir(dirp_avhrrice)) != NULL) {
	if (strncmp(dirl_avhrrice->d_name,BASEFNAME,strlen(BASEFNAME)) == 0 &&
		strstr(dirl_avhrrice->d_name,".hdf") != NULL && strlen(dirl_avhrrice->d_name) >= MINLENFNAME) {
	    sret = strncpy(datestr,&dirl_avhrrice->d_name[10],12);
	    datestr[12] = '\0';
	    sprintf(datestr_ymdhms,"%s00",datestr);
	    ftime = ymdhms2fmsec1970(datestr_ymdhms,0);
	    if (ftime >= stime && ftime <= etime) {
		numf ++;
	    }
	}
    }
    rewinddir(dirp_avhrrice);

    infile_avhrrice = (char **) malloc(numf*sizeof(char *));
    if (! infile_avhrrice) {
	fmerrmsg(where,"Could not allocate infile_avhrrice");
	exit(FM_MEMALL_ERR);
    }
    infile_sorted   = (char **) malloc(numf*sizeof(char *));
    if (! infile_sorted) {
	fmerrmsg(where,"Could not allocate infile_sorted");
	exit(FM_MEMALL_ERR);
    }
    i = 0;

    /*
     * Loop through files in input directory.
     */
    while ((dirl_avhrrice = readdir(dirp_avhrrice)) != NULL) {

	/*
	 * Check input filename
	 */
	if (strncmp(dirl_avhrrice->d_name,BASEFNAME,strlen(BASEFNAME)) 
		== 0 &&	strstr(dirl_avhrrice->d_name,".hdf") != NULL && 
		strlen(dirl_avhrrice->d_name) >= MINLENFNAME) {

	    /*
	     * Check time of file
	     */
	    sret = strncpy(datestr,&dirl_avhrrice->d_name[10],12);
	    datestr[12] = '\0';
	    sprintf(datestr_ymdhms,"%s00",datestr);
	    ftime = ymdhms2fmsec1970(datestr_ymdhms,0);
	    if (ftime < stime || ftime > etime) {
		continue;
	    }

	    /*
	     * Check that file can be opened (this removes files of size
	     * zero!
	     */
	    checkfile = (char *) malloc(256*sizeof(char));
	    if (! checkfile) {
		fmerrmsg(where,"Could not allocate checkfile");
		exit(FM_MEMALL_ERR);
	    }
	    sprintf(checkfile,"%s/%s",dir_avhrrice,dirl_avhrrice->d_name);
	    init_osihdf(&checkfileheader);
	    ret = read_hdf5_product(checkfile,&checkfileheader,1);
	    free(checkfile);
	    if (ret != 0) {
		fmerrmsg(where,"Could not open %s, skipping file",dirl_avhrrice->d_name);
		free_osihdf(&checkfileheader);
		continue;
	    }

	    /*
	     * Check satellite name. satname is not part of input
	     * filename, must  check the header to see if file should be
	     * kept.
	     */
	    if (lflg || tflg) { 
		satfound = 0;
		if (tflg) {
		    if (strstr(checkfileheader.h.source,procsat)) {
			satfound++;
		    } 
		}
		else {
		    for (j=0;j<numsat;j++) { 
			if (strstr(checkfileheader.h.source,satlist[j])) {
			    satfound++;
			}
		    }
		} 
		free_osihdf(&checkfileheader);
		if (!satfound) {
		    continue;
		}
	    }

	    /*
	     * Check tile, this returns the element number.
	     */
	    ind = find_sat_area_index(arealist,numarea,dirl_avhrrice->d_name);
	    if (ind < 0) {
		fprintf(stdout,"\tFile not in area tile list, skipping %s\n",
			dirl_avhrrice->d_name);
		continue;
	    }

	    /*
	     * If all tests passed, add to list of ok files :-)
	     */
	    infile_avhrrice[i] = (char *) malloc(256*sizeof(char));
	    if (! infile_avhrrice[i]) {
		fmerrmsg(where,"Could not allocate infile_avhrrice[%d]", i);
		exit(FM_MEMALL_ERR);
	    }
	    infile_sorted[i]   = (char *) malloc(256*sizeof(char));
	    if (! infile_sorted[i]) {
		fmerrmsg(where,"Could not allocate infile_sorted[%d]", i);
		exit(FM_MEMALL_ERR);
	    }
	    sprintf(infile_avhrrice[i],"%s/%s",dir_avhrrice,dirl_avhrrice->d_name);
	    i ++;
	    num_files_area[ind]++;
	}
    }
    free(dir_avhrrice);
    nrInput = i;

    if (nrInput == 0) {
	fmerrmsg(where,"No files to be processed.");
	exit(FM_OK);
    }

    fmlogmsg(where,"Found a number of files to integrate over. %d files fulfilled the time constraints, %d files fulfilled the satid constraints.", numf, nrInput);
    for (i=0;i<nrInput;i++) {
	fprintf(stdout,"\t%2d %s\n",i,infile_avhrrice[i]);
    }

    fmlogmsg(where,"Sorting the file to process.");
    for (t=0;t<numarea;t++) {
	fprintf(stdout,
		"\t\tTile: %s  Number of files: %d\n",
		arealist[t],num_files_area[t]);
	num_files_area_counter[t] = 0;
    }
    for (i=0;i<nrInput;i++) {
	for (t=0;t<numarea;t++) {
	    sprintf(checkstr,"fmsnow_%2s",arealist[t]);
	    if (strstr(infile_avhrrice[i],checkstr)) {
		sorted_index = 0;
		for (j=0;j<t;j++) {
		    sorted_index += num_files_area[j];
		}
		sorted_index+=num_files_area_counter[t];
		sprintf(infile_sorted[sorted_index],infile_avhrrice[i]) ;
		num_files_area_counter[t]++;
	    }
	}
    }

    /*
     * Check the time to be associated with the integrated product (i.e.
     * the time of the integration.
     */
    if (timedate.fm_hour > 23 || timedate.fm_min > 59 || 
	timedate.fm_mday > 31 || timedate.fm_mon > 12 || 
	timedate.fm_year > 3000 || timedate.fm_year < 1978)  {
	fmerrmsg(where,"Incorrect input date and time.");
	exit(FM_IO_ERR);
    }

    /*Dette kan hentes fra første fil. Tar utgangspunkt i SAF-tilene. Om
      data fra andre sensorer skal innlemmes seinere får de heller
      tilpasse seg AVHRR.*/
    /*kan dette som følger herfra settes i loop for flere tiler?*/

    for (tile=0;tile<numarea;tile++){

	if (num_files_area[tile]==0) {
	    fprintf(stdout,"\t No input files for tile %s, continuing on list\n",
		    arealist[tile]);
	    continue;
	}

	index_offset=0;
	for (t=0;t<tile;t++) {
	    index_offset += num_files_area[t];
	}

	/*
	 * Creating list of files for this tile
	 */
	infile_currenttile = (char **) malloc(num_files_area[tile]*sizeof(char *));
	if (! infile_currenttile) {
	    fmerrmsg(where,"Could not allocate infile_currenttile");
	    exit(FM_MEMALL_ERR);
	}
	for (f=0;f<num_files_area[tile];f++) {
	    infile_currenttile[f] = (char *) malloc(256*sizeof(char));
	    if (! infile_currenttile[f]) {
		fmerrmsg(where,"Could not allocate infile_currenttile[%d]", f);
		exit(FM_MEMALL_ERR);
	    }
	    sprintf(infile_currenttile[f],infile_sorted[f+index_offset]);
	}
	/*for (f=0;f<num_files_area[tile];f++){
	  fprintf(stdout,"%s\n",infile_currenttile[f]);
	  }*/

	inputhdf = (osihdf *) malloc(num_files_area[tile]*sizeof(osihdf));
	if (! inputhdf) {
	    fmerrmsg(where,"Could not allocate inputhdf");
	    exit(FM_MEMALL_ERR);
	}

	init_osihdf(&snowprod); /*Dette init. kun headeren, ikke datafeltet!*/

	/*reading headers of all files on list to compare + collect info*/
	for (f=0;f<num_files_area[tile];f++) {
	    init_osihdf(&inputhdf[f]);
	    ret=read_hdf5_product(infile_currenttile[f],&inputhdf[f], 1);
	    if (ret != 0) {
		fmerrmsg(where,"Could not read header of %s",infile_currenttile[f]);
	    }
	    if (inputhdf[f].h.iw != inputhdf[0].h.iw ||
		    inputhdf[f].h.ih != inputhdf[0].h.ih ||
		    inputhdf[f].h.Ax != inputhdf[0].h.Ax ||
		    inputhdf[f].h.Ay != inputhdf[0].h.Ay ||
		    inputhdf[f].h.Bx != inputhdf[0].h.Bx ||
		    inputhdf[f].h.By != inputhdf[0].h.By ) {
		fmerrmsg(where,"Input files for tile %s, are from different tiles.",arealist[tile]);
		exit(FM_IO_ERR);
	    }
	}

	snowprod.h.iw = inputhdf[0].h.iw;
	snowprod.h.ih = inputhdf[0].h.ih;
	snowprod.h.z  = inputhdf[0].h.z;
	snowprod.h.Ax = inputhdf[0].h.Ax;
	snowprod.h.Ay = inputhdf[0].h.Ay;
	snowprod.h.Bx = inputhdf[0].h.Bx;
	snowprod.h.By = inputhdf[0].h.By;
	snowprod.h.year = timedate.fm_year;
	snowprod.h.month = timedate.fm_mon;
	snowprod.h.day = timedate.fm_mday;
	snowprod.h.hour = timedate.fm_hour;
	snowprod.h.minute = timedate.fm_min;
	sprintf(snowprod.h.area, "%s", inputhdf[0].h.area);
	sprintf(snowprod.h.source, "%s", inputhdf[0].h.source);
	sprintf(snowprod.h.product, "%s", inputhdf[0].h.product);
	sprintf(snowprod.h.projstr, "%s", inputhdf[0].h.projstr);

	safucs.Ax = snowprod.h.Ax;
	safucs.Ay = snowprod.h.Ay;
	safucs.Bx = snowprod.h.Bx;
	safucs.By = snowprod.h.By;
	safucs.iw = snowprod.h.iw;
	safucs.ih = snowprod.h.ih;


	class = (unsigned char *) malloc(safucs.iw*safucs.ih*sizeof(char));
	if (! class) {
	    fmerrmsg(where,"Could not allocate class");
	    exit(FM_MEMALL_ERR);
	}
	snowclass = (unsigned char *) malloc(safucs.iw*safucs.ih*sizeof(char));
	if (! snowclass) {
	    fmerrmsg(where,"Could not allocate snowclass");
	    exit(FM_MEMALL_ERR);
	}
	probsnow = (float *) malloc(safucs.iw*safucs.ih*sizeof(float));
	if (! probsnow) {
	    fmerrmsg(where,"Could not allocate probsnow");
	    exit(FM_MEMALL_ERR);
	}
	probclear = (float *) malloc(safucs.iw*safucs.ih*sizeof(float));
	if (! probclear) {
	    fmerrmsg(where,"Could not allocate probclear");
	    exit(FM_MEMALL_ERR);
	}

	if (!class || !snowclass || !probsnow || !probclear) {
	    fmerrmsg(where,"Could not allocate memory for class, snowclass, probsnow or probclear");
	    exit(FM_MEMALL_ERR);
	}

	for (i=0;i<safucs.iw*safucs.ih;i++) {
	    class[i]    = C_UNDEF; 
	    snowclass[i] = 0;
	    probsnow[i]  = PROB_MISVAL;
	    probclear[i]= PROB_MISVAL;
	}

	/*Her skjer midlinga!*/

	if (num_files_area[tile] > 0) {
	    fprintf(stdout,"\n\tNow averaging tile %s (%d files)..\n",
		    arealist[tile],num_files_area[tile]);
	    ret = average_merge_files(infile_currenttile, num_files_area[tile],
		    safucs, class, snowclass, probsnow, probclear, cloudlim); 
	    if (ret != 0) {
		fmerrmsg(where,"Could not finish average_merge_files");
		exit(FM_OTHER_ERR);
	    }
	}

	for (f=0;f<num_files_area[tile];f++) {
	    free(infile_currenttile[f]);
	    /*free_osihdf(&inputhdf[f]); Unødvendig siden ikke_datafeltet_ er allokert uansett. Denne freer kun datafeltet.*/
	}
	free(infile_currenttile);


	/* Allocate hdf5 product */
	ret = malloc_osihdf(&snowprod,prod_ft,prod_desc);
	if (ret != 0) {
	    fmerrmsg(where,"Could not run malloc_osihdf");
	    return(FM_MEMALL_ERR); 
	}

	for (i=0;i<safucs.iw*safucs.ih;i++){
	    ((int*)snowprod.d[0].data)[i] = class[i];
	    ((float*)snowprod.d[1].data)[i] = probsnow[i];
	    ((float*)snowprod.d[2].data)[i] = probclear[i];
	}

	fprintf(stdout,"\tCreating output files for tile %s:\n",arealist[tile]);

	/* create hdf product file */
	outfHDF = (char *) malloc(FILELEN+5);
	if (!outfHDF) exit(FM_MEMALL_ERR);
	sprintf(outfHDF,"%s/%s_%s_%04d%02d%02d-%dhours_%s.hdf",path_outf,
		pref_outf,arealist[tile],snowprod.h.year,snowprod.h.month,
		snowprod.h.day ,period,satstring);
	ret = store_hdf5_product(outfHDF, snowprod);
	if (ret != 0)  {
	    fmerrmsg(where,"Could not create HDF file %s", outfHDF);
	    exit(FM_IO_ERR);
	}    


	/* create mitiff for classified image */
	outfMITIFF_class = (char *) malloc(FILELEN+5);
	if (!outfMITIFF_class) exit(FM_MEMALL_ERR);
	sprintf(outfMITIFF_class,"%s/%s-%s_%s_%04d%02d%02d-%dhours_%s.mitiff",
		path_outf,pref_outf,pref_cl,arealist[tile],snowprod.h.year,
		snowprod.h.month,snowprod.h.day,period,satstring);
	ret = store_snow(outfMITIFF_class, snowprod.h, class, 
		CLASSLIMITS, class_desc,satstring);
	if (ret != 0)  { 
	    fmerrmsg(where,"Could not create MITIFF file %s", outfMITIFF_class);
	    exit(FM_IO_ERR);
	} 


	/* create mitiff for classified ice probability image */
	outfMITIFF_psnow = (char *) malloc(FILELEN+5);
	if (!outfMITIFF_psnow) exit(FM_MEMALL_ERR);
	sprintf(outfMITIFF_psnow,"%s/%s-%s_%s_%04d%02d%02d-%dhours_%s.mitiff",
		path_outf,pref_outf,pref_ps,arealist[tile],snowprod.h.year, 
		snowprod.h.month,snowprod.h.day,period,satstring);
	ret = store_snow(outfMITIFF_psnow, snowprod.h, snowclass, 
		PROBLIMITS, snow_desc,satstring);
	if (ret != 0)  {
	    fmerrmsg(where,"Could not create MITIFF file %s", outfMITIFF_psnow);
	    exit(FM_IO_ERR);
	}


	fprintf(stdout,"\t%s\n",outfHDF);
	fprintf(stdout,"\t%s\n",outfMITIFF_class);
	fprintf(stdout,"\t%s\n\n",outfMITIFF_psnow);


	/* Deallocate memory */
	free(outfHDF); 
	free(outfMITIFF_class);
	free(outfMITIFF_psnow);
	free(class);
	free(snowclass);
	free(probsnow);
	free(probclear);
	free_osihdf(&snowprod);
    } /*end for-loop for tile*/

    /* Free more memory */
    free(date_start);
    if (aflg) { free(pref_outf); }
    free(path_outf);
    free(satstring);

    for (i=0;i<numsat;i++) { free(satlist[i]); }
    free(satlist);

    if (mflg) {
	for (i=0;i<numarea;i++) { free(arealist[i]); }
	free(arealist);
    }

    free(num_files_area);
    free(num_files_area_counter);
    for (f=0;f<nrInput;f++) {
	free(infile_avhrrice[f]);
	free(infile_sorted[f]);
    }

    free(infile_avhrrice);
    free(infile_sorted);

    if (lflg) {free(satlistfile);}
    if (mflg) {free(arealistfile);}

    fprintf(stdout,"\t=================================================\n");


    exit(FM_OK);
}

void usage() {
    fprintf(stdout,"\n  SYNTAX: \n");
    fprintf(stdout,
	    "  accusnow -s <dir_avhrrice> -d <date_end> -p <period>\n");
    fprintf(stdout,"\t  -a <pref_outf> -o <path_outf> (-t <satellite name>\n");
    fprintf(stdout,
	    "\t  -l <satlist> -c <cloudlimit> -m <arealist> -z) \n\n");
    fprintf(stdout,
	    "  <dir_avhrrice> : Directory with hdf5 files with avhrr ice data.\n");
    fprintf(stdout,"  <date_end>   : End date of merging period\n");
    fprintf(stdout,"  <period>       : Merging period in hours.\n");
    fprintf(stdout,"  <pref_outf>    : Prefix for output filename ");
    fprintf(stdout,"(currently not in use!).\n");
    fprintf(stdout,"  <path_outf>    : Path for output file.\n");
    fprintf(stdout,"  <satellite>    : Name of one satellite ");
    fprintf(stdout,"if processing only one (optional).\n");
    fprintf(stdout,
	    "  <satlist>      : File with satellites to use (optional).\n");
    fprintf(stdout,"  <arealist>     : File with tile areas to use.\n");
    fprintf(stdout,
	    "  <cloudlimit>   : Probability limit for class cloud (optional).\n");
    fprintf(stdout,"  -z             : Use threshold on satellite ");
    fprintf(stdout,"zenith angle (not in use!).\n\n");
    exit(FM_OK);
}
