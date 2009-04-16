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
 * NAME: accusnow (working name)
 * 
 * PURPOSE: Average (AVHRR) ice product files to compound a cloudfree
 *          snowmap.
 *          Based on /osisaf/OSI_HL_SST/product/binsrc/average_merge_SST.c
 *          
 * 
 * SYNTAX: accusnow -s <dir_avhrrice> -d <date_start> 
 *         -p <period> -a <pref_outf> -o <path_outf>
 *         (-t <satellite> -l <satlist> -m <arealist> -z)
 *
 *    <dir_avhrrice> : Directory with hdf5 files with avhrr ice data.
 *    <date_start>   : Start date of merging period.
 *    <period>       : Merging period in hours.
 *    <pref_outf>    : Prefix for output filename.
 *    <path_outf>    : Path for output file.
 *    <satellite>    : Name of one satellite if processing only one (optional).
 *    <satlist>      : File with satellites to use (optional).
 *    <arealist>     : File with tile areas to use (optional).
 *    -z             : Use threshold on satellite zenith angle (value from header file).
 * 
 *
 * AUTHOR: Steinar Eastwood, DNMI, 21.08.2000
 * MODIFIED: 
 * SE, DNMI, 19.03.2001, 11.06.2001
 * SE, DNMI/FOU, 18.07.2001  Removed DNMI FELT file output.
 * SE, met.no, 26.08.2002  Change input arguments.
 * SE, met.no, 10.01.2003  Change reading of SST/QSST fields and
 *                         function average_merge_filesQF.
 * SE, met.no, 09.07.2003  Time SST to HDF file.
 * SE, met.no, 28.11.2003  New handling of input arguments and reading of
 *                         of input files with dirent.
 * SE, met.no, 09.12.2003  New merging function as standard.
 * SE, met.no, 13.01.2004  Introduced int main.
 * SE, met.no, 13.02.2004  Introduced -t option: to process only one satellite.
 * SE, met.no, 18.02.2004  Changed check on file name length.
 * SE, met.no, 24.02.2004  Bug in product time.
 * SE, met.no, 04.08.2004  No writing of grib files.
 * SE, met.no, 10.08.2004  HR SST files now contain SST, QFLG and SZA.
 * SE, met.no, 04.02.2005  Added use of satlist and arealist.
 * SE, met.no, 20.04.2005  Optional use of satellite senith angle in selecting
 *                         obs for averaging.
 * SE, met.no, 06.06.2005  Including proj string in HDF5 output header.
 * SE, met.no, 08.03.2006  Changed interface to sec1970toDateStr..
 * SE, met.no, 01.11.2006  No longer QSST and TSST HDF output, all fields in on file.
 * SE, met.no, 28.11.2006  Assume no ice if icefile is not provided.
 *
 * MAK,met.no, 08.01.2009  Using average_merge_SST.c as starting point 
 *                         for averaging files for OSI SAF avhrrice.
 */ 

 
#include <accusnow.h>
#include <unistd.h>

#define LISTLEN 100
#define DEFAULTCLOUD 0.2
#define TOTAREAS 8
#define MAXSAT 20
#define MAXAREA 20

void usage(); 


int main(int argc, char *argv[]) 
{    
  char *errmsg="\n\tERROR(accusnow):";
  extern char *optarg;
  char *dir_avhrrice, *date_start, *date_prod, *date_end;
  int sflg, dflg, pflg, aflg, oflg, tflg, lflg, mflg, zflg, cflg;
  int period, i, j, f, t, tile, nrInput, ret, ind, numf, ctr;
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
  char *prod_desc[AVHRRICEPROD_LEVELS] = {"class","P(snow)","P(clear)"};
  osi_dtype prod_ft[AVHRRICEPROD_LEVELS] = {CLASS_DT,PROB_DT,PROB_DT};
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
  fprintf(stdout,"\t|                  ACCUSNOW                     |\n");
  fprintf(stdout,"\t=================================================\n");
  fprintf(stdout,"\n");

  /* Interprete commandline arguments */
  sflg=dflg=pflg=aflg=oflg=tflg=lflg=mflg=zflg=cflg=0;
  while ((ret = getopt(argc, argv, "s:d:p:a:o:t:l:m:c:z")) != EOF) {
    switch (ret) {
    case 's':
      dir_avhrrice = (char *) malloc(strlen(optarg)+1);
      if (!strcpy(dir_avhrrice, optarg)) exit(1);
      sflg++;
      break;
    case 'd':
      date_start = (char *) malloc(strlen(optarg)+1);
      if (!strcpy(date_start, optarg)) exit(1);
      dflg++;
      break;
    case 'p':
      period = atoi(optarg);
      pflg++;
      break;
    case 'a':
      pref_outf = (char *) malloc(strlen(optarg)+1);
      if (!strcpy(pref_outf, optarg)) exit(1); 
      aflg++; 
      break;
    case 'o':
      path_outf = (char *) malloc(strlen(optarg)+1);
      if (!strcpy(path_outf, optarg)) exit(1);
      oflg++;
      break;
    case 't':
      procsat = (char *) malloc(strlen(optarg)+1);
      if (!strcpy(procsat, optarg)) exit(1);
      tflg++;
      break;
    case 'l':
      satlistfile = (char *) malloc(strlen(optarg)+1);
      if (!strcpy(satlistfile, optarg)) exit(1);
      lflg++;
      break;
    case 'm':
      arealistfile = (char *) malloc(strlen(optarg)+1);
      if (!strcpy(arealistfile, optarg)) exit(1);
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
    exit(0);
  }
  
 
  if (strlen(date_start) ==8 ) strcat(date_start,"00");
  if (strlen(date_start) != 10) {
    fprintf(stdout,
	    "\n ERROR: date_start has wrong format, yyyymmddhh needed\n\n");
    exit(0);
  }

  stime    = ymdh2fmsec1970(date_start,0);
  etime    = stime + period*3600;
  date_end = (char *) malloc(20+1);
  fmsec19702isodatetime(etime, date_end);
  prodtime = stime + period*3600/2;
  date_prod= (char *) malloc(20+1);
  fmsec19702isodatetime(prodtime, date_prod);
  ret      = tofmtime(prodtime,&timedate);

  ctr = 0;
  fprintf(stdout,"\n\tAVHRR ice files dir:   %s \n", dir_avhrrice);
  fprintf(stdout,"\tDate start:            %c", date_start[ctr]);
  for (ctr=1;ctr<8;ctr++) {fprintf(stdout,"%c",date_start[ctr]);}
  fprintf(stdout,"\n\tPeriod:                %d \n", period);
  fprintf(stdout,"\tDate end:              %s \n", date_end);
  fprintf(stdout,"\tDate product:          %s \n", date_prod);
  fprintf(stdout,"\tPeriod start:          %d \n", (int)stime);
  fprintf(stdout,"\tPeriod end:            %d \n", (int)etime);
  fprintf(stdout,"\tPath outfile:          %s \n", path_outf);
  if (aflg) {
    fprintf(stdout,"\tPrefix outfile:        %s \n", pref_outf);
  }
  else {
    pref_outf = defpref;
  }

  if (tflg) {
    fprintf(stdout,"\n\tProcess only satellite: %s \n", procsat);
  }
  if (lflg) {
  fprintf(stdout,"\n\tUsing list of satellites from file: %s \n", satlistfile);
  }
  if (!cflg) {
    cloudlim = DEFAULTCLOUD;
  }
  fprintf(stdout,"\n\tUsing cloud probability limit: %.1f \n", cloudlim);
  if (mflg) {
    fprintf(stdout,"\n\tUsing area tiles from file: %s \n", arealistfile);
  }
  
  satlist = (char **) malloc(MAXSAT*sizeof(char *)); /*flytte til etter if lfgl?*/
 
  /* Read satellite name list */
  numsat = 0;
  if (lflg) {
    numsat = read_sat_area_list(satlistfile,satlist);
    if (numsat == 0) {
      fprintf(stderr,"\n\tFound NO satellites on file %s\n\texiting..\n\n",
	      satlistfile);
    }
    else if (numsat < 0) {
      fprintf(stderr,"%s Could not read file %s\n",errmsg,satlistfile);
      fprintf(stderr,"\t       Processing all available satnames.\n");
      lflg = 0;
    }
    else {
      fprintf(stdout,"\n\tProcessing satellites (%d): \n\t  ",numsat);
      satstring = (char *) malloc(LISTLEN*sizeof(char));
      for (i=0;i<numsat;i++) { 
	fprintf(stdout,"%s ",satlist[i]);
	if (i != 0) strcat(satstring,"_");
	strcat(satstring,satlist[i]);
      }
      fprintf(stdout,"\n");
    }
  } 
  else if (tflg) {
    satstring = (char *) malloc(strlen(procsat)*sizeof(char));
    sprintf(satstring,"%s",procsat);
  }
  else { /* No sat.name list given, need componet for filename */
    satstring = (char *) malloc(strlen("allsats")*sizeof(char));
    sprintf(satstring,"allsats");
    fprintf(stdout,"\n\tProcessing all available satnames.\n");
  }

  numarea = 0;
  arealist = (char **) malloc(MAXAREA*sizeof(char *));
  if (mflg) {/* Read tile list */
    numarea = read_sat_area_list(arealistfile,arealist);
    if (numarea == 0) {
      fprintf(stderr,"\n\tWARNING! Found NO area tiles on file %s\n",
	      arealistfile);
      exit(0);
    }
    else if (numarea < 0) {
      fprintf(stderr,"%s Could not read file %s\n",errmsg,arealistfile);
      fprintf(stderr,
	      "\t       Processing all available area tile names instead.\n");
      mflg = 0;
    }
    else {
      fprintf(stdout,"\n\tProcessing area tiles (%d): \n\t  ",numarea);
      for (i=0;i<numarea;i++) { 
	fprintf(stdout,"%s ",arealist[i]);
      }
      fprintf(stdout,"\n");
    }
  }
  if (mflg == 0) { /*hard code a list of all known areas*/
    numarea = TOTAREAS;
    arealist = default_arealist;
    fprintf(stdout,
    "\n\tProcessing files for tiles on default tilelist (%d):\n\t ",numarea);
    for (i=0;i<numarea;i++) {
      printf("%s ",arealist[i]);
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n");


  /* Allocate and init. variables to handle sorting of files based on tile*/
  num_files_area = (int *)malloc(numarea*sizeof(int));
  num_files_area_counter = (int *)malloc(numarea*sizeof(int));
  for (i=0;i<numarea;i++){
    num_files_area[i] = num_files_area_counter[i] = 0;
  }


  /* Reading directory, find all avhrrice files to process. */
  dirp_avhrrice = opendir(dir_avhrrice);
  if (!dirp_avhrrice) {
    fprintf(stderr,"%s Could not open dir %s\n",errmsg,dir_avhrrice);
    return(10);
  }

  numf = 0;
  while ((dirl_avhrrice = readdir(dirp_avhrrice)) != NULL) {
    if (strncmp(dirl_avhrrice->d_name,BASENMAVHRRICE,strlen(BASENMAVHRRICE)) == 0 &&
	strstr(dirl_avhrrice->d_name,".hdf") != NULL && strlen(dirl_avhrrice->d_name) >= MINLENNMAVHRRICE) {
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
  infile_sorted   = (char **) malloc(numf*sizeof(char *));
  i = 0;

  /*Loop through files in input directory*/
  while ((dirl_avhrrice = readdir(dirp_avhrrice)) != NULL) {

    /*Check input filename*/
    if (strncmp(dirl_avhrrice->d_name,BASENMAVHRRICE,strlen(BASENMAVHRRICE)) 
	== 0 &&	strstr(dirl_avhrrice->d_name,".hdf") != NULL && 
	strlen(dirl_avhrrice->d_name) >= MINLENNMAVHRRICE) {

      /*Check time of file*/
      sret = strncpy(datestr,&dirl_avhrrice->d_name[10],12);
      datestr[12] = '\0';
      sprintf(datestr_ymdhms,"%s00",datestr);
      ftime = ymdhms2fmsec1970(datestr_ymdhms,0);
      if (ftime < stime || ftime > etime) {
	continue;
      }

      /*Check that file can be opened (this removes files of size zero!*/
      checkfile = (char *) malloc(256*sizeof(char));
      sprintf(checkfile,"%s/%s",dir_avhrrice,dirl_avhrrice->d_name);
      init_osihdf(&checkfileheader);
      ret = read_hdf5_product(checkfile,&checkfileheader,1);/*header only*/
      free(checkfile);
      if (ret != 0) {
	fprintf(stderr,"%s while trying to open file %s, \n",
		errmsg,dirl_avhrrice->d_name);
	fprintf(stderr,"\t skipping file.\n\n");
	free_osihdf(&checkfileheader);
	continue;
      }

      /*Check satellite name*/
      if (lflg || tflg) { /* satname is not part of input filename, must */
	satfound = 0;     /* check the header to see if file should be kept */
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
	  /*fprintf(stdout,"\tFile not in satname list, skipping %s\n",
	    dirl_avhrrice->d_name);*/
	  continue;
	}
      }
      /* else { /\*avhrrice currently not nice for metop, removing metop!*\/ */
      /* checkfile = (char *) malloc(256*sizeof(char)); */
      /* sprintf(checkfile,"%s/%s",dir_avhrrice,dirl_avhrrice->d_name); */
      /* init_osihdf(&checkfileheader); */
      /* ret = read_hdf5_product(checkfile,&checkfileheader,1);/\*header only*\/ */
      /* free(checkfile); */
      /* if (strcmp(checkfileheader.h.source,"NOAA-p02")== 0){ */
      /*   fprintf(stdout,"\tSkipping METOP for now..(%s)\n", */
      /* 	  dirl_avhrrice->d_name); */
      /*  continue; */
      /*  } */
      /* } /\*this else-part can be removed when metop files are ok*\/ */

   
      /*Check tile, this returns the element number*/
      ind = find_sat_area_index(arealist,numarea,dirl_avhrrice->d_name);
      if (ind < 0) {
	fprintf(stdout,"\tFile not in area tile list, skipping %s\n",
		dirl_avhrrice->d_name);
	continue;
      }
      
	      
      /*If all tests passed, add to list of ok files :-)*/
      infile_avhrrice[i] = (char *) malloc(256*sizeof(char));
      infile_sorted[i]   = (char *) malloc(256*sizeof(char));
      sprintf(infile_avhrrice[i],"%s/%s",dir_avhrrice,dirl_avhrrice->d_name);
      i ++;
      /*MAK adding this here*/
      num_files_area[ind]++;
    }
  }
  free(dir_avhrrice);
  nrInput = i;

  if (nrInput == 0) {
    fprintf(stdout,"\n\tNo files to be processed, exiting\n\n");
    return(0);
  }

  fprintf(stdout,"\n\tFiles to be processed: (%d,%d) \n",nrInput,numf);
  for (i=0;i<nrInput;i++) {
    fprintf(stdout,"\t%2d %s\n",i,infile_avhrrice[i]);
  }

  fprintf(stdout,"\n\tSorting input files..\n");
  for (t=0;t<numarea;t++) {
    fprintf(stdout,"\t Tile: %s  Number of files: %d\n",
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


  if (timedate.fm_hour > 23 || timedate.fm_min > 59 || timedate.fm_mday > 31 || timedate.fm_mon > 12 || 
      timedate.fm_year > 3000 || timedate.fm_year < 1978)  {
    fprintf(stderr,"%s In input time/date. Exiting. \n",errmsg);
    exit(2);
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

    /*creating list of files for this tile*/
    infile_currenttile = (char **) malloc(num_files_area[tile]*sizeof(char *));
    for (f=0;f<num_files_area[tile];f++) {
      infile_currenttile[f] = (char *) malloc(256*sizeof(char));
      sprintf(infile_currenttile[f],infile_sorted[f+index_offset]);
    }
    /*for (f=0;f<num_files_area[tile];f++){
      fprintf(stdout,"%s\n",infile_currenttile[f]);
      }*/




    inputhdf = (osihdf *) malloc(num_files_area[tile]*sizeof(osihdf));

    init_osihdf(&snowprod); /*Dette init. kun headeren, ikke datafeltet!*/

    /*reading headers of all files on list to compare + collect info*/
    for (f=0;f<num_files_area[tile];f++) {
      init_osihdf(&inputhdf[f]);
      ret=read_hdf5_product(infile_currenttile[f],&inputhdf[f], 1);
      if (ret != 0) {
	fprintf(stderr,"%s while trying to read header of %s \n",errmsg,
		infile_currenttile[f]);
      }
      if (inputhdf[f].h.iw != inputhdf[0].h.iw ||
	  inputhdf[f].h.ih != inputhdf[0].h.ih ||
	  inputhdf[f].h.Ax != inputhdf[0].h.Ax ||
	  inputhdf[f].h.Ay != inputhdf[0].h.Ay ||
	  inputhdf[f].h.Bx != inputhdf[0].h.Bx ||
	  inputhdf[f].h.By != inputhdf[0].h.By ) {
	fprintf(stderr,
		"%s Input files for tile %s are from different tiles!?\n",
		errmsg,arealist[tile]);
	exit(2);
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
    snowclass = (unsigned char *) malloc(safucs.iw*safucs.ih*sizeof(char));
    probsnow = (float *) malloc(safucs.iw*safucs.ih*sizeof(float));
    probclear = (float *) malloc(safucs.iw*safucs.ih*sizeof(float));
    
    if (!class || !snowclass || !probsnow || !probclear) {
      fprintf(stderr,"%s Could not allocate memory.\n",errmsg);
      exit(3);
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
	fprintf(stderr,"%s while running \"average_merge_files\" \n",errmsg);
	exit(10);
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
      fprintf(stderr,"%s malloc_osisaf encountered an error.\n",errmsg);
      return(3); /*rappa fra func_sst_hdf, annen return value?*/
    }

    for (i=0;i<safucs.iw*safucs.ih;i++){
      ((int*)snowprod.d[0].data)[i] = class[i];
      ((float*)snowprod.d[1].data)[i] = probsnow[i];
      ((float*)snowprod.d[2].data)[i] = probclear[i];
    }
    
    

    fprintf(stdout,"\tCreating output files for tile %s:\n",arealist[tile]);

    /* create hdf product file */
    outfHDF = (char *) malloc(FILELEN+5);
    if (!outfHDF) exit(3);
    sprintf(outfHDF,"%s/%s_%s_%04d%02d%02d-%dhours_%s.hdf",path_outf,
	    pref_outf,arealist[tile],snowprod.h.year,snowprod.h.month,
	    snowprod.h.day ,period,satstring);
    ret = store_hdf5_product(outfHDF, snowprod);
    if (ret != 0)  {
      fprintf(stderr,"%s while creating HDF file.\n",errmsg);
      exit(10);
    }    
    

    /* create mitiff for classified image */
    outfMITIFF_class = (char *) malloc(FILELEN+5);
    if (!outfMITIFF_class) exit(3);
    sprintf(outfMITIFF_class,"%s/%s-%s_%s_%04d%02d%02d-%dhours_%s.mitiff",
	    path_outf,pref_outf,pref_cl,arealist[tile],snowprod.h.year,
	    snowprod.h.month,snowprod.h.day,period,satstring);
    ret = store_snow(outfMITIFF_class, snowprod.h, class, 
			       CLASSLIMITS, class_desc,satstring);
    if (ret != 0)  { 
      fprintf(stderr,"%s while creating MITIFF file.\n",errmsg);
      exit(10); 
    } 
    

    /* create mitiff for classified ice probability image */
    outfMITIFF_psnow = (char *) malloc(FILELEN+5);
    if (!outfMITIFF_psnow) exit(3);
    sprintf(outfMITIFF_psnow,"%s/%s-%s_%s_%04d%02d%02d-%dhours_%s.mitiff",
	    path_outf,pref_outf,pref_ps,arealist[tile],snowprod.h.year, 
	    snowprod.h.month,snowprod.h.day,period,satstring);
    ret = store_snow(outfMITIFF_psnow, snowprod.h, snowclass, 
			       PROBLIMITS, snow_desc,satstring);
    if (ret != 0)  {
      fprintf(stderr,"%s while creating MITIFF file.\n",errmsg);
      exit(10);
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




  /* Deallocate more memory */
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


  exit(0);
}




void usage() 
{
    fprintf(stdout,"\n  SYNTAX: \n");
    fprintf(stdout,
    "  accusnow -s <dir_avhrrice> -d <date_start> -p <period>\n");
    fprintf(stdout,"\t  -a <pref_outf> -o <path_outf> (-t <satellite name>\n");
    fprintf(stdout,
    "\t  -l <satlist> -c <cloudlimit> -m <arealist> -z) \n\n");
    fprintf(stdout,
    "  <dir_avhrrice> : Directory with hdf5 files with avhrr ice data.\n");
    fprintf(stdout,"  <date_start>   : Start date of merging period\n");
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
    exit(1);
}
