/*
 * NAME:
 * nwp_read.c
 *
 * PURPOSE:
 * To provide an interface between libusenwp and libfmcol.
 *
 * REQUIREMENTS:
 * o libsatimg
 * o libusenwp
 *
 * INPUT:
 * o filenames to search (not used at present)
 * o requested time
 * o requested output grid (as specified by remote sensing software)
 *
 * OUTPUT:
 * data structure
 * return values:
 * 0-OK
 * 1-command syntax error
 * 2-I/O problem
 * 3-Memory trouble
 *
 * NOTES:
 * Filenames are hardcoded at present, input filenames are not used due to
 * constraints in the C/Fortran interface (all filenames have to be equal
 * lengths) and the initialisation has to be predefined. Dynamic
 * allocation should be implemented but it needs to be contiguous in
 * memory...
 *
 * Run nwpice_init before nwpice_read!!
 *
 * BUGS:
 * NA
 *
 * AUTHOR:
 * �ystein God�y, met.no/FOU, 18.10.2004 
 *
 * MODIFIED:
 * NA
 *
 * CVS_ID:
 * $Id: getnwp.c,v 1.2 2009-03-01 21:24:15 steingod Exp $
 */

#include <getnwp.h>

int nwpice_read(char *filenames, fmtime reqtime, fmucsref refucs, 
	nwpice *nwp) {

    char *where="nwpice_read";
    int i, imgsize;
    /* HIRLAM/NWP variables */
    int nfiles=FFNS, iunit=10, interp=1, itime[5];
    int nparam1=NOFIELDS1, nparam2=NOFIELDS2;
    int ierror, iundef, len1, len2;
    float satgrid[10], *nwpfield1, *nwpfield2;
    /*
    char nwpfiles[FFNS][28];
    */
    /* Files to access NWP data at surface lev */
    static char feltfilepresl[FFNS][30] = { 
	    "/opdata/hirlam20/grdn00.dat", 
	    "/opdata/hirlam20/grdn06.dat", 
	    "/opdata/hirlam20/grdn12.dat", 
	    "/opdata/hirlam20/grdn18.dat" 
    };
    /* Files to access NWP data at model lev*/
    static char feltfilemodl[FFNS][31] = { 
	    "/opdata/hirlam20/grdnm00.dat", 
	    "/opdata/hirlam20/grdnm06.dat", 
	    "/opdata/hirlam20/grdnm12.dat", 
	    "/opdata/hirlam20/grdnm18.dat" 
    };
    /*
    enum nwp_par {T0M,T2M,T950,T800,T750,T500,MSLP,PW,TOPO};
    */
    int iparam1[NOFIELDS1][4]={ /* Spec of what to get */
	{ 30, 2, 1000, 0},  /* Temp. at 0m */
    }; 
    /*
    int iparam2[NOFIELDS2][4]={ 
	{  0, 0,  0, 7}   
    }; 
    */
    int icontrol[6]={ /* Accep. spec */
	    88,  /* Producer 88 -> MI */
	-32767,  /* Model grid -32767-> first found */
	     3,  /* Min allowed forecast length in hours */
	   +24,  /* Max allowed forecast length in hours */
	    -3,  /* Min allowed offset in hours from image time */
	    +3   /* Max allowed offset in hours from image time */
    };
    fmtime nwptime;

    /* Requested grid description, only one supported yet */
    satgrid[0] = 60.;
    satgrid[1] = 0.;
    satgrid[2] = 1000.;
    satgrid[3] = 1000.;
    satgrid[4] = 0.;
    satgrid[5] = 0.;
    satgrid[6] = refucs.Ax;
    satgrid[7] = refucs.Ay;
    satgrid[8] = refucs.Bx;
    satgrid[9] = refucs.By;
    /* Requested valid time */
    itime[0] = reqtime.fm_year;
    itime[1] = reqtime.fm_mon;
    itime[2] = reqtime.fm_mday;
    itime[3] = reqtime.fm_hour;
    itime[4] = reqtime.fm_min;

    /*
     * The fields are returned in the vector "field[MAXIMGSIZE]" which is
     * organized like this: "field[nx*ny*nparam]", where nx and ny are
     * the numbers of gridpoints in x and y direction and nparam is the
     * number of fields/parameters collected. 
     */
    printf(" Collecting HIRLAM data ");
    printf("and interpolating to satellite grid...\n"); 
    len1 =  strlen(feltfilepresl[0])+1;
    len2 =  strlen(feltfilemodl[0])+1;
    printf(" len1: %d len2: %d\n", len1, len2);
    /*
    len1 =  strlen(nwpfiles[0])+strlen(cfg.nwp_p)+2;
    for (j=0; j<FFNS;j++) {
	nwpfiles[j] = (char *) malloc(FILENAMELEN);
	if (!nwpfiles[j]) {
	    error(errmsg,"memory allocation feltfile");
	    return(3);
	}
	sprintf(nwpfiles[j],"%s/%s",cfg.nwp_p,feltfile[j]);
	printf(" %s\n", nwpfiles[j]);
    }
    */
    imgsize = refucs.iw*refucs.ih;
    if (imgsize > FMIO_MAXIMGSIZE) {
	fprintf(stderr," ERROR: Not enough space to hold NWP data.\n");
	fprintf(stderr," Required image size: %d, available size: %d\n", 
		imgsize, FMIO_MAXIMGSIZE);
	return(2);
    }
    nwpfield1 = (float *) malloc(NOFIELDS1*imgsize*sizeof(float));
    if (!nwpfield1) {
	printf(" Could not allocate nwpfield1..\n");
	return(3);
    }
    /*
    nwpfield2 = (float *) malloc(NOFIELDS2*imgsize*sizeof(float));
    if (!nwpfield2) {
	printf(" Could not allocate nwpfield2..\n");
	return(3);
    }
    */

    /*
     * Collect data from the files containing surface, pressure and
     * parameter fields
     */
    /*
    for (i=0;i<FFNS;i++) {
	printf(" feltfilepresl %s\n",feltfilepresl[i]);
    }
    */
    printf(" First parameter, surface and pressure fields are searched...\n");
    getfield_(&nfiles, feltfilepresl[0], &iunit, &interp, 
	    satgrid, &refucs.iw, &refucs.ih, itime, 
	    &nparam1, iparam1, icontrol, nwpfield1, &iundef, &ierror, len1);

    if (ierror) {
	fmerrmsg(where,"error reading surface and parameter fields");
	return(2);
    }

    /*
     * Collect data from the files containing model level fields
     */
    /*
    for (i=0;i<FFNS;i++) {
	printf(" feltfilemodl %s\n",feltfilemodl[i]);
    }
    */
    /* Not needed...
    printf(" Then model fields are searched...\n");
    getfield_(&nfiles, feltfilemodl[0], &iunit, &interp, 
	    satgrid, &refucs.iw, &refucs.ih, itime, 
	    &nparam2, iparam2, icontrol, nwpfield2, &iundef, &ierror, len2);

    if (ierror) {
	errmsg(where,"error reading model level fields");
	return(2);
    }
    */

    /*
     * Print status information about the model field which was found.
     * if "iundef" is different from 0, undefined values were found
     * in the HIRLAM fields. The code for undefined is +1.e+35.
     */
    printf(" Date: %02d/%02d/%d", 
	    itime[2], itime[1], itime[0]);
    printf(" and time: %02d:00 UTC\n", itime[3]);
    printf(" Forecast length: %d\n", itime[4]);
    printf(" Status of 'getfield':");
    printf(" ierror=%d iundef=%d\n\n", ierror, iundef);

    /*
     * Transfer the time information to the nwpice structure
     */
    nwp->leadtime = itime[4];
    nwptime.fm_year = itime[0];
    nwptime.fm_mon = itime[1];
    nwptime.fm_mday = itime[2];
    nwptime.fm_hour = itime[3];
    nwptime.fm_min = 0;
    nwptime.fm_sec = 0;
    nwp->validtime = tofmsec1970(nwptime);

    /*
     * Transfer the UCS information
     */
    nwp->refucs = refucs;

    /*
     * Allocate the data structure that will contain the NWP data
     */
    /*
    if (!nwp->t950hpa) {
	nwp->t950hpa = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->t950hpa) return(3);
    }
    if (!nwp->t800hpa) {
	nwp->t800hpa = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->t800hpa) return(3);
    }
    if (!nwp->t700hpa) {
	nwp->t700hpa = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->t700hpa) return(3);
    }
    if (!nwp->t500hpa) {
	nwp->t500hpa = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->t500hpa) return(3);
    }
    if (!nwp->t2m) {
	nwp->t2m = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->t2m) return(3);
    }
    */
    if (!nwp->t0m) {
	nwp->t0m = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->t0m) return(3);
    }
    /*
    if (!nwp->ps) {
	nwp->ps = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->ps) return(3);
    }
    if (!nwp->pw) {
	nwp->pw = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->pw) return(3);
    }
    if (!nwp->rh) {
	nwp->rh = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->rh) return(3);
    }
    if (!nwp->topo) {
	nwp->topo = (float *) malloc(imgsize*sizeof(float));
	if (!nwp->topo) return(3);
    }
    */

    /*
     * Transfer the NWP data from the local temporary array to the nwpice
     * structure
     */
    for (i=0;i<imgsize;i++) {
	/*
	 * Surface, parameter and pressure fields
	 */
	nwp->t0m[i] = nwpfield1[0*imgsize+i];
	/*
	nwp->t2m[i] = nwpfield1[1*imgsize+i];
	nwp->t950hpa[i] = nwpfield1[2*imgsize+i];
	nwp->t800hpa[i] = nwpfield1[3*imgsize+i];
	nwp->t700hpa[i] = nwpfield1[4*imgsize+i];
	nwp->t500hpa[i] = nwpfield1[5*imgsize+i];
	nwp->ps[i] = nwpfield1[6*imgsize+i];
	nwp->topo[i] = nwpfield1[7*imgsize+i];
	nwp->rh[i] = nwpfield1[8*imgsize+i];
	*/
	
	/*
	 * Model level fields
	 */
	/*
	nwp->pw[i] = nwpfield2[0*imgsize+i];
	*/
    }

    /*
     * Free allocated memory
     */
    free(nwpfield1);
    /*
    free(nwpfield2);
    */

    return(0);
}

/*
 * NAME:
 * nwpice_init
 *
 * PURPOSE:
 * To initialize the data structure and ease memory management.
 *
 * REQUIREMENTS:
 *
 * INPUT:
 *
 * OUTPUT:
 *
 * NOTES:
 *
 * BUGS:
 *
 * AUTHOR:
 * �ystein God�y, met.no/FOU, 21.10.2004 
 *
 * MODIFIED:
 * NA
 */

int nwpice_init(nwpice *nwp) {

    nwp->t950hpa = NULL;
    nwp->t800hpa = NULL;
    nwp->t700hpa = NULL;
    nwp->t500hpa = NULL;

    nwp->t0m = NULL;
    nwp->t2m = NULL;

    nwp->ps = NULL;
    nwp->pw = NULL;
    nwp->rh = NULL;

    nwp->topo = NULL;

    return(0);
}

/*
 * NAME:
 * nwpice_free
 *
 * PURPOSE:
 * To free resources allocated within nwpice_read.
 *
 * REQUIREMENTS:
 *
 * INPUT:
 *
 * OUTPUT:
 *
 * NOTES:
 *
 * BUGS:
 *
 * AUTHOR:
 * �ystein God�y, met.no/FOU, 21.10.2004 
 *
 * MODIFIED:
 * NA
 */

int nwpice_free(nwpice *nwp) {

    if (nwp->t950hpa) free(nwp->t950hpa);
    if (nwp->t800hpa) free(nwp->t800hpa);
    if (nwp->t700hpa) free(nwp->t700hpa);
    if (nwp->t500hpa) free(nwp->t500hpa);

    if (nwp->t0m) free(nwp->t0m);
    if (nwp->t2m) free(nwp->t2m);

    if (nwp->ps) free(nwp->ps);
    if (nwp->pw) free(nwp->pw);
    if (nwp->rh) free(nwp->rh);

    if (nwp->topo) free(nwp->topo);

    return(0);
}
