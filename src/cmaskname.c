/*
 * PURPOSE:
 * To generate filenames for cloud mask as nobody seems to take care of
 * consistency of such stupid things... 
 *
 * NOTES:
 *
 * BUGS:
 *
 * RETURN VALUES:
 * 0 - normal and correct ending
 * 1 - i/o problem
 * 2 - memory problem
 *
 * DEPENDENCIES:
 *
 * AUTHOR:
 * Øystein Godøy, DNMI/FOU, 20/09/2000
 * MODIFICATIONS:
 * Øystein Godøy, met.no/FOU, 27.09.2004
 * Modification for full Bayes approach started, that would imply that
 * this function is obsolete...
 *
 * CVS_ID:
 * $Id: cmaskname.c,v 1.2 2009-03-01 21:24:15 steingod Exp $
 */

#include <fmutil.h>
#include <fmio.h>

short cmaskname(char *rname, struct mihead info, char *cmname) {

    char sname[7];
    
    if (strncmp(info.satellite,"NOAA-14",7) == 0) {
	sprintf(sname,"%s","noaa14");
    } else if (strncmp(info.satellite,"NOAA-12",7) == 0) {
	sprintf(sname,"%s","noaa12");
    } else if (strncmp(info.satellite,"NOAA-15",7) == 0) {
	sprintf(sname,"%s","noaa15");
    } else if (strncmp(info.satellite,"NOAA-16",7) == 0) {
	sprintf(sname,"%s","noaa16");
    } else {
	return(1);
    }

    sprintf(cmname,"%4d%02d%02d%02d%02d.dn%c%c.%s.acmg.mitiff",
	info.year,info.month,info.day,info.hour,info.minute,
	rname[20],rname[21],sname);

    return(0);
}
