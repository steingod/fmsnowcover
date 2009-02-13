/* *****************************************************************
 * COPYRIGHT:
 * EUMETSAT
 *
 * PRODUCED BY:
 * Norwegian Meteorological Institute (met.no)
 * Research and Development Department
 * P.O.BOX 43 - Blindern, N-0313 OSLO, NORWAY
 *
 * This SW was developed by met.no and the Danish Meteorological
 * Institute (DMI) within the context of the Co-operation Agreement
 * for the development of a pilot SAF on Ocean and Sea Ice.
 * *****************************************************************/

/*
 * PURPOSE:
 * To generate filenames for cloud mask as nobody seems to take care of
 * consistency of such stupid things... Belongs to OSI HL pap ice package.
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
 * $Id: cmaskname.c,v 1.1 2009-02-13 23:23:13 steingod Exp $
 */

#include <satimg.h>

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
