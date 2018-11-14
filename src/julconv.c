/*        Function "jday"
 *        Converts a Gregorian calendar date to the
 *        appropriate Buoy Group Julian calendar date.
 *        Version 1.0        Feb., 1989        R. Goldsmith 
 */
/*        The whole thing starts off 1200 Jan 01, 4713 BC.   */

#include  <stdio.h>
#include <string.h>

#define   dpy4         1461        /* days in four years  */
#define   dpc4       146097        /* days in four centuries  */
#define   jdbase    1721119        /* Julian base days BC. */

jday (day, mon, yearc, julian)


int       day;            /* input day of the month.  */
int       mon;            /* input month of the year.  */
int       yearc;          /* input year including century.  */
int      *julian;         /* returned Julian day number passed to routine.   */

{
int       d, m, y, c;

d = day;
y = yearc;
if (mon > 2) 
  m = mon - 3;
 else  {
  m = mon + 9;
  y = y - 1;
  }

c = y/100;
y = y - c*100;
*julian = (dpc4*c)/4 + (dpy4*y)/4 + (153*m + 2)/5 + d + jdbase;

return;
}


/* ==================================================================== */
/*        Function "jdate".
 *        Convert Julian day numbers to corresponding Gregorian calendar
 *        dates.   Converted to C from VAX BUOY Group routine.
 *        Version 1.0         Feb., 1989        R. Goldsmith 
 */
 
jdate (julian, day, mon, yearc)

int       julian;         /* Julian day number passed to routine.   */
int      *day;            /* returned day of the month.  */
int      *mon;            /* returned month of the year.  */
int      *yearc;          /* returned year including century.  */

{
int       j, in;
int       d, m, y;


j = julian - jdbase;          /* Get things back to about 2000 bc.  */
in = 4*j - 1;
y = in/dpc4;
j = in - dpc4*y;
in = j/4;
in = 4*in +3;
j = in/1461;
d = ((in - dpy4*j) + 4)/4;
in = 5*d - 3;
m = in/153;
d = ((in - 153*m) + 5)/5;
y = y*100 + j;
if (m < 10)  
  m = m + 3;
 else  {
  m = m - 9;
  y = y + 1;
  }

*day = d;
*mon = m;
*yearc = y;

return;
}


/* ==================================================================== */
/* Function gdate2dday convert the Gregorian date/time to the decimal
 * Julian day and fraction thereof.
 */

gdate2dday (gdate, gtime, dday)
int     gdate[3], gtime[3];   /* YMD HMS format */
double  *dday;

{
int     tem_jday;

double  dtem;


     /* Compute the Julian day. */
jday (gdate[2], gdate[1], gdate[0], &tem_jday);
     /* Compute the fraction. */
hms2dd (gtime[0], gtime[1], gtime[2], &dtem);

*dday = (double) tem_jday + dtem;
return (0);
}


/* =================================================================*/
/*   Program date_c0  is a generic date string conversion routine.
/*   It makes some assumptions about the format of a date string and
/*   then calls something relevant to parse it into a three element
/*   integer array.
 */

int date_c0 (cdate, idate)
char    cdate[];
int     idate[];


{
int     lenc;
int     iopstat;


lenc = strlen (cdate);

     /* Is it the 971231 form? */
if (lenc == 6)  {
  iopstat = date_c6 (cdate, idate);
  }

     /* Is it the 31Dec97 form? */
 else if (lenc == 7)  {
  iopstat = date_c7 (cdate, idate);
  }

     /* Is it the 19971231 form? */
 else if (lenc == 8)  {
  iopstat = date_c8 (cdate, idate);
  }

     /* Is it the 31Dec1997 form? */
 else if (lenc == 9)  {
  iopstat = date_c9 (cdate, idate);
  }

return (iopstat);
}


/* ==================================================================== */
/*   Program date_c6  takes a date string of the form 971231
/*   and returns a three element integer array.
 */

int date_c6 (cdate, idate)
char    cdate[];
int     idate[];


{
int     iopstat;
int     iy, im, id;


iopstat = sscanf (cdate, "%2d%2d%2d", &iy, &im, &id);
if (iopstat < 3)  {
  fprintf (stderr, "  *** date_c6 error:  converting date string %s\n",
           cdate);
  return (-1);
  }

if (iy > 70  && iy < 100)  iy = iy + 1900;
 else  iy = iy + 2000;

     /* Check for date validity. */
iopstat = date_chk (iy, im, id);
if (iopstat != 0)  {
  fprintf (stderr, "  *** date_c6:  error validating date %d %d %d\n",
           iy, im, id);
  return (iopstat);
  }

idate[0] = iy;
idate[1] = im;
idate[2] = id;
iopstat = 0;

return (iopstat);
}


/* ==================================================================== */
/*   Program date_c7  takes a date string of the form 31Dec97  
/*   and returns a three element integer array.
 */

int date_c7 (cdate, idate)
char    cdate[];
int     idate[];


{
char    cm[4];
char    cMON[] = {"JAN:FEB:MAR:APR:MAY:JUN:JUL:AUG:SEP:OCT:NOV:DEC\0"};

int     i, iopstat;
int     iy, im, id;


iopstat = sscanf (cdate, "%2d%3s%2d", &id, cm, &iy);
if (iopstat < 3)  {
  fprintf (stderr, "  *** date_c7 error:  converting date string %s\n",
           cdate);
  return (-1);
  }
 
im = 0;
for (i=0; i<strlen (cMON); i=i+4)  {
  if (strncasecmp (cm, cMON+i, 3) == 0)  {
    im = i/4 + 1;
    break;
    }
  }

if (iy > 70  && iy < 100)  iy = iy + 1900;
 else  iy = iy + 2000;

     /* Check for date validity. */
iopstat = date_chk (iy, im, id);
if (iopstat != 0)  {
  fprintf (stderr, "  *** date_c7:  error validating date %d %d %d\n",
           iy, im, id);
  return (iopstat);
  }

idate[0] = iy;
idate[1] = im;
idate[2] = id;
iopstat = 0;

return (iopstat);
}


/* ==================================================================== */
/*   Program date_c8  takes a date string of the form 19971231
/*   and returns a three element integer array.
 */

int date_c8 (cdate, idate)
char    cdate[];
int     idate[];


{
int     iopstat;
int     iy, im, id;


iopstat = sscanf (cdate, "%4d%2d%2d", &iy, &im, &id);
if (iopstat < 3)  {
  fprintf (stderr, "  *** date_c8 error:  converting date string %s\n",
           cdate);
  return (-1);
  }

     /* Check for date validity. */
iopstat = date_chk (iy, im, id);
if (iopstat != 0)  {
  fprintf (stderr, "  *** date_c8:  error validating date %d %d %d\n",
           iy, im, id);
  return (iopstat);
  }

idate[0] = iy;
idate[1] = im;
idate[2] = id;
iopstat = 0;

return (iopstat);
}


/* ==================================================================== */
/*   Program date_c9  takes a date string of the form 31Dec1997
/*   and returns a three element integer array.
 */

int date_c9 (cdate, idate)
char    cdate[];
int     idate[];


{
char    cm[4];
char    cMON[] = {"JAN:FEB:MAR:APR:MAY:JUN:JUL:AUG:SEP:OCT:NOV:DEC\0"};

int     i, iopstat;
int     iy, im, id;


iopstat = sscanf (cdate, "%2d%3s%4d", &id, cm, &iy);
if (iopstat < 3)  {
  fprintf (stderr, "  *** date_c9 error:  converting date string %s\n",
           cdate);
  return (-1);
  }
 
im = 0;
for (i=0; i<strlen (cMON); i=i+4)  {
  if (strncasecmp (cm, cMON+i, 3) == 0)  {
    im = i/4 + 1;
    break;
    }
  }

     /* Check for date validity. */
iopstat = date_chk (iy, im, id);
if (iopstat != 0)  {
  fprintf (stderr, "  *** date_c9:  error validating date %d %d %d\n",
           iy, im, id);
  return (iopstat);
  }

idate[0] = iy;
idate[1] = im;
idate[2] = id;
iopstat = 0;
return (iopstat);
}


/* ==================================================================== */
/*   Check validity of integer date array values.
 */

int  date_chk (iy, im, id)

int     iy, im, id;
 

{
int     md[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};


if (iy < 1900)  return (-1);
if (im < 1  ||  im > 12)  return (-2);
if (id < 1  ||  id > md[im-1])  {
  if (im == 2  &&  id == 29  &&  (iy % 4) == 0)  {
    if (iy == 2000)  return (-3);
    }
   else  return (-3);
  }

return (0);
}


/* ==================================================================== */
/* Function hms2dd converts integer {hours, minutes, seconds} to
 * decimal days.
 */

hms2dd (h, m, s, dd)
int     h, m, s;
double *dd;

{
*dd = ((double) h +
       (double) m/60.0 +
       (double) s/3600.0)/24.0;

return (0);
}


/* ==================================================================== */
/* Function dd2hms converts decimal days to 
 * integer {hours, minutes, seconds}
 */

dd2hms (dd, h, m, s)
double  dd;
int     *h, *m, *s;

{
int    ih, im , is;

double dtem;


dtem = dd*24.0*60.0*60.0 + 0.001;
ih = (int) dtem/3600;
im = (int) (dtem - (double) ih*3600.0)/60;
is = (dtem - (double) ih*3600.0 - (double) im*60.0);

*h = ih;
*m = im;
*s = is;

return (0);
}


