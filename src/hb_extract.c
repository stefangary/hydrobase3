/* hb_extract.c
................................................................................
                              *  HydroBase 3 *
................................................................................
                          
                    author:  Ruth Curry
                             Woods Hole Oceanographic Inst
                             1993
			     Updated to ANSI-C Feb 2000
			     Updated to HydroBase 3: Dec 2008
			     Supports multiple -E options Dec 2009
			     Supports variance,count,quality arrays 
................................................................................
................................................................................
.  Extracts stations from the HydroBase files by a variety of criteria 
.
................................................................................
................................................................................

*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "hydrobase.h"

#define    DIR     ""
#define    EXTENT   ""
#define    MAXEXT   10  /* maximum # of different file extents specified */
#define    MAXCRIT  12   /* maximum # of criteria for any of the options */
#define    PRINT_MSG 1

/* prototypes for locally defined functions */

   int     checkpos(struct HYDRO_HDR *, float, float, float, float);
   int     checkmonth(struct HYDRO_HDR *, int);
   int     checkyear(struct HYDRO_HDR *, int, int);
   int     checkship(struct HYDRO_HDR *, char **, int); 
   int     checkcruise(struct HYDRO_HDR *, int *, int); 
   int     checknation(struct HYDRO_HDR *, char **, int);
   int     checkdepth(struct HYDRO_HDR *, double *, int, int); 
   int     checkpdr(struct HYDRO_HDR *, int, int); 
   int     checklevs(struct HYDRO_HDR *, struct HYDRO_DATA *, int, int, float, float);
   int     checkprop(struct HYDRO_HDR *, int *, int, int, int);
   int     checkinstrument(struct HYDRO_HDR *, char *, int); 
   void    print_usage(char *);

main (int argc, char **argv)
{
   int     curfile = 1, nfiles = 0;
   int     curext, n_extents; 
   int     i,  m, index;
   int     g_flag=0, h_flag=0, m_flag=0, y_flag=0, d_flag=0, e_flag=0;
   int     n_flag=0, s_flag=0, c_flag=0, r_flag=0, i_flag=0;
   int     ld_flag=0, lp_flag=0, b_flag=0, xdateline = 0;
   int     p_flag=0, and_flag=0, or_flag=1;
   int     clist[MAXCRIT], ncru = 0;
   int     proplist[MAXCRIT], nprops = 0;
   int     ncountry=0, nship = 0, ninstr = 0;
   char    **ship, **country, instr[MAXCRIT], cprop[4];
   int     status, staOK, topt = 0;
   int     nout = 0, nskip = 0;
   int     mask, month, error = 0;
   int     minsta, maxsta;
   int     minyr, maxyr;
   int     mindepth, maxdepth;
   int     minpdr, maxpdr;
   FILE *infile, *rejectfile, *outfile;
   float   top, bot, left, right;
   float   minlev, maxlev;
   struct HYDRO_HDR h;
   struct HYDRO_DATA data;
   char   *dir, **extent_list;
   char   *s;

/* are there command line arguments? */

   if (argc < 2) {
      print_usage(argv[0]);
      exit(1);
   }
   
/* allocate space for lists */

    ship = (char **) calloc((size_t)MAXCRIT, (size_t) sizeof(char *));
    country = (char **) calloc((size_t)MAXCRIT, (size_t) sizeof(char *));
    for (i = 0; i < MAXCRIT; ++i) {
       ship[i] = (char *) calloc(3, (size_t) sizeof(char));
       country[i] = (char *) calloc(3, (size_t) sizeof(char));
    }
    extent_list = (char **) calloc((size_t)MAXEXT, (size_t) sizeof(char *));

/*  set these default values */

    dir = DIR;
    extent_list[0] = EXTENT;
    n_extents = 0;
    mask = 0;
    outfile = stdout;

/*  parse the command line arguments */

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-') {
            switch (argv[i][1]) {
               case 'D':
                        dir = &argv[i][2];
                        break;
               case 'E':
                        extent_list[n_extents] = &argv[i][2];
			++n_extents;
                        break;
               case 'H':
                        h_flag = 1;
                        break;

               case 'h':
                        print_usage(argv[0]);
			exit(0);
                        break;

               case 'O':
                       outfile = create_hydro_file(&argv[i][2], OVERWRITE);
                        if (outfile == NULL) {
                          fprintf(stderr,"\nUnable to open %s for writing.", &argv[i][2]);
                          exit(1);
                        }
                        break;
               case 'R':
                        r_flag = 1;
                        rejectfile = create_hydro_file(&argv[i][2], NOCLOBBER);
                        if (rejectfile == NULL) {
                          fprintf(stderr,"\nUnable to open %s for writing.", &argv[i][2]);
                          fprintf(stderr,"\nIt may exist already?\n");
                          exit(1);
                        }
                        break;

               case 'T':
                        topt = 1;
                        s = &argv[i][2];
                        do {
                          switch (*s) {
                            case 'b' :
                                     b_flag = 1;
                                     ++s;  /* move past b */
                                     if (*s == '^') {
                                          b_flag = -1;
                                          ++s;
                                     }
                                     if (*s == '/')
                                          ++s;

                                     error += (sscanf(s,"%d", &minpdr) != 1);
                                     s = strchr(s,'/');  /* check for another delimiter*/
                                     error = (s == NULL);
                                     if (s != NULL) {
                                        ++s;  /* move past delimiter */
                                        error += (sscanf(s,"%d", &maxpdr) != 1);
                                     }
                                     if (error) break;
                                     if (maxpdr < minpdr) {
                                        fprintf(stderr,"\nmin_seafloor_depth exceeds maxdepth\n");
                                        exit(1);
                                     }
                                     if ((s = strchr(s,'/')) == NULL)  /* check for another delimiter */
                                          break;

                                     if (*s == '/')
                                          ++s;
                                     break;
                            case 'c' :
                                     c_flag = 1;
                                     ++s;
                                     if (*s == '^') {
                                          c_flag = -1;
                                          ++s;
                                     }
                                     if (*s == '/')
                                          ++s;
                                     
                                     do {
                                       error = sscanf(s,"%d", &clist[ncru++]) != 1;
                                       if (error)   
                                           break;
                                       if ((s = strchr(s,'/')) == NULL) 
                                          break;
                                       if (*s == '/')
                                         ++s;
                                     } while ((ncru < MAXCRIT) && (isdigit(*s)));
                                     error = (ncountry >= MAXCRIT);
                                     break;
                            case 'd' :
                                     d_flag = 1;
                                     ++s;  /* move past d */
                                     if (*s == '^') {
                                          d_flag = -1;
                                          ++s;
                                     }
                                     if (*s == '/')
                                          ++s;

                                     error += (sscanf(s,"%d", &mindepth) != 1);
                                     s = strchr(s,'/');  /* check for another delimiter*/
                                     error = (s == NULL);
                                     if (s != NULL) {
                                        ++s;  /* move past delimiter */
                                        error += (sscanf(s,"%d", &maxdepth) != 1);
                                     }
                                     if (error) break;
                                     if (maxdepth < mindepth) {
                                        fprintf(stderr,"\nmindepth exceeds maxdepth\n");
                                        exit(1);
                                     }
                                     if ((s = strchr(s,'/')) == NULL)  /* check for another delimiter */
                                          break;

                                     if (*s == '/')
                                          ++s;
                                     break;
                            case 'e' :
                                     e_flag = 1;
                                     ++s;
                                     if (*s == '^') {
                                          e_flag = -1;
                                          ++s;
                                     }
                                     if (*s == '/')
                                          ++s;
                                     
                                     do {
                                       error = sscanf(s,"%c", &instr[ninstr++]) != 1;
                                       if (error)   
                                           break;
                                       ++s;
				       if (*s == '/')
				           break;
				       if (*s == '\0')
				           break;
                                     } while ((ninstr < MAXCRIT) );
                                     error = (ninstr >= MAXCRIT);
                                     break;
                            case 'g':
                                      g_flag = 1;
                                     ++s;  /* move past p */
                                     if (*s == '^') {
                                           g_flag = -1;
                                          ++s;
                                     }

                                     if (*s == '/')
                                          ++s;
                                     error = (sscanf(s,"%f", &left) != 1);
                                     while (*(++s) != '/')
                                            ;  
                                     ++s;  /* move past delimiter */
                                     error += (sscanf(s,"%f", &right) != 1);
                                     while (*(++s) != '/')
                                            ;  
                                     ++s;
                                     error += (sscanf(s,"%f", &bot) != 1);
                                     while (*(++s) != '/')
                                            ;  
                                     ++s;  /* move past delimiter */
                                     error += (sscanf(s,"%f", &top) != 1);
                                     
                                     if (left > 0 && right < 0)
                                         right += 360.;
                                     if (right > 180)
                                         xdateline = 1;

                                     if ((s = strchr(s,'/')) == NULL) 
                                          break;

                                     if (*s == '/')
                                          ++s;
                                     break;

                            case 'i' :
                                     i_flag = 1;
                                     ++s;  /* move past i */
                                     if (*s == '^') {
                                          i_flag = -1;
                                          ++s;
                                     }

                                     if (*s == '/')
                                          ++s;

                                     error += (sscanf(s,"%d", &minsta) != 1);
                                     s = strchr(s,'/');  /* check for another delimiter*/
                                     error = (s == NULL);
                                     if (s != NULL) {
                                        ++s;  /* move past delimiter */
                                        error += (sscanf(s,"%d", &maxsta) != 1);
                                     }
                                     if (error) break;
                                     if (maxsta < minsta) {
                                        fprintf(stderr,"\nminsta exceeds maxsta\n");
                                        exit(1);
                                     }
                                     if ((s = strchr(s,'/')) == NULL)  /* check for another delimiter */
                                          break;

                                     if (*s == '/')
                                          ++s;
                                     break;
                            case 'l' :
                                     ++s;  /* move past l */
                                     if (*s == 'd')
                                         ld_flag = 1;
                                     else if (*s == 'p')
                                         lp_flag = 1;
                                     else {
                                         fprintf(stderr, "\n Specify p or d after the 'l' extraction type.\n");
                                         exit(1);
                                     }
                                     ++s; 
                                     if (*s == '^') {
                                          if (ld_flag)
                                            ld_flag = -1;
                                          else
                                            lp_flag = -1;
                                          ++s;
                                     }
                                     if (*s == '/')
                                          ++s;
                                     error += (sscanf(s,"%f", &minlev) != 1);
                                     s = strchr(s,'/');  /* check for another delimiter*/
                                     error = (s == NULL);
                                     if (s != NULL) {
                                        ++s;  /* move past delimiter */
                                        error += (sscanf(s,"%f", &maxlev) != 1);
                                     }
                                     if (error) break;
                                     if (maxlev < minlev) {
                                        fprintf(stderr,"\nminlev exceeds maxlev\n");
                                        exit(1);
                                     }
                                     if ((s = strchr(s,'/')) == NULL)  /* check for another delimiter */
                                          break;

                                     if (*s == '/')
                                          ++s;

                                     break;
                            case 'n' :
                                     n_flag = 1;
                                     ++s;  /* move past n */
                                     if (*s == '^') {
                                          n_flag = -1;
                                          ++s;
                                     }
                                     if (*s == '/')
                                          ++s;
                                     do {
                                        error = (sscanf(s,"%2s", country[ncountry++]) != 1);
                                        if (error)   
                                           break;
                                        if ((s = strchr(s,'/')) == NULL) 
                                           break;
                                        if (*s == '/')
                                          ++s;
                                     }  while ((ncountry < MAXCRIT) && (isdigit(*s)));
                                     error += (ncountry >= MAXCRIT);
                                     break;
                            case 'p' :
                                     p_flag = 1;
                                     ++s;  /* move past p */
				     
				     while (*s == '^' || *s == '+' ) {
				     
                                        if (*s == '^') 
                                          p_flag = -1;
				     
                                        if (*s == '+') {
                                          and_flag = 1;
					  or_flag = 0;
					}
					  
                                       ++s;
				     } /* end while */
				     
                                     if (*s == '/')
                                          ++s;
				     nprops = 0;
				     do {
				        cprop[0] = cprop[1]=cprop[2]=cprop[3]= '\0';
					
				        sscanf(s,"%[^'/']",cprop);
				        index = get_prop_indx(cprop);
					if (index >= 0) {
					   proplist[nprops++] = index;
					}
					s += strlen(cprop);   
                                        if (*s == '/')
                                          ++s;
                                     } while (index >= 0 && !( *s == '\0'));
				     
                                     if (nprops <= 0)
                                        error = 1;
					
				     break;
                            case 's' :
                                     s_flag = 1;
                                     ++s;  /* move past s */
                                     if (*s == '^') {
                                          s_flag = -1;
                                          ++s;
                                     }

                                     if (*s == '/')
                                          ++s;
                                     do {
                                        error = (sscanf(s,"%2s", ship[nship++]) != 1);
                                        if (error)   
                                           break;
                                        if ((s = strchr(s,'/')) == NULL) 
                                           break;
                                        if (*s == '/')
                                          ++s;
                                     } while ((nship < MAXCRIT) && (*s < 'a'));
                                     error = (nship >= MAXCRIT);
                                     break;
                            case 'm' :
                                     m_flag = 1;
                                     ++s;  /* move past m */
                                     if (*s == '^') {
                                          m_flag = -1;
                                          ++s;
                                     }

                                     if (*s == '/')
                                          ++s;
                                     do {
                                       error = sscanf(s,"%d", &month) != 1;
                                       if (error)   
                                           break;
                                       m = 1;
                                       mask = (m << (month-1)) | mask;  /* bitwise or */
                                       if ((s = strchr(s,'/')) == NULL) 
                                          break;
                                       if (*s == '/')
                                         ++s;
                                     } while (isdigit(s[0]));
                                     break;
                            case 'y' :
                                     y_flag = 1;
                                     ++s;  /* move past y */
                                     if (*s == '^') {
                                          y_flag = -1;
                                          ++s;
                                     }

                                     if (*s == '/')
                                          ++s;

                                     error += (sscanf(s,"%d", &minyr) != 1);
                                     s = strchr(s,'/');  /* check for another delimiter*/
                                     error = (s == NULL);
                                     if (s != NULL) {
                                        ++s;  /* move past delimiter */
                                        error += (sscanf(s,"%d", &maxyr) != 1);
                                     }
                                     if (error) break;
                                     
                                     if (maxyr < minyr) {
                                        fprintf(stderr,"\nminyr exceeds maxyr\n");
                                        exit(1);
                                     }
                                     if ((s = strchr(s,'/')) == NULL)  /* check for another delimiter */
                                          break;

                                     if (*s == '/')
                                          ++s;
                                     break;
                             default :
                                     error = 1;

                          }  /* end switch */
  
                          if (error)  {
                              fprintf(stderr,"\nError parsing command line");
                              fprintf(stderr,"\n in particular: %s\n", argv[i]);
                              exit(1);          
                          }  
                         
                           
                        } while ((s != NULL) && (*s != '\0'));
                        break;
	        
               default  :
                        fprintf(stderr,"\nError parsing command line");
                        fprintf(stderr,"\n in particular: %s\n", argv[i]);
                        exit(1);
            }  /* end switch */
       }  /* end if */
       else  {
           ++nfiles;
       }
   }  /* end for */

   if ( !nfiles) {
       print_usage(argv[0]);
       fprintf(stderr,"\n\nYou must specify input file(s)!\n");
       exit(1);
   }
   
   if (! topt)  {
     y_flag = 1;
     minyr = 0;
     maxyr = 9999;
   }
   

/* initialize some variables */
   for (i = 0; i < MAXPROP; ++i) {
     data.observ[i] = (double *) NULL;
     data.count[i] = (double *) NULL;
     data.variance[i] = (double *) NULL;
     data.quality[i] = (double *) NULL;
   }
   h.prop_id = (int *)NULL;

   if (n_extents == 0) n_extents = 1;


 /* loop for each input file */

   do {
     curext = 0;
     do {
     infile = open_hydro_file(dir, argv[curfile], extent_list[curext], PRINT_MSG);
     if (infile == NULL) 
       goto NEXTFILE;
     

     /* loop for each station */

     while ((status = get_station(infile, &h, &data)) == 0) { 

       if (xdateline && h.lon < 0)
          h.lon += 360.;
          
       staOK = 1;
       if ( g_flag && staOK){
               staOK = checkpos(&h, top, bot, left, right);
               if ( g_flag < 0)
                   staOK = !staOK;
        }
       if (m_flag && staOK) {
               staOK = checkmonth(&h, mask);
               if (m_flag < 0)
                   staOK = !staOK;
        }

       if (i_flag && staOK) {
               staOK = checksta(&h, minsta, maxsta);
               if (i_flag < 0)
                   staOK = !staOK;
        }
       if (y_flag && staOK) {
               staOK = checkyear(&h,minyr,maxyr);
               if (y_flag < 0)
                   staOK = !staOK;
        }
       if (c_flag && staOK) {
               staOK = checkcruise(&h, clist, ncru);
               if (c_flag < 0)
                   staOK = !staOK;
        }
       if (d_flag && staOK) {
               staOK = checkdepth(&h, data.observ[(int)DE], mindepth, maxdepth);
               if (d_flag < 0)
                   staOK = !staOK;
        }
        if (e_flag && staOK) {
               staOK = checkinstrument(&h, instr, ninstr);
               if (e_flag < 0)
                   staOK = !staOK;
        }
       if (b_flag && staOK) {
               staOK = checkpdr(&h, minpdr, maxpdr);
               if (b_flag < 0)
                   staOK = !staOK;
        }
       if (s_flag && staOK) {
               staOK = checkship(&h, ship, nship);
               if (s_flag < 0)
                   staOK = !staOK;
        }
       if (n_flag && staOK) {
               staOK = checknation(&h, country, ncountry);
               if (n_flag < 0)
                   staOK = !staOK;
        }

       if (p_flag && staOK) {
               staOK = checkprop(&h, proplist, nprops, and_flag, or_flag);
               if (p_flag < 0)
                   staOK = !staOK;
       }
       if ((ld_flag || lp_flag) && staOK) {
          staOK = checklevs(&h, &data, ld_flag, lp_flag, minlev, maxlev);
       }

       if (staOK) {
          if (h_flag) {
             error = write_hydro_hdr(outfile, &h);
          }
          else {
             error = write_hydro_station(outfile, &h, &data);
          }  
          if (error) {
             fprintf(stderr,"\nError in write_hydro_station(): error code = %d.\n", error);
             exit(1);
          }
          ++nout;
       }
       else {
          ++nskip;
          if (r_flag) {
              if (h_flag) {
                error = write_hydro_hdr(rejectfile, &h);
              }
              else {
                error = write_hydro_station(rejectfile, &h, &data); 
              }
          }
       }

     }  /* end while */

     report_status(status, stderr);
     fclose(infile);

NEXTFILE:
     ;
     
    } while (++curext < n_extents);  
   } while (curfile++ < nfiles);

   fprintf(stderr,"\n\n%d stations extracted   %d stations skipped", 
           nout, nskip);
   fprintf(stderr,"\n\nEnd of hb_extract.\n\n");
   exit(0);

} /* end main() */




/****************************************************************************/

void print_usage(char *program)
{
   fprintf(stderr,"\n%s stations from HydroBase station files by", program);
   fprintf(stderr,"\na variety of criteria.\n");

   fprintf(stderr,"\nUsage: %s filename_root(s)", program);
   fprintf(stderr," [-Ddirname] [-Eextent] [-H] [-Ooutcast_file] ");
   fprintf(stderr,"-T[b[^]/bathmin/bathmax]/[c[^]/cruise1/cruise2/...]/[d[^]/depthmin/depthmax]");
   fprintf(stderr,"/[l<d_or_p>[^]/minlev/maxlev]");
   fprintf(stderr,"/[m[^]/m1/m2/...]/[n[^]/code1/code2/...]/");
   fprintf(stderr,"[g[^]/minlon/maxlon/minlat/maxlat]/");
   fprintf(stderr,"[i[^]/minsta/maxsta]/");
   fprintf(stderr,"[p[^][+]/p1[/p2/p3]]/");
   fprintf(stderr,"[s[^]/ship1/ship2/...]/[y[^]/minyear/maxyear]\n");
   fprintf(stderr,"\n    -D  : specifies dirname (default is ./) ");
   fprintf(stderr,"\n        ex: -D../data/ ");
   fprintf(stderr,"\n    -E  : specifies file extent (default is no extent)");  
   fprintf(stderr,"\n        ex: -E.dat ");
   fprintf(stderr,"\n    -H  : output header only");  
   fprintf(stderr,"\n    -O  : specifies output file for stations that meet criteria");  
   fprintf(stderr,"\n        ex: -Onpac_post2003.sta");
   fprintf(stderr,"\n    -R  : specifies file to which rejected stations are written");  
   fprintf(stderr,"\n        ex: -Rnpac_pre2003.sta");
   fprintf(stderr,"\n    -T  : specifies type of extraction and criteria");
   fprintf(stderr,"\n        b = if bathymetric depth falls between min/max ex:  -Tb/4000/10000");
   fprintf(stderr,"\n        c = NODC cruise #  ex:  -Tc/8103/7202/83");
   fprintf(stderr,"\n        d = if deepest level sampled falls between min/max ex:  -Td/4000/10000");
   fprintf(stderr,"\n        e = by equipment (instrument code). Do not separate codes with a slash (/)");
   fprintf(stderr,"\n            ex:  -Tebc  will extract bottle and ctd profiles");
   fprintf(stderr,"\n            ex:  -Te^b  extracts all profiles except bottle data");
   fprintf(stderr,"\n        g = geographic position  ex:  -Tg/-80/-10/0/60");
   fprintf(stderr,"\n        i = by individual station    ex:  -Ti/1/71");
   fprintf(stderr,"\n        ld = depth level  ex:  -Tld1000/2000");
   fprintf(stderr,"\n        lp = pressure level  ex:   -Tlp^1000/2000");
   fprintf(stderr,"\n        m = by month  ex:  -Tm/1/2/3/12");
   fprintf(stderr,"\n        n = nation (country code)ex:  -Tn31/90");
   fprintf(stderr,"\n        p = properties. Use '+' after p to extract station only");
   fprintf(stderr,"\n            if it contains ALL properties in list (logical AND)");
   fprintf(stderr,"\n            otherwise it will be extracted if it contains ANY properties in list (logical OR)");
   fprintf(stderr,"\n            ex:  -Tp+sa/si  extract if salt AND silicate present");
   fprintf(stderr,"\n            ex:  -Tp^ox/o2  extract if  NOT (ox OR o2)");
   fprintf(stderr,"\n        s = ship  ex:  -TsVE/6N");
   fprintf(stderr,"\n        y = by year   ex:  -Ty/1965/1980");
   fprintf(stderr,"\n        combination of above:  -Tg/-80/-10/0/160 -Tm/1/2/3/12");
   fprintf(stderr,"\n\n  Use ^ after the extraction type to extract stations");
   fprintf(stderr,"\n   which DO NOT meet the criteria. ");
   fprintf(stderr,"\n   ex: -Ts^/VE/y1982  to extract all stations from 1982");
   fprintf(stderr,"\n       except those collected by the ship VE.");
   fprintf(stderr,"\n -h  help ... prints this message.");
   fprintf(stderr,"\n\n");  
   return;
}

/****************************************************************************/

int checkmonth(struct HYDRO_HDR *hptr, int mask)

/*  Returns a nonzero value if the station was done during one of the 
    specified months  or 0 if it does not */
{
   int  m = 1;

   return ((m << (hptr->month -1)) & mask);    /* bitwise and */

} /* end checkmonth() */

/****************************************************************************/
int checkprop(struct HYDRO_HDR *hptr, int *plist, int n, int aflag, int oflag)
{
   int i, navail;
   
   navail = 0;
   i = -1;
   while (++i < n) {
      if (available((enum property) plist[i], hptr) ) {
         if (oflag)
            return(1);
	 ++navail;
      }
   }
   if (aflag && (navail == n)) {
     return(1);
   }
   if (!aflag) {
    if (navail > 0)
      return(1);
   }  
   
   return(0);

} /* end checkprop() */


/****************************************************************************/


int checkpos(struct HYDRO_HDR *hptr, float top, float bot, float left, float right)
{

   if ((hptr->lat > top) || (hptr->lat < bot))
        return (0);

   if ((hptr->lon > right) || (hptr->lon < left))
        return (0);

   return (1);

} /* end checkpos() */

/****************************************************************************/
int checkpdr(struct HYDRO_HDR *hptr, int min, int max)
{

   return ((hptr->pdr <= max) && (hptr->pdr >= min));

} /* end checkpdr() */

/****************************************************************************/
int checkyear(struct HYDRO_HDR *hptr, int min, int max)
{

   return ((hptr->year <= max) && (hptr->year >= min));

} /* end checkyear() */

/****************************************************************************/
int checksta(struct HYDRO_HDR *hptr, int min, int max)
{

   return ((hptr->station <= max) && (hptr->station  >= min));

} /* end checksta() */


/****************************************************************************/
int checknation(struct HYDRO_HDR *hptr, char **list, int n)
  /*   list:  list of country codes
       n:      number of country codes in list 
  */
{
   int i;

   i = -1;
   while (++i < n) {
        if (strncmp(hptr->country,list[i],2) == 0)
            return(1);
   }
   return(0);

} /* end checknation() */
/****************************************************************************/
int checkship(struct HYDRO_HDR *hptr, char **list, int n)
  /*  list:    list of ship codes 
         n:    number of ship codes in list 
  */
{
   int i;

   i = -1;
   while (++i < n){
        if (strncmp(hptr->ship,list[i],2) == 0)
            return(1);
   }
   return(0);

} /* end checkship() */
/****************************************************************************/
int checkinstrument(struct HYDRO_HDR *hptr, char *list, int n)
  /*  list:    list of instrument codes 
         n:    number of  codes in list 
  */
{
   int i;

   i = -1;
   while (++i < n){
        if (hptr->instrument == list[i])
            return(1);
   }
   return(0);

} /* end checkship() */
/****************************************************************************/
int checkdepth(struct HYDRO_HDR *hptr, double *z, int min, int max)

   /* Returns a 1 if the bottommost observation is at a level between
      min and max (inclusive); or a 0 if not. */
{

   return ((z[hptr->nobs-1] <= max) && (z[hptr->nobs-1] >= min));

} /* end checkdepth() */

/****************************************************************************/
int checkcruise(struct HYDRO_HDR *hptr, int *list, int n)
    /*  list:      list of NODC cruise #s 
        n:         number of items in list 
    */
{
     int i;

   i = -1;
   while (++i < n) {
      if (list[i] == hptr->cruise)
            return(1);
   }
   return(0);

} /* end checkcruise() */

/****************************************************************************/
int checklevs(struct HYDRO_HDR *hptr, struct HYDRO_DATA *dptr, int d_flag, int p_flag, float min, float max)

  /* d_flag, p_flag:  -1, 0, or 1 for depth/pressure levels */

  /*  Determines which pressure or depth levels meet the min/max criteria
       and orders the data accordingly, adjusting nobs if necessary.
       Returns a 1 if any levels meet the criteria or 0 if no levels do. */

{
   int i, j, k, ii, not_flag, levOK;
   double *prop;

   not_flag = 0;

   if (d_flag) {
      if (d_flag < 0)  
         not_flag = 1;
      prop = dptr->observ[(int)DE];
   }
   else {
      prop = dptr->observ[(int)PR];
      if (p_flag < 0)
         not_flag = 1;
   }

   k = 0;
   for (j = 0; j < hptr->nobs; ++j) {
      levOK = (prop[j] <= max) && (prop[j] >= min);
      if (not_flag)
           levOK = !levOK;
      if (levOK) {
        for (i = 0; i < hptr->nprops; ++i) {
            ii = hptr->prop_id[i];
            dptr->observ[ii][k] = dptr->observ[ii][j];
        }
        ++k;
      }
   }
   hptr->nobs = dptr->nobs = k;

   if (k == 0) 
       return(0);

   return(1);

} /* end checklevs() */
