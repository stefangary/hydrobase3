/* hb_smoothsurf2d.c
 ................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth G Curry
                             Woods Hole Oceanographic Instition
                             2000  ANSI compliant
			     updated Dec 2010 for HB3
...............................................................................
--------------------------------------------------------------------
.
 * 
 * Smooths a 2D gridded field using a gaussian weighted average
 * and a specified ellipse -S<xradius>/<yradius>.
         weight[n] = e^[- alpha * dist^2 ]
	 alpha = (pi/radius)^2 /4.5 (or a value specified with -A)
	 dist = sqrt(dx*dx + dy*dy) where dx, dy are earth distances computed 
	         from lat/lon coordinates. 
 * If a masking file is supplied, it will be treated as topographic info --
 * smoothing will not cross topographic ridges.
 *  If -F<xdist>/<ydist> is specified, it will attempt to fill empty (non-masked) nodes 
 *  using a nearest-neighbor algorithm and the x, y length scales (cartesian
 *  distances) specified.

 *-------------------------------------------------------------------- */

#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "hydrobase.h"
#include "hb_grids.h"
#include "hb_memory.h"

/* globally defined variables */

double xlen, ylen;        /* length scale limits for interpolating */
double xradius, yradius;  /* for smoothing ellipse */

/* Define internal representation of masked / flagged values    
   The value for empty/masked nodes written out is defined 
   in hb_grids.h as TOOBIG */
   
double empty_val = -TOOBIG;         /* large negative value */
double mask_val  = TOOBIG;          /* large positive value */
double flagged = -TOOBIG +1;        /* check for input missing values  */
double masked = TOOBIG -1;        /* check for input masked values  */

/*  prototypes for locally declared functions */

void print_usage(char *);

main (int argc, char **argv)
{
  int bflag, iflag, sflag, fill_flag;
  int gflag, fgflag;
  int i, jrow, icol,  n, nread, n_fields;
  int ix, iy, n_set, n_mask, n_empty;
  int nfiles = 0, curfile = 1;
  int nsq, row, col, sq;
  int rad, found;
  int *zcount;
  int nwsq, wsq,wrow, wcol;
  int xrad, yrad;
 

  BOOLEAN error;
  BOOLEAN  yfirst;
  BOOLEAN  skip;
  BOOLEAN use_mask;

  double in[3], x_left, x_right, y_top, y_bottom; 
  double xoffset, yoffset;
  double lat, lon, dydist, dxdist;
  double x1, x2, y1, y2, z1, z2;
  double x, y, zout, wx, wy;
  double sd, aw;
  double *z, *zwork; 
  double *weights, *wsave;   
  double alpha;


  char line[BUFSIZ];
  char *st, **sptr;
  char *mask;
  
  FILE *outfile;
  FILE *infile;
  FILE *maskfile;
  
  struct GRID_INFO h, hwork;
  
  /* initialize these values */

  use_mask  = FALSE;
  error = 0;
  yfirst = 0;
  bflag = iflag = sflag  = fill_flag = 0;
  gflag = fgflag = 0;
  n_mask = 0;
  xlen  = ylen = 100.0;  /* length scales for filling empty nodes */
  outfile = stdout;
  infile = stdin;
  h.node_offset = 1;  /* default pixel registration */
  mask = (char *)NULL;  
  xradius = yradius = 100.0;
  alpha = 0.0;
  
  
  if (argc < 1) {
    print_usage(argv[0]);
    exit(1);
  }
  
  /* parse command line arguments */

  
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
        case 'A':
         error = sscanf(&argv[i][2],"%lf", &alpha) != 1;
          break;
      
        case 'B':  /* get grid bounds */
	  bflag = 1;
          st = &argv[i][2];
          if (*st == '/')
             ++st;
          error = (sscanf(st,"%lf", &h.x_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &h.x_max) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &h.y_min) != 1);
          while (*(st++) != '/')
                           ;  
          error += (sscanf(st,"%lf", &h.y_max) != 1);
                        	     
	  if (h.x_min > h.x_max) {
	    fprintf(stderr,"\nxmin bound must be numerically < xmax bound");
	    exit(1);
	  }
	  
	  if (h.y_min > h.y_max) {
	    fprintf(stderr,"\nymin bound must be numerically < ymax bound");
	    exit(1);
	  }
          break;
	  
        case 'F':
	  fill_flag = 1;
	  st = &argv[i][2];
	  if (*st == 'g') {
	      fgflag = 1;
	      ++st;
	  }
	  
          error = sscanf(st,"%lf", &xlen) != 1;
          st = strchr(st,'/');
	  ylen = xlen;
          if (st != NULL) {
             sscanf(++st,"%lf", &ylen);
          }
          break;

	  
        case 'G':
	  h.node_offset = 0;
	  break;
        case 'I':
	  iflag = 1;
          error = (sscanf(&argv[i][2], "%lf", &h.x_inc) == 1) ? 0 : 1;
	  h.y_inc = h.x_inc;
          st = strchr(&argv[i][2],'/');
          if (st != NULL) {
             sscanf(++st,"%lf", &h.y_inc);
          }
          break;
	  
        case 'M':
          maskfile = fopen(&argv[i][2],"r");
          if (maskfile == NULL) {
                fprintf(stderr,"\nError opening %s for reading.\n",&argv[i][2]);
                exit(1);
          }
          use_mask = TRUE;
          break;
	  
        case 'O':
          outfile = fopen(&argv[i][2],"w");
          if (outfile == NULL) {
                fprintf(stderr,"\nError opening %s for output.\n",&argv[i][2]);
                exit(1);
          }
          break;

        case 'P':
	  h.node_offset = 1;
	  break;
	  
        case 'S':
	  sflag = 1;
	  st = &argv[i][2];
	  if (*st == 'g') {
	      gflag = 1;
	      ++st;
	  }
          error = sscanf(st,"%lf", &xradius) != 1;
          st = strchr(st,'/');
	  yradius = xradius;
          if (st != NULL) {
             sscanf(++st,"%lf", &yradius);
          }
          break;

	case 'h':  
	   print_usage(argv[0]);
	   exit(0);
	   
	case ':':
	   yfirst = 1;
	   break;

        default:
          error = TRUE;
          break;
        } /* end switch */
	  
       if (error ) {
         fprintf(stderr,"\nError parsing command line args.\n");
         fprintf(stderr,"     in particular: '%s'\n", argv[i]);
         exit(1);
       }

     }
     else 
        ++nfiles;
	
    }  /* end for  */
    
    
    if (!bflag || !iflag ) {
       fprintf(stderr,"\nYou must specify -B<bounds> and -I<gridspacing>");
       fprintf(stderr,"\nUse -h for complete usage info. \n");
      exit(1);
    
    }
    
   error = 0;
   
  if (xlen < 0 || ylen < 0) {
    fprintf (stderr, "SYNTAX ERROR -F option.  Search radius must be positive distance in km or in # of gridpoints.\n");
    error++;
    }
  if (error) 
        exit(1);
 
  ix = yfirst; 
  iy = 1 - ix;
        
  
   xoffset = 0.5*h.x_inc;
   yoffset = 0.5*h.y_inc;
     
  /*  for default pixel registration */
  
   h.nx = (int) NINT((h.x_max - h.x_min) / h.x_inc);
   h.ny = (int) NINT((h.y_max - h.y_min) / h.y_inc);
   x_left = h.x_min;    /* data boundaries */ 
   x_right = h.x_min  + h.nx * h.x_inc;
   y_bottom = h.y_min;
   y_top = h.y_min + h.ny * h.y_inc;  

  
 /* for gridnode registration , adjust bounds to be a pixel grid*/
 
   if ( h.node_offset == 0) {        
     ++h.nx;
     ++h.ny;
     x_left = h.x_min - xoffset;   
     x_right = h.x_min + (h.nx-1) * h.x_inc + xoffset;
     y_bottom = h.y_min - yoffset;
     y_top = h.y_min + (h.ny-1) * h.y_inc + yoffset;  
  }


  fprintf (stderr, "Grid dimensions are nx = %d, ny = %d\n", h.nx, h.ny);
  if (h.node_offset)
     fprintf (stderr, "\n using pixel registration. ");
  else
     fprintf (stderr, "\n using gridnode registration. ");

   
/* allocate space for z arrays */

   nsq = h.nx * h.ny;
   
   z = (double *) get_memory((void *)NULL, nsq, sizeof(double));
   zcount = (int *) get_memory((void *)NULL, nsq, sizeof(int));
   
   nread = 0;  /* counts number of points read in */

 /* loop for each input file :  attach each point read in to its
    nearest grid point.  Sum values and find average in case multiple points fall on a grid node*/
   do {
  
     if (! nfiles)
        fprintf(stderr,"\nExpecting input from stdin ...");
      
     else {
        infile = fopen(argv[curfile],"r");
        if (infile == NULL) {
         fprintf(stderr, "\nUnable to open %s for reading.", argv[curfile]);
         goto NEXTFILE;
        }
        fprintf(stderr,"\nOpened %s ", argv[curfile]);
     }
   

     while (fscanf(infile,"%[^\n]", line) != EOF) { 
        error = getc(infile);   /* move past newline */
	
	if ((st = strchr(line,'N')) != NULL) continue; /* check for NaN */
	
        n_fields = sscanf(line,"%lf %lf %lf", &in[0], &in[1], &in[2]);

        if (n_fields != 3) {
           fprintf(stderr, "Mismatch between actual (%d) and expected (3) fields near line %d\n", n_fields, nread);
           exit(1);
        }
	
        skip = FALSE;
	x = in[ix];
	y = in[iy];
        if (x < x_left || x > x_right) skip = TRUE;
        if (y < y_bottom || y > y_top) skip = TRUE;
        if (ABS(in[2]) > masked) skip = TRUE;   /* weed out empty/mask flags */
    
        if (!skip) {
	   error = xy2ij(&h, x, y, &col, &row);
	   sq = row * h.nx + col;
	   z[sq] += in[2];
	   ++zcount[sq];
           ++nread;
         }  /* End if skip. */
	
    }  /* End while  */
    
    
    if (nfiles) fclose (infile);

      
NEXTFILE:
     ;

  } while (++curfile < nfiles);      /* End input phase */

/*  Find average values for each gridnode that has data attached to it */
      n_empty = 0; 
   for (sq = 0; sq < nsq; ++sq) {
      if (zcount[sq] > 0)
          z[sq] /= zcount[sq];
      else
          z[sq] = empty_val;
   }  
 
  if (use_mask) {
     mask = (char *) get_memory((void *)NULL, (size_t)(nsq), sizeof(char *));
     fprintf(stderr,"\nReading mask file....");
     n_mask = get_mask(maskfile, &h, mask, yfirst);
     for (i = 0; i < nsq; i++)  {
       if (mask[i] )
          z[i] = mask_val;
     }
     free((void *)mask);
  }
  
/* At this point, every element of z is either a value, empty or masked */

/*****************************/
    

   if (fill_flag) {
   
      fprintf (stderr, "Filling empty nodes...\n");
    fprintf (stderr, "length scales:  xlen = %lf  ylen = %lf ", xlen, ylen);
    if ( fgflag)
      fprintf (stderr, " gridpoints\n");
    else
       fprintf (stderr, " km\n");
       
    /* visit each gridpoint in big grid and attempt to fill empty nodes */

      hwork.x_inc = h.x_inc;
      hwork.y_inc = h.y_inc;
      hwork.node_offset = 0;  /* work grids are gridline registered */
      
      yrad = NINT (ylen);  /* Search radius in gridpoints */
      xrad = NINT (xlen);
      
      if (! fgflag) {   /* Search radius was specified in km */
         dydist = h.y_inc * RADperDEG * EarthRadius * KMperNM;
         yrad = (int) floor(ylen / ABS(dydist));
      }
      
      for (sq = 0; sq < nsq; ++sq) {
   
         if (z[sq] < flagged) {
           jrow = (int) floor(sq / h.nx);
	   icol = sq - jrow * h.nx;
	   ij2xy (&h, icol, jrow, &lon, &lat);
	   
           if (! fgflag) {
              dxdist = h.x_inc * RADperDEG * cos(RADperDEG * ABS(lat)) * EarthRadius * KMperNM;
	      xrad = (int) floor(xlen / dxdist);
	   }
	   
           hwork.nx = xrad * 2 + 1;
           hwork.ny = yrad * 2 + 1; 
           nwsq = hwork.nx * hwork.ny;
	   
	   hwork.x_min =  lon - xrad * hwork.x_inc;
	   hwork.x_max =  lon + xrad * hwork.x_inc;
	   hwork.y_min =  lat - yrad * hwork.y_inc;
	   hwork.y_max =  lat + yrad * hwork.y_inc;
	   
           zwork = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
	   weights = (double *)  get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
	   
	   if (fgflag) 
	      get_weights_g(weights, &hwork, alpha, xrad, yrad);
	   
	   else 
	      get_weights_c(weights, &hwork,alpha, xlen, ylen, xrad, yrad); 
	  
	   
           /* fill zwork with elements of z */
	   
	   for (wrow = 0; wrow < hwork.ny; ++wrow) {
	       for (wcol = 0; wcol < hwork.nx; ++wcol) {
 	          wsq = wrow * hwork.nx + wcol;
	          error = ij2xy(&hwork, wcol, wrow, &wx, &wy);
	          error = xy2ij(&h, wx, wy, &ix, &iy);

	          if (error < 0)
		     zwork[wsq] = mask_val;
	          else 
		     zwork[wsq] = z[iy * h.nx + ix];
	       }  /* for wcol */
	   } /* for wrow */
	   	   
	   nn_interp2d(zwork, weights, empty_val, mask_val, xrad, yrad, &hwork, &z[sq], &sd, &aw, &n);
	   
	   free(weights);
	   free(zwork);
     
         } /* end if flagged */
      }  /* end for sq */
    } /* end if fill flag */

/*****************************/
   if (sflag) {
      fprintf(stderr,"\nNow smoothing over ellipse with radii %lf/%lf ", xradius, yradius); 
      if (gflag)
         fprintf(stderr," pts\n");
      else
         fprintf(stderr," km\n");
      
      n_empty = 0; 

      xrad = NINT (xradius);  
      yrad = NINT (yradius);

      
      if (!gflag) {
         dydist = h.y_inc * RADperDEG * EarthRadius * KMperNM;
         yrad = NINT (floor(yradius / dydist));
      }
      
      hwork.x_inc = h.x_inc;
      hwork.y_inc = h.y_inc;
      hwork.node_offset = 0;  /* work space is always grid registered */
      
      
  /*visit each square in big grid, smooth nodes with actual values and output x,y,z triplets */
  
      for (sq = 0; sq < nsq; ++sq ) {
      
         if (z[sq] < flagged) {
	   ++n_empty;
         }
      
         else if (z[sq] < masked ) {   /* this node not empty or masked */
            jrow = NINT (floor(sq / h.nx));
            icol = sq - jrow * h.nx;
            error = ij2xy (&h, icol, jrow, &lon, &lat);
      
            if (!gflag) {
               dxdist = h.x_inc * RADperDEG * cos(RADperDEG * lat)* EarthRadius * KMperNM;
               xrad = NINT (floor(xradius / dxdist));
            }
      
     /* set up work grid and weights*/
     
	     hwork.x_min =  lon - xrad * hwork.x_inc;
	     hwork.x_max =  lon + xrad * hwork.x_inc;
	     hwork.y_min =  lat - yrad * hwork.y_inc;
	     hwork.y_max =  lat + yrad * hwork.y_inc;
             hwork.nx = xrad * 2 + 1;
             hwork.ny = yrad * 2 + 1; 
             nwsq = hwork.nx * hwork.ny;

             weights = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
             zwork = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
      
            if (gflag)
               get_weights_g(weights, &hwork, alpha, xrad, yrad); 
            else 
               get_weights_c(weights, &hwork, alpha, xradius, yradius, xrad, yrad); 
   
 	    for (wrow = 0; wrow < hwork.ny; ++wrow) {
	       for (wcol = 0; wcol < hwork.nx; ++wcol) {
	         
		  wsq = wrow * hwork.nx + wcol;
	          error = ij2xy(&hwork, wcol, wrow, &wx, &wy);
	          error = xy2ij(&h, wx, wy, &ix, &iy);
		  
	          if (error < 0)
		     zwork[wsq] = mask_val;
	          else 
		     zwork[wsq] = z[iy * h.nx + ix];
	        }  /* for wcol */
	     } /* for wrow */
	
              wsq = yrad * hwork.nx + xrad;
	      if (ABS(zwork[wsq])< masked){
	         zero_weights_2d(zwork, weights, empty_val, mask_val, &hwork, xrad, yrad);
	         zout = weighted_mean(zwork, weights, empty_val, mask_val, nwsq);
	   
	         if (ABS(zout) < masked)
	            fprintf(outfile,"%lf %lf %lf\n", lon, lat, zout);
	      }
	   
	      free(weights);
	      free(zwork);
         } /* end else if */  
       }  /* end for sq*/
    } /* end if sflag */

  fprintf(stderr,"\nNumber read in: <%d>  filled: <%d>   empty: <%d>   masked: <%d>\n", nread, n_set, n_empty, n_mask );
  fprintf(stderr,"\nEnd of %s.\n", argv[0]);
  
  exit(0);
}

/******************************************************************/
void print_usage(char *program)
{
  
    fprintf (stderr, "%s -  smoothes the 2D field using distance-weighted gaussian averaging with variable x/y radii \n", program);
    fprintf(stderr, "Outputs lon,lat,zsmooth triplets \n\n");
    fprintf(stderr, "USAGE: %s [xyzfile(s)] [-A<alpha>] -B<xmin/xmax/ymin/ymax> -I<dx>[/<dy>]", program);
    fprintf(stderr, "[-F[g]<xlen>[/<ylen>]] [-G] [-M<mask_file> ] [-O<output_file>]");
    fprintf(stderr, "  [-S[g]<xradius>[/<yradius>]]  [-:] [-h]\n\n");
    fprintf(stderr, "-B   sets the grid bounds in user units: xmin/xmax/ymin/ymax.\n");
    fprintf(stderr, "-I   sets the grid spacing for the x/y directions.\n");
    fprintf(stderr, "\n\tOPTIONS:\n");
    fprintf(stderr, "-A specify weighting factor alpha: e^[-alpha * dist^2] .\n");
    fprintf(stderr, "      default:  alpha = (pi/radius)^2 /4.5   \n");
    fprintf(stderr, "-F fill empty (non-masked) nodes using x,y length scales (in km) \n");
    fprintf(stderr, "      If no data are within this range of a node, it  remains empty.\n");
    fprintf(stderr, "      Default is [-F%lf/%lf]\n", xlen, ylen);
    fprintf(stderr, "   To specify length scales in grid points instead of km, append 'g' after -F\n");       fprintf(stderr, "     ex: -Fg<xlength>[/<ylength>]\n");
    fprintf(stderr, "-G  force gridnode registration (default is pixel registration).\n");
    fprintf(stderr, "-M name of file containing masking info.\n");
    fprintf(stderr, "      Each line of file may specify individual points to mask\n");
    fprintf(stderr, "      or multiple polygons separated by a '>' character. \n");
    fprintf(stderr, "      Specify 'I' or 'O' immediately after the '>' character\n");
    fprintf(stderr, "      to mask INSIDE or OUTSIDE of polygon.\n");
    fprintf(stderr, "      No '>' or a 'C' immediately after the '>'  \n");
    fprintf(stderr, "      signifies cell mode masking.  In this mode\n");
    fprintf(stderr, "      an {x,y} pair OR {x,y,mask} triplet is given.\n");
    fprintf(stderr, "      The <mask> value can be 1 to mask the point\n");
    fprintf(stderr, "      or 0 (zero) to unmask the point.\n");
    fprintf(stderr, "-O  name of output file  (default is stdout).\n");
    fprintf(stderr, "-S set xradius/yradius of smoothing ellipse (distance in km)\n");
    fprintf(stderr, "      default is: [-S%lf/%lf].\n", xradius, yradius);
     fprintf(stderr, "   To specify radius in grid points instead of km, use -Sg<xradius>[/<yradius>]\n");

    fprintf(stderr, "-: input data are ordered y x z  \n");
    fprintf(stderr, "       [default order is x y z]   \n");
    fprintf(stderr, "-h help....prints this message.  \n");
    return;
}  /* end print_usage() */

