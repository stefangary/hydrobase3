/* hb_gridsection.c
 ................................................................................
                          *******  HydroBase 3 *******
               ...................................................
                          
                    author:  Ruth G Curry
                             Woods Hole Oceanographic Instition
                             2000  ANSI compliant
			     updated to HB3M Jan 2010
...............................................................................
--------------------------------------------------------------------
.
 * Interpolates xyz values onto a grid by linear interpolation in the vertical (y)
 * direction first, followed by interpolation in the horizontal (x) direction according
 * to x and y length scales specified (as integer # of gridnodes). 
 * Optionally smooths the gridded values with a gaussian weighted average of values within
 * a specified ellipse -Sxradius/yradius.
         weight[n] = e^[- (pi/radius)^2 /4.5 * dist^2 ]
	 dist = sqrt(dx*dx + dy*dy) where dx, dy are distance from gridnode.  

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

int xlen, ylen;        /* length scale limits for interpolating */
int xradius, yradius;   /* for smoothing ellipse */

/* Define internal representation of masked / flagged values    
   The value for empty/masked nodes written out is defined 
   in hb_grids.h as TOOBIG */
   
double empty_val = -TOOBIG;         /* large negative value */
double mask_val  = TOOBIG;          /* large positive value */
double flagged = -TOOBIG +1;        /* check for input missing values  */
double masked = TOOBIG -1;        /* check for input masked values  */

/*  prototypes for locally declared functions */

void print_usage(char *);
void get_weights( double *, struct GRID_INFO *, int, int);

main (int argc, char **argv)
{
  int bflag, iflag, sflag;
  int i, j,  nread, n_fields;
  int ix, iy, n_set, n_mask, n_empty;
  int nfiles = 0, curfile = 1;
  int nsq, row, col, sq;
  int rad, found1, found2;
  int *zcount;
  int nwsq, wsq,wrow, wcol;
  int npts;

  BOOLEAN error;
  BOOLEAN  yfirst;
  BOOLEAN  skip;
  BOOLEAN use_mask;

  double in[3], x_left, x_right, y_top, y_bottom; 
  double xoffset, yoffset;
  double x1, x2, y1, y2, z1, z2;
  double x, y, zout, wx, wy;
  double stddev, wave;
  double *z, *zwork, *xwork, *ywork; 
  double *weights, *wsave;   


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
  bflag = iflag = sflag  = 0;
  n_mask = 0;
  xlen  = 0;
  ylen = 0;
  outfile = stdout;
  infile = stdin;
  h.node_offset = 0;  /* default gridnode registration */
  mask = (char *)NULL;  
  xradius = yradius = 0;
  outfile = stdout;
  
  
  if (argc < 1) {
    print_usage(argv[0]);
    exit(1);
  }
  
  /* parse command line arguments */

  
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      
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
	  
        case 'I':
	  iflag = 1;
          error = (sscanf(&argv[i][2], "%lf", &h.x_inc) == 1) ? 0 : 1;
	  h.y_inc = h.x_inc;
          st = strchr(&argv[i][2],'/');
          if (st != NULL) {
             sscanf(++st,"%lf", &h.y_inc);
          }
          break;
	  
        case 'L':
          error = sscanf(&argv[i][2],"%d", &xlen) != 1;
          st = strchr(&argv[i][2],'/');
	  ylen = xlen;
          if (st != NULL) {
             sscanf(++st,"%d", &ylen);
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
          error = sscanf(&argv[i][2],"%d", &xradius) != 1;
          st = strchr(&argv[i][2],'/');
	  yradius = xradius;
          if (st != NULL) {
             sscanf(++st,"%d", &yradius);
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
   
  if (xlen == 0 ) {
    fprintf (stderr, "Must specify search lengthscale(s) (as # of integer gridnodes) with -L.\n");
    exit(1);
    }
  if (xlen < 0 || ylen < 0) {
    fprintf (stderr, "SYNTAX ERROR -L option.  Search radius must be positive.\n");
    exit(1);
    }
 
  ix = yfirst; 
  iy = 1 - ix;
  
  /* ensure x axis is multiple of x_inc 
  
  h.nx = ceil(h.x_max / h.x_inc);
  h.x_max = h.nx * h.x_inc;  */
        
   fprintf (stderr, "length scales:  xlen = %d   ylen = %d\n", xlen, ylen);
   xoffset = 0.5*h.x_inc;
   yoffset = 0.5*h.y_inc;
     
  /*  for default gridnode registration */
  
   h.nx = 1 + (int) NINT((h.x_max - h.x_min) / h.x_inc);
   h.ny = 1 + (int) NINT((h.y_max - h.y_min) / h.y_inc);

      /* compute data bounds */
   
   x_left = h.x_min - xoffset;   
   x_right = h.x_min + (h.nx-1) * h.x_inc + xoffset;
   y_bottom = h.y_min - yoffset;
   y_top = h.y_min + (h.ny-1) * h.y_inc + yoffset;  
  
 /* for pixel registration */
 
   if ( h.node_offset == 1) {        
     --h.nx;
     --h.ny;
     x_left = h.x_min;    /* data boundaries */ 
     x_right = h.x_min  + h.nx * h.x_inc;
     y_bottom = h.y_min;
     y_top = h.y_min + h.ny * h.y_inc;  
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
      
NEXTFILE:
    if (nfiles) fclose (infile);

  } while (++curfile < nfiles);      /* End input phase */

/*  Find average values for each gridnode that has data attached to it */
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

  ywork = (double *) calloc(h.ny, sizeof(double));
  zwork = (double *) calloc(h.ny, sizeof(double));
  n_set = 0;
  n_empty = 0;

/* explicitly generate y-vector corresponding to grid*/  
  for (row = 0; row < h.ny; ++ row) {
     ij2xy(&h, 0, row, &x, &y);
     ywork[row] = (double) y;
  }
/* interpolate vertically column by column... */
  
  for (col = 0; col < h.nx; ++col) {
      for (row = 0; row < h.ny; ++row) {
         sq = row * h.nx + col;
	 zwork[row] = z[sq];
      }
      
      for (row= 0; row < h.ny; ++row) {
          if (zwork[row] <= flagged ) {
	     /* search backward and forward in zwork array to find neighbors within ylen;
	     Suspend search if a masked neighbor is encountered.*/
	     found1 = 0;
	     rad = 1;
	     i = row -1;
	     while (i >= 0 && rad <= ylen && !found1) {

	         if (zwork[i] >= masked) {
		     rad = ylen +1;   /* suspend search */
		 }
		 else if (zwork[i] <= flagged) {
		   --i;
		   ++rad;
		 }
		 else  {
		   found1 = 1;
		 }
	     } /* end while */
	     
	     if (found1) {
	         z1 = zwork[i];
	         y1 = ywork[i];
	     
	        found2 = 0;
	        rad = 1;
	        i = row + 1;
	        while (i < h.ny && rad <= ylen && !found2) {
	            if (zwork[i] >= masked) {
		        rad = ylen +1;   /* suspend search */
		     }
		     else if (zwork[i] <= flagged) {
		         ++i;
			 ++rad;
		     }
		     else  {
		       found2 = 1;
	               z2 = zwork[i];
	               y2 = ywork[i];
		       zwork[row] = z1 + (z2 - z1) * (ywork[row] - y1) / (y2 - y1);
		       ++n_set;
		    }
	        }  /* end while */
	     }  /* end if found2 */
	  } /*end if zwork[row] */
      }/* end for row */
      
      /* load zwork back into z */
      
       for (row = 0; row < h.ny; ++row) {
         sq = row * h.nx + col;
	 z[sq] = zwork[row];
      }
     
  }  /* end for col */

  free(zwork);
  
/* interpolate horizontally row by row ... */
  xwork = (double *) calloc(h.nx, sizeof(double));
  zwork = (double *) calloc(h.nx, sizeof(double));
  
/* explicitly generate x-vector corresponding to grid*/  
  for (col = 0; col < h.nx; ++ col) {
    ij2xy(&h, col, 0, &x, &y);
    xwork[col] = x;
  }
  
  for (row = 0; row < h.ny; ++row) {
    for (col = 0; col < h.nx; ++col) {
      sq = row * h.nx + col;
      zwork[col] = z[sq];
    } 
       
    for (col = 0; col < h.nx; ++col) {
      if (zwork[col] <= flagged) {  
	/* search backward and forward in array to find neighbors within xlen.
	   Interpolate if 2 neighbors are found, Extrapolate if only 1 */
	    
	found1 = 0;
	rad = 1;
	i = col-1;
	while (i >= 0 && rad <= xlen && !found1) {
	  if (zwork[i] >= masked)  {
	    rad = xlen + 1;   /* suspend search */
	  }
	  else if (zwork[i] <= flagged) {
	    --i;
	    ++rad;
	  }
	  else {
	    found1 = 1;
	  }
	} /* end while */
	        
	z1 = zwork[i];
	x1 = xwork[i];
	       
	found2 = 0;
	rad = 1;
	i = col + 1;
	while (i < h.nx && rad <= xlen && !found2 ) {
	  if (zwork[i] >= masked)  {
	    rad = xlen + 1;  /* suspend search */
	  }
	  else if (zwork[i] <= flagged) {
	    ++i;
	    ++rad;
	  }
	  else  {
	    found2 = 1;
	  }
	  
	  
	  z2 = zwork[i];
	  x2 = xwork[i];
	  
	  if (found1 && found2) {
	    zwork[col] = z1 + (z2 - z1) * (xwork[col] - x1) / (x2 - x1);
	    ++n_set;
	  }
	  else if (found1){
	    zwork[col] = z1;
	    ++n_set;
	  }
	  else if (found2) {
	    zwork[col] = z2;
	    ++n_set;
	  } 
	  else {
	    ++n_empty;  
	  }
	  
	}  /* end while */ 
	
      }  /* end if zwork[col] */
    }/* end for col */
       
    /* load zwork back into z */
      
    for (col = 0; col < h.nx; ++col) {
      sq = row * h.nx + col;
      z[sq] = zwork[col];
    }

  } /*end for row */ 
  
  free(zwork);

  fprintf(stderr,"\nNumber of gridnodes: <%d>", nsq);
      

/****************************************************/
  if ( !sflag) {
  /* output xyz triplets */
   
      for (row = 0; row < h.ny; ++row) {
         for (col = 0; col < h.nx; ++col) {
            sq = row * h.nx + col;
            if ( z[sq] > flagged && z[sq] < masked) {
	       fprintf(outfile, "%8.3lf %8.3lf %.9g \n", xwork[col], ywork[row], z[sq]);
	    }

         } /* end for col */
      } /* end for row */
      
      fprintf(stderr,"\nNumber read in: <%d>  filled: <%d>   empty: <%d>   masked: <%d>\n", nread, n_set, n_empty, n_mask );
      fprintf(stderr,"\nNo smoothing has been performed."); 
      fprintf(stderr,"\nEnd of %s.\n", argv[0]); 
      exit(0);
   } /* end if */
 
  
/* set up work grid and weights*/
 
   fprintf(stderr,"\nNow smoothing with ellipse radius %1d/%1d pts ...", xradius, yradius); 
   
   hwork.nx = xradius * 2 + 1;
   hwork.ny = yradius * 2 + 1; 
   nwsq = hwork.nx * hwork.ny;
   hwork.x_inc = h.x_inc;
   hwork.y_inc = h.y_inc;
   hwork.node_offset = 0;  /* work grid is always gridline registered */
   
   zwork = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
   wsave = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
   get_weights(wsave, &hwork, xradius, yradius);
   
 
  /*visit each square in big grid, smooth nodes with actual values and output x,y,z triplets */
  
  for (i = 0; i < h.nx; ++i)  {
     for (j = 0; j < h.ny; ++j ) {
	sq = j * h.nx + i;
	
	if (z[sq] < flagged) {
	   ++n_empty;
	}
	else if (z[sq] < masked ) {
	   /* this node not empty or masked */
           error = ij2xy(&h, i, j, &x, &y);
	   hwork.x_min =  x - xradius * hwork.x_inc;
	   hwork.x_max =  x + xradius * hwork.x_inc;
	   hwork.y_min =  y - yradius * hwork.y_inc;
	   hwork.y_max =  y + yradius * hwork.y_inc;
           weights = (double *) get_memory((void *) NULL, (size_t)nwsq, sizeof(double));
	
	   for (wrow = 0; wrow < hwork.ny; ++wrow) {
	       for (wcol = 0; wcol < hwork.nx; ++wcol) {
	          wsq = wrow * hwork.nx + wcol;
	          weights[wsq] = wsave[wsq];
	          error = ij2xy(&hwork, wcol, wrow, &wx, &wy);
	          error = xy2ij(&h, wx, wy, &ix, &iy);
	          if (error < 0)
		     zwork[wsq] = mask_val;
	          else 
		     zwork[wsq] = z[iy * h.nx + ix];
	       }  /* for wcol */
	   } /* for wrow */
	
           wsq = yradius * hwork.nx + xradius;
	   if (ABS(zwork[wsq])< masked){
	      zero_weights_2d(zwork, weights, empty_val, mask_val, &hwork, xradius, yradius);
	      zout = weighted_mean(zwork, weights, empty_val, mask_val, nwsq);
	   
	      if (ABS(zout) < masked)
	         fprintf(outfile,"%lf %lf %lf\n", x, y, zout);
	   }
	   free(weights);
	 } /* end else if */  
     }  /* end for j */
  }  /* end for i*/

  fprintf(stderr,"\nNumber read in: <%d>  filled: <%d>   empty: <%d>   masked: <%d>\n", nread, n_set, n_empty, n_mask );
  fprintf(stderr,"\nEnd of %s.\n", argv[0]);
  
  exit(0);
  
}


/******************************************************************/
void print_usage(char *program)
{
  
    fprintf (stderr, "%s - Interpolates vertical section profiles onto a regular grid and optionally smoothes the 2D field using distance-weighted gaussian averaging with variable x/y radii \n\n", program);
    fprintf(stderr, "USAGE: %s [xyzfile(s)] -B<xmin/xmax/ymin/ymax> -I<dx>[/<dy>]", program);
    fprintf(stderr, "-L<xlen>[/<ylen>] [-M<mask_file> ] [-O<output_file>] [-P]");
    fprintf(stderr, "  [-S<xradius>/<yradius>]  [-:] [-h]\n\n");
    fprintf(stderr, "-B   sets the grid bounds in user units: xmin/xmax/ymin/ymax.\n");
    /*    fprintf(stderr, "-E   Allows horizontal extrapolation if only one grid node is found.\n");
	  fprintf(stderr, "     If not specified, then only interpolation is allowed.\n");*/
    fprintf(stderr, "-I   sets the grid spacing for the x/y directions.\n");
    fprintf(stderr, "-L   set x,y search radii in integer gridnodes \n");
    fprintf(stderr, "      If no data are within this range of a node, it is set to empty. \n");
    fprintf(stderr, "\n\tOPTIONS:\n");
    fprintf(stderr, "-M name of file containing masking info. These points will not be filled.\n");
    fprintf(stderr, "      Each line of file may specify individual points to mask\n");
    fprintf(stderr, "      or multiple polygons separated by a '>' character. \n");
    fprintf(stderr, "      Specify 'I' or 'O' immediately after the '>' character\n");
    fprintf(stderr, "      to mask INSIDE or OUTSIDE of polygon.\n");
    fprintf(stderr, "-O  name of output file  (default is stdout).\n");
    fprintf(stderr, "-P  force pixel registration (default is gridnode registration).\n");
    fprintf(stderr, "-S set xradius/yradius of smoothing ellipse (integer grid nodes)\n");
    fprintf(stderr, "      default is no smoothing: [-S%d/%d].\n", xradius, yradius);
 
    fprintf(stderr, "-: input data are ordered y x z  \n");
    fprintf(stderr, "       [default order is x y z]   \n");
    fprintf(stderr, "-h help....prints this message.  \n");
    return;
}  /* end print_usage() */

/**********************************************************/
void get_weights(double *weight, struct GRID_INFO *hptr, int xcntr, int ycntr)
  /* sets weight based on distance of grid_node from center node
        weight[n] = e^[- (pi/halfwidth)^2 /4.5 * dist^2 ].  
*/
{

int ix, iy, sq, dx, dy;
double phi, dist, halfwidth, alpha, pi;

   pi = 3.141592654;
   halfwidth = xcntr;
   if (ycntr > xcntr)
      halfwidth = ycntr;
      
   alpha = pi / halfwidth;
   alpha  *= alpha;
   alpha /= 4.5;
  
   for (ix = 0; ix < hptr->nx; ++ix) {
      for (iy = 0; iy < hptr->ny; ++iy) {
          sq = iy * hptr->nx + ix;
	  dy = iy-ycntr;
	  dx = ix-xcntr;
	  dist = sqrt(dy*dy + dx*dx);
	  weight[sq] = exp( -alpha * dist * dist);
      }
   }

   return;
} /* end get_weights() */
