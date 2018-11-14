/*  hb_sigbins.h
.
.  definitions of structures and sigma bins used for least squares statistical
.  check of hydrographic data.
.
*/

struct binrec {
   int    nxy;
   double xsum;
   double ysum;
   double xysum;
   double xxsum;
   double yysum;
};

