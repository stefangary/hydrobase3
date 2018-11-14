Some notes about the differences between HB2 and HB3.
Also included here are notes about little bits of
experience picked up with both HB2 and HB3.

I added the ./HB3/scripts folder to keep all my HB3 scripts
in the same place.

1) The isopycnal bins in lists/*.sig?list are the same
   from HB2 to HB3.  Direct comparison revealed an
   error in 7501.sig0list in HB2 (extra space).  This
   has been removed in the HB2 version that sits at
   ~/src/HB2/lists

2) hb_grid3d in HB2 has been essentially replaced with
   hb_bin3d in HB3.  Similar syntax, looks like major
   addition is 3D variance in addition to 3D obsevation
   counts.  Will also split data into monthly bins
   for creating a monthly seasonal climatology.  This
   feeds into the objective interpolation routines
   later that work in space and time.

3) hb_fit3d in HB2 has been essentially replaced with
   hb_ncfg3d in HB3.  The -U option () is not available
   in hb_ncfg3d because the topography file in hb_ncfg3d
   is now a netcdf .grd file, so that information is
   embedded.  Also, looks like the -W option in hb_fit3d
   has been removed and replaced with an automatic look
   up to whatever topography file is specified in the
   hb_paths.h header.  Updated the header to reflect
   local installation, see below.

4) The topography is roughly an order of magnitude
   higher resolution in HB3.

5) hb_smooth3d is replaced with hb_ncsmooth3d.  Looks
   like a direct one-to-one replacement, all the options
   appear to be the same.

6) If the variance of a variable is present in the
   netCDF file, then the way to extract a section of
   the variance is with -PVsa/sa to get both the
   salinity and its variance.  The gridding will
   simply not allow ht, s0, etc. (and gives warnings
   for it to - just skips over those derived values) 
   so the only variables we can get errors for are
   the basic ones - pr, sa, te, etc.

7) HB2 had only T68 but HB3 now uses both T68 and
   T90.  The CISRO sw routines use T90, so be certain
   to specify t90 and th9 for output.  During hb_bin3d,
   can only have one temperature variable, so work
   with T90.

Modifications for installation:

1) Specified the netcdf 3.6.3 library and include
   file locations in the makefile.

2) Specified the topography file location in
   hb_paths.h along with all other paths to
   sigma lists, etc.  Had to manually delete
   all the binaries, make clean, and then
   go through make and make install.

3) Something is wrong with -Pht - refuses to
   accept reference levels for dynamic height.
   For now, skip this and just provide T and S
   gridded profiles for the transport calculations.

4) Added ability to update station number, lon, and
   lat to hb_updatehdr.  Clarified the documentation
   built into that code (what is displayed with -h).

=======================================================

Some small experiences to keep in mind:
1) When making any manual edits to HB files, very
   important to ensure that there are the right number
   of white spaces around each number.  For example,
   I changed a GPS coordinate for a bad glider GPS
   fix and then when running hb_counsta, got
   "stack smashing" errors -> this was because I had
   changed a coordinate from -12 to -9 and neglected
   to include an extra bit of white space in front
   of the -9.

2) hb_prseries has a bug in it that may allow the newly
   created pressure series to overrun the deepest pressure,
   thus causing spurious HB_MISSING values because they
   outside the range of the interpolation.  Now fixed.
   Also, note that hb_prseries does not check for HB_MISSING
   during the interpolation proceedure, so it may be possible
   to include HB_MISSING values during the interpolation!!!
   This is different from hb_interpolate (HB2) which does
   explicitly filter out any HB_MISSING values.


