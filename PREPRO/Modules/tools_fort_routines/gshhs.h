/* Include file defining structures used in gshhs.c  Paul Wessel, SOEST */

#define _POSIX_SOURCE 1		/* GSHHS code is POSIX compliant */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* For byte swapping if needed */

#define swabi2(i2) (((i2) >> 8) + (((i2) & 255) << 8))
#define swabi4(i4) (((i4) >> 24) + (((i4) >> 8) & 65280) + (((i4) & 65280) << 8) + (((i4) & 255) << 24))

struct GSHHS {	/* Global Self-consistant Hierarchical High-resolution Shorelines */
	int id;				/* Unique polygon id number, starting at 0 */
	int n;				/* Number of points in this polygon */
        int flag;  /* = level + version << 8 + greenwich << 16 + source << 24 + river << 25 + p << 26 */
        /* flag contains 6 items, as follows:
         * low byte:    level = flag & 255: Values: 1 land, 2 lake, 3 island_in_lake,
         *                                      4 pond_in_island_in_lake
         * 2nd byte:    version = (flag >> 8) & 255: Values: Should be 9 for GSHHG release 9.
         * 3rd byte:    greenwich = (flag >> 16) & 3: Values: 0 if Greenwich nor Dateline are crossed,
         *              1 if Greenwich is crossed, 2 if Dateline is crossed, 3 if both is crossed.
         * 4th byte:    source = (flag >> 24) & 1: Values: 0 = CIA WDBII, 1 = WVS
         * 4th byte:    river = (flag >> 25) & 1: Values: 0 = not set, 1 = river-lake and
         *                                       GSHHG level = 2 (or WDBII class 0)
         * 4th byte:    area magnitude scale p (as in 10^p) = flag >> 26.  We divide area by 10^p.
         */
	int west, east, south, north;	/* min/max extent in micro-degrees */
	int area;			/* Area of polygon in 1/10 km^2 */
        int area_full;  /* Area of original full-resolution polygon in 1/10 km^2 */
        int container;  /* Id of container polygon that encloses this polygon (-1 if none) */
        int ancestor;   /* Id of ancestor polygon in the full resolution set that was
                                         the source of this polygon (-1 if none) */
};

struct	POINT {
	int	x;
	int	y;
};
