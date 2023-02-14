/* PURPOSE: To extract binary shoreline data file as described in the 1996
 Wessel & Smith JGR Data Analysis Note. Adapted from the original code gshhs.c
 written by Paul Wessel (wessel@soest.hawaii.edu) by converting it into a
 function which can be called from FORTRAN program and transmits all the data
 back into the calling program via arguments instead of creating an ASCII file.
*/

#include "gshhs.h"
char *trimwhitespace(char *str)
{
  char *end;

  // Trim leading space
  while(isspace((unsigned char)*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace((unsigned char)*end)) end--;

  // Write new null terminator character
  end[1] = '\0';

  return str;
}





int read_gshhs_(int *nsegm, int *isegm, struct GSHHS *hdr, struct POINT *p,char *name)
{
  FILE *fp;
  int j, k, m;
  size_t count;
  char *fname;

  fname=trimwhitespace(name);
  fprintf(stdout,"fname %s ... %ld .\n",fname,strlen(fname));
  if ((fp = fopen(fname,"rb")) == NULL ) {
     fprintf(stderr, "read_gshhs :: Cannot find file %s.\n", fname);
     return 1;
  }
  fprintf(stdout, "read_gshhs :: reading data file %s.\n", fname);

   j=0; m=0; isegm[0]=1;
		
  while (fread((void *)&hdr[j], (size_t)sizeof(struct GSHHS), (size_t)1, fp) == 1) {

#ifdef FLIP
    hdr[j].id = swabi4 ((unsigned int)hdr[j].id);
    hdr[j].n = swabi4 ((unsigned int)hdr[j].n);
    hdr[j].flag = swabi4 ((unsigned int)hdr[j].flag);
    hdr[j].west = swabi4 ((unsigned int)hdr[j].west);
    hdr[j].east = swabi4 ((unsigned int)hdr[j].east);
    hdr[j].south = swabi4 ((unsigned int)hdr[j].south);
    hdr[j].north = swabi4 ((unsigned int)hdr[j].north);
    hdr[j].area = swabi4 ((unsigned int)hdr[j].area);
    hdr[j].area_full = swabi4 ((unsigned int)hdr[j].area_full);
    hdr[j].container = swabi4 ((unsigned int)hdr[j].container);
    hdr[j].ancestor = swabi4 ((unsigned int)hdr[j].ancestor);
#endif

//   fprintf(stdout, "j = %d, m = %d, id = %d, n = %d\n", j,m,hdr[j].id, hdr[j].n);

     count=hdr[j].n;
     if (fread((void *)&p[m], (size_t)sizeof(struct POINT), count, fp) != count) {
       fprintf(stderr, "### ERROR: read_gshhs :: cannot read polygon %d from file %s.\n",
                                hdr[j].id, fname); return 1;
     }
#ifdef FLIP
     for (k = 0; k < hdr[j].n; k++) {
       p[m+k].x = swabi4((unsigned int)p[m+k].x);
       p[m+k].y = swabi4((unsigned int)p[m+k].y);
     }
#endif
   m=m+hdr[j].n; j++; isegm[j]=m+1; *nsegm=j;
  }
  fclose (fp);
  return 0;
}

int decode_flag_(int *flag, int*level, int *version, int *greenwich,
                        int *source, int  *river,  int *area_scale){

// flag = level + version << 8 + greenwich << 16 + source << 24 + river << 25 + p << 26
// hence contains 6 items, as follows:

    *level = *flag & 255 ;          // Values: 1 land, 2 lake, 3 island_in_lake,
                                    //           4 pond_in_island_in_lake
    *version = (*flag >> 8) & 255;  // Values: Should be 9 for GSHHG release 9.
    *greenwich = (*flag >> 16) & 3; // Values: 0 neither Greenwich nor Dateline are
                                    //     crossed; 1 if Greenwich is crossed; 2 if
                                    //    Dateline is crossed, 3 if both are crossed.
    *source = (*flag >> 24) & 1;    // Values: 0 = CIA WDBII, 1 = WVS
     *river = (*flag >> 25) & 1;    // Values: 0 = not set, 1 = river-lake and
                                    //        GSHHG level = 2 (or WDBII class 0)
    *area_scale = *flag >> 26;      //  area magnitude scale p (as in 10^p)
                                    //       We divide area by 10^p.
  return 0;
}
