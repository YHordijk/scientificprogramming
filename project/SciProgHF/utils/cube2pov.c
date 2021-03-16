/*
 *     Copyright (c) 2018 by the authors of DIRAC.
 *     All Rights Reserved.
 *
 *     This source code is part of the DIRAC program package.
 *     It is provided under a written license and may be used,
 *     copied, transmitted, or stored only in accordance to the
 *     conditions of that written license.
 *
 *     In particular, no part of the source code or compiled modules may
 *     be distributed outside the research group of the license holder.
 *     This means also that persons (e.g. post-docs) leaving the research
 *     group of the license holder may not take any part of Dirac,
 *     including modified files, with him/her, unless that person has
 *     obtained his/her own license.
 *
 *     For information on how to get a license, as well as the
 *     author list and the complete list of contributors to the
 *     DIRAC program, see: http://www.diracprogram.org
 */

/*
  Read a gaussian cube file and output povray files to render an isosurface.
  Ulf Ekstrom 2008
*/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <limits.h>

struct grid
{
  int size[3];
  double *data;
  double origin[3];
  double frame[3][3];
};

void write_bigendian_uint16(FILE *dst, unsigned short int value)
{
  unsigned char *data = (unsigned char *)&value;
#ifdef BIG_ENDIAN
#error Big endian not supported right now. If you are on i386 with gcc make sure you use -std=c99
#else
  fwrite(&data[1],1,1,dst);
  fwrite(&data[0],1,1,dst);
#endif
}



void write_bigendian_uint32(FILE *dst, unsigned int value)
{
  unsigned char *data = (unsigned char *)&value;
#ifdef BIG_ENDIAN
#error Big endian not supported right now. If you are on i386 with gcc make sure you use -std=c99
#else
  fwrite(&data[3],1,1,dst);
  fwrite(&data[2],1,1,dst);
  fwrite(&data[1],1,1,dst);
  fwrite(&data[0],1,1,dst);
#endif
}

void write_bigendian_uint32_float(FILE *dst, double x)
{
  if (x < 0)
    x = 0;
  else if (x > 1)
    x = 1;
  write_bigendian_uint32(dst,UINT_MAX*x);
}

void skip_line(FILE *src)
{
  while ((!feof(src)) && (fgetc(src) != '\n'))
    ;
}

struct grid *read_cube_data(FILE *src)
{
  int res, natom, i, j, np;
  struct grid *cube;
  cube = malloc(sizeof*cube);
  skip_line(src);
  skip_line(src);
  assert(!feof(src) && "Unexpected end of file looking for nr atoms");
  res = fscanf(src,"%i",&natom);
  assert(res == 1 && "Error reading nr atoms");
  fprintf(stderr,"Nr atoms: %i\n",natom);
  for (i=0;i<3;i++)
    {
      res = fscanf(src,"%lf",&cube->origin[i]);
      assert(res == 1 && "Error reading origin");
    }
  skip_line(src);
  for (i=0;i<3;i++)
    {
      res = fscanf(src,"%i",&cube->size[i]);
      assert(res == 1 && "Error reading grid size");
      for (j=0;j<3;j++)
	{
	  res = fscanf(src,"%lf",&cube->frame[i][j]);
	  assert(res == 1 && "Error reading frame");
	}
      skip_line(src);
    }
  for (i=0;i<natom;i++)
    skip_line(src);
  np = cube->size[0]*cube->size[1]*cube->size[2];
  cube->data = malloc(sizeof*cube->data*np);
  assert(cube->data && "Out of memory");
  for (i=0;i<np;i++)
    {
      res = fscanf(src,"%lf",&cube->data[i]);
      assert(res == 1 && "Error reading grid value");
    }
  return cube;
}

char *strsuffix(const char *s, const char *oldsuffix, const char *newsuffix)
{
  char *ns;
  int len = strlen(s), lensuf = strlen(oldsuffix);
  if (lensuf <= len && strcmp(oldsuffix,s+len-lensuf) == 0)
    {
      ns = malloc(len-lensuf+strlen(newsuffix)+1);
      strcpy(ns,s);
      strcpy(ns+len-lensuf,newsuffix);
      return ns;
    }
  else
    {
      ns = malloc(len+strlen(newsuffix)+1);
      strcpy(ns,s);
      strcpy(ns+len,newsuffix);
    }
  return ns;
}

double dgetopt(int argc, char * const argv[], const char *opt, double defaul)
{
  int i;
  for (i=0;i<argc;i++)
    {
      char *s = strstr(argv[i],opt);
      if (s == argv[i])
	{
	  double val;
	  int res = sscanf(argv[i]+strlen(opt),"%lf",&val);
	  if (res)
	    return val;
	  else
	    fprintf(stderr,"Error parsing option %s, using default\n",opt);
	}
    }
  return defaul;
}

double normalize_on_01(double x, double x0, double k)
{
  return atan((x-x0)*k)/3.14159265358979323 + 0.5;
}

int main(int argc, char *argv[])
{
  double datamax,isovalue,avg;
  int i, np, nclamp, nlower;
  const char *cubepath;
  const char *df3path;
  const char *povpath;
  struct grid *grid;
  FILE *gridfile, *df3file, *povfile;

  if (argc < 2)
    {
      fprintf(stderr,"Usage: cube2iso [--value=ISOVALUE] CUBEFILE\n");
      return EXIT_FAILURE;
    }
  cubepath = argv[argc-1];
  df3path = strsuffix(cubepath,".cube",".df3");
  povpath = strsuffix(cubepath,".cube",".pov");
  isovalue = dgetopt(argc,argv,"--value=",0.1);
  printf("isovalue=%lf\n",isovalue);
  
  gridfile = fopen(cubepath,"r");
  assert(gridfile && "Error opening cube file");
  grid = read_cube_data(gridfile);
  fclose(gridfile);
  assert(grid && "Error reading cube file");

  fprintf(stderr,"Normalizing data.\n");
  np = grid->size[0]*grid->size[1]*grid->size[2];
  datamax = 0;
  nclamp = 0;
  avg = 0;
  nlower = 0;
  for (i=0;i<np;i++)
    {
      if (grid->data[i] > datamax)
	datamax = grid->data[i];
      if (grid->data[i] < 0)
	{
	  grid->data[i] = 0;
	  nclamp++;
	}

      /* Normalizing to put the iso value at 0.5 -> later at ~16000*/
      if (grid->data[i] < isovalue)
	nlower ++;
      grid->data[i] = normalize_on_01(grid->data[i], isovalue, 1000);
      avg += grid->data[i];
    }
  fprintf(stderr,"Largest value before normalization: %lf\n",datamax);
  if (nclamp > 0)
    fprintf(stderr,"Warning, %i negative gridpoints set to 0.\n",nclamp);
  fprintf(stderr,"Average output density: %lf\n",avg/np);
  fprintf(stderr,"Total of %.2g%% of points smaller than iso value.\n",nlower/(double)np*100);

  df3file = fopen(df3path,"w");
  assert(df3file && "Could not open df3 file for writing");
  for (i=0;i<3;i++)
    write_bigendian_uint16(df3file,grid->size[i]);
  for (i=0;i<np;i++)
    write_bigendian_uint32_float(df3file,grid->data[i]);
  fclose(df3file);

  povfile = fopen(povpath,"w");
  fprintf(povfile,"#include \"colors.inc\"\n"
"camera{\n"
"	location <3, 1, -20>\n"
"	look_at <0,0,0>\n"
"	angle 45\n"
"}\n"
"light_source{ <-4, 10, -10> color rgb <1, 1, 1> shadowless}\n"
	  "background{ color rgb <0.7,0.7,1> }\n");
  fprintf(povfile,"#declare isofun=function{ pattern{\n"
	  "density_file df3 \"%s\" interpolate 3}}\n",df3path);

  fprintf(povfile,"isosurface{\n"
	  "function{ 1 - isofun(x,y,z) } \n"
	  "threshold 0.5\n"
	  "contained_by{ box{ 0,1 }}\n"
	  "max_gradient 5\n"
	  "texture{ pigment{ color White }} scale 3 }\n");
  fclose(povfile);

  if (errno)
    {
      perror("Error");
      return EXIT_FAILURE;
    }
  else
    {
      return EXIT_SUCCESS;
    }
}
