#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#define BUFSIZE (1 << 20)

uint64_t read_marker(FILE *src, int marker_bytes)
{
  size_t siz;
  uint32_t u32;
  uint64_t u64;
  if (marker_bytes == 4)
    {      
      siz = fread(&u32,marker_bytes,1,src);
      u64 = u32;
    }
  else
    {
      siz = fread(&u64,marker_bytes,1,src);
    }
  assert(siz == 1);
  return u64;
}
void write_marker(FILE *dst, uint64_t marker, int marker_bytes)
{
  size_t siz;
  if (marker_bytes == 4)
    {
      uint32_t m = marker;
      siz = fwrite(&m,marker_bytes,1,dst);
      assert(siz == 1 && "32 bit case");
    }
  else
    {
      siz = fwrite(&marker,marker_bytes,1,dst);
      assert(siz == 1 && "64 bit case");
    }
}

void echo_bytes(FILE *in, FILE *out, size_t len)
{
  unsigned char buf[BUFSIZE];
  size_t s, siz;
  while (len > 0)
    {
      if (len > BUFSIZE)
	s = BUFSIZE;
      else
	s = len;
      siz = fread(buf,s,1,in);
      assert(siz = 1);
      siz = fwrite(buf,s,1,out);
      assert(siz = 1);
      len -= s;
    }
}

int main(int argc, const char *argv[])
{
  FILE *infile, *outfile;
  uint32_t i32;
  uint64_t i64, marker;
  int i;
  size_t siz;
  int ibytes, nfsym, marker_bytes, output_ibytes, output_marker_bytes, ncmotq, norbt;
  printf("DFCOEF conversion program, Ulf Ekstrom 2011\n");
  if (argc != 3)
    {
      fprintf(stderr,"Usage: fortran_convert INPUT OUTPUT\n");
      fprintf(stderr,"Reads input, a 32 or 64 bit fortran file, convert to the other bit format and and write result to OUTPUT\n");
      fprintf(stderr,"Endianness is currently _not_ converted, but could be implemented.\n");
      return EXIT_FAILURE;
    }
  infile = fopen(argv[1],"rb");
  if (!infile)
    {
      fprintf(stderr,"Cannot open %s for reading, quitting.\n",argv[1]);
      return EXIT_FAILURE;
    }
  outfile = fopen(argv[2],"wb");
  if (!outfile)
    {
      fprintf(stderr,"Cannot open %s for writing, quitting.\n",argv[2]);
      return EXIT_FAILURE;
    }

  
  i32 = read_marker(infile,4);
  rewind(infile);
  i64 = read_marker(infile,8);
  marker_bytes = 4;
  switch (i32)
    {
    case 98:
      ibytes = 4;
      nfsym = 1;
      break;
    case 110:
      ibytes = 4;
      nfsym = 2;
      break;
    case 114:
      ibytes = 8;
      nfsym = 1;
      break;
    case 138:
      ibytes = 8;
      nfsym = 2;
      break;
    default:
      marker_bytes = 8;
    }
  if (marker_bytes == 8)
    switch (i64)
      {
      case 98:
	ibytes = 4;
	nfsym = 1;
	break;
      case 110:
	ibytes = 4;
	nfsym = 2;
	break;
      case 114:
	ibytes = 8;
	nfsym = 1;
	break;
      case 138:
	ibytes = 8;
	nfsym = 2;
	break;
      default:
	fprintf(stderr,"Could not make any sense of the input, quitting\n");
	return EXIT_FAILURE;
    }
  if (ibytes == 8)
    output_ibytes = 4;
  else
    output_ibytes = 8;
  output_marker_bytes = marker_bytes;
  printf("Integer bits: %i, NFSYM: %i, marker bits: %i, output integer bits: %i\n",
	 ibytes*8,nfsym,marker_bytes*8,output_ibytes*8);
  /* Now do the conversion */
  
  rewind(infile);
  /* Header */
  marker = read_marker(infile,marker_bytes);
  write_marker(outfile,74+output_ibytes*(1+3*nfsym)+8,output_marker_bytes);
  echo_bytes(infile,outfile,74);
  for (i=0;i<1+3*nfsym;i++)
    write_marker(outfile,read_marker(infile,ibytes),output_ibytes);
  echo_bytes(infile,outfile,8);
  i64 = read_marker(infile,marker_bytes);
  printf("%llu %llu\n",marker,i64);
  assert(marker == i64 && "Header: Start marker not maching end marker, cannot proceed");
  write_marker(outfile,74+output_ibytes*(1+3*nfsym)+8,output_marker_bytes);

  /* Orbitals */
  marker = read_marker(infile,marker_bytes);
  write_marker(outfile,marker,output_marker_bytes);
  ncmotq = marker/8;
  printf("NCMOTQ: %i\n",ncmotq);
  echo_bytes(infile,outfile,marker);
  i64 = read_marker(infile,marker_bytes);
  assert(marker == i64 && "CMO: Start marker not maching end marker, cannot proceed");
  write_marker(outfile,i64,output_marker_bytes);


  /* Eigenvalues */
  marker = read_marker(infile,marker_bytes);
  write_marker(outfile,marker,output_marker_bytes);
  norbt = marker/8;
  printf("NORBT: %i\n",norbt);
  echo_bytes(infile,outfile,marker);
  i64 = read_marker(infile,marker_bytes);
  assert(marker == i64 && "EIG: Start marker not maching end marker, cannot proceed");
  write_marker(outfile,i64,output_marker_bytes);

  /* IBEIG */
  marker = read_marker(infile,marker_bytes);
  write_marker(outfile,norbt*output_ibytes,output_marker_bytes);
  assert(marker == norbt*ibytes);
  for (i=0;i<norbt;i++)
    write_marker(outfile,read_marker(infile,ibytes),output_ibytes);
  i64 = read_marker(infile,marker_bytes);
  assert(marker == i64 && "IBEIG: Start marker not maching end marker, cannot proceed");
  write_marker(outfile,norbt*output_ibytes,output_marker_bytes);

  return EXIT_SUCCESS;
}
