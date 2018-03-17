/**
 *   \file super_seq.c
 *   \brief Takes a fsn or fasta file and returns an outfile with 
 *          with a single line concatening all sequences
 *
 *  Detailed description
 *  
 *   \Author: Danny Múnera - Parallel Programing UdeA 
 *   \Date: 17/3/2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SQ 5000
#define MAX_LINE 100

//#define DEBUG

int main(int argc, char *argv[])
{
  if (argc != 3)
    {
      fprintf(stderr, "ERROR - usage: superseq <file> <outfile>\n");
      exit(1);
    }
  char in_file[20];
  char out_file[20];
  //char sq_buffer[MAX_SQ];
  char temp_buf[MAX_LINE];
  size_t t_len, ln_len;
  
  strcpy(in_file, argv[1]);
  strcpy(out_file, argv[2]);
  
  // Open infile to load sequences 
  FILE *infp = fopen(in_file, "r");
  // Open outfile
  FILE *outfp = fopen(out_file, "w");
  t_len=0;
  while (fgets(temp_buf,MAX_LINE,infp) != NULL)
    {
      // check if line is a properties line (strating with '>')
      if(temp_buf[0] == '>')
	continue;
      
      // to not include newline char
      ln_len = strlen(temp_buf) - 1;
      t_len += ln_len;

      fwrite(temp_buf, 1, ln_len, outfp);
      //strncpy(sq_buffer, temp_buf, ln_len);
      //fprintf(outfp,"%s",sq_buffer);
    } // End While

  //print a new line at the end of file
  fprintf(outfp,"\n");
  fclose(infp);
  fclose(outfp);

  printf("total size of the super seq: %ld\n", t_len);
  
  return 0;
}
