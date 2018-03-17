/**
 *   \file histo.c
 *   \brief A Documented file.
 *
 *  Detailed description
 *  This programs creates subsequences from a given ADN sequences file
 *  
 *   \Author: Danny Múnera - Parallel Programing UdeA 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_SQ 5000
#define MAX_LINE 100

//#define DEBUG

void process_sq (char* sq, size_t sq_len, int num_items, char* histogram);
size_t get_index(char* sq, size_t sz);
  
int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      fprintf(stderr, "ERROR - usage: histo <file> num_items <outfile>\n");
      exit(1);
    }
  char in_file[20];
  char out_file[20];
  int num_items;
  char sq_buffer[MAX_SQ], temp_buf[MAX_LINE];
  size_t sq_len, ln_len, i;
  
  strcpy(in_file, argv[1]);
  num_items = strtol(argv[2], NULL, 10);
  strcpy(out_file, argv[3]);
  
  // create vector
  // using (8-bits)characters to keep the frequency of the histogram
  size_t max_ent = pow(4, num_items);
  char* histogram = (char *) calloc(max_ent, sizeof(char));

  /* Open file to load a sequence */
  FILE *infp = fopen(in_file, "r");
  int n_seq = 0;
  while (fgets(temp_buf,MAX_LINE,infp) != NULL)
    {
      // check if line is a properties line (strating with '>')
      if(temp_buf[0] == '>')
	{
	  if(n_seq > 0)
	    process_sq(sq_buffer, sq_len, num_items, histogram);
	  strcpy(sq_buffer,"");
	  sq_len = 0;
	  n_seq++;
	  continue;
	}
      // to not include newline char
      ln_len = strlen(temp_buf) - 1;
      sq_len += ln_len; 
      strncat(sq_buffer, temp_buf, ln_len);
      
      //sq_len = strlen(sq_buffer);
#     ifdef DEBUG  
      printf("Sequence of size %ld is %s \n", sq_len, sq_buffer);
#     endif
    } // End While
  fclose(infp);
  
  // create an output file
  FILE *outfp = fopen(out_file, "w");
  char fq;
  for (i = 0; i < max_ent; i++)
    {
      if((fq = histogram[i])!=0)
	fprintf(outfp,"%ld\t%d\n", i, fq);
    }
  fclose(outfp);
  free(histogram);
  return 0;
}


void process_sq (char* sq, size_t sq_len, int num_items, char* histogram)
{
  int i;
  char sub_sq[num_items + 1]; // including '\0' char
  for(i = 0; i <= sq_len - num_items; i++)
    {
      memcpy(sub_sq, &sq[i], num_items);
      sub_sq[num_items] = '\0';
      size_t in = get_index(sub_sq, num_items );
      histogram[in]++;
#     ifdef DEBUG
      printf("sub sq %s , index = 0x%.8lX \n",sub_sq, in);
#     endif
    }
}

size_t get_index(char* sq, size_t sz)
{
  size_t i, index = 0;
  for(i = 0; i < sz; i++)
    {
      // compare string from last position to first
      //printf("char i %ld = %c\n",i, sq[sz-(i+1)] );
      switch (sq[sz-(i+1)]) {
      case 'A': 
	//index += 0;   
	break;
      case 'C':
	index += 1 << 2 * i;   
	break;
      case 'G': 
	index += 2 << 2 * i;      
	break;
      case 'T': 
	index += 3 << 2 * i;
	break;
      default:
	break;
      }
    }
  return index;  
}
