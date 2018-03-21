/**
 *   \file histo-vector.c
 *   \brief Creates a histogram from a "fna" or "fasta" file.
 *
 *  Detailed description
 *  This program creates an histogram form a "fna" or "fasta" file 
 *  using a vector to represent the histogram array.  
 *  
 *  Compile: gcc -Wall -o histo-vector histo-vector.c -lm
 *  Usage: ./histo-vector Test_Bancomini.fna 15 out.dat
 *
 *   \Author: Danny Múnera - Parallel Programing UdeA 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define MAX_SQ 5000
#define MAX_LINE 100

//#define DEBUG

void process_all_sq (char** all, size_t sq_num, int k_mers, unsigned short* histogram);
size_t get_index(char* sq, size_t sz);
  
int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      fprintf(stderr, "ERROR - usage: histo <file> k_mers <outfile>\n");
      exit(1);
    }
  char in_file[200];
  char out_file[200];
  int k_mers;
  char sq_buffer[MAX_SQ], temp_buf[MAX_LINE];
  size_t sq_len, ln_len, i;
  struct timeval t1, t2;
  double elapsedTime;  
  
  strcpy(in_file, argv[1]);
  k_mers = strtol(argv[2], NULL, 10);
  strcpy(out_file, argv[3]);
  
  // create vector
  // using (8-bits)characters to keep the frequency of the histogram
  size_t max_ent = pow(4, k_mers);
  unsigned short* histogram = (unsigned short*) calloc (max_ent, sizeof(unsigned short));
  if(histogram == NULL)
    {
      fprintf(stderr, "Calloc error while assigning memory to vector\n");
      exit(1);
    }

  // Data structure for sequences
  int all_sq_sz = MAX_SQ;
  char** all_sq = (char **) malloc(all_sq_sz * sizeof(char*));
  if(all_sq == NULL)
    {
      fprintf(stderr, "Malloc error while assigning memory to seq array\n");
      exit(1);
    }
  
  /* Open file to load a sequence */
  FILE *infp = fopen(in_file, "r");
  int n_seq = 0;
  
  while (fgets(temp_buf,MAX_LINE,infp) != NULL)
    {
      // check if line is a properties line (strating with '>')
      if(temp_buf[0] == '>')
	{
	  if(n_seq > 0)
	    {
#             ifdef DEBUG  
	      printf("Sequence %d of size %ld is %s \n", n_seq, sq_len, sq_buffer);
#             endif
	      all_sq[n_seq - 1] = (char*) malloc ((sq_len + 1)*sizeof(char));
	      if(all_sq[n_seq - 1] == NULL)
		{
		  fprintf(stderr, "Malloc error while assigning memory to char array\n");
		  exit(1);
		}
	      strcpy(all_sq[n_seq - 1], sq_buffer);
	      if ( n_seq % MAX_SQ == 0)
		{
		  all_sq_sz += MAX_SQ;
		  all_sq = realloc(all_sq,(all_sq_sz*sizeof(char*)));
		  if(all_sq == NULL)
		    {
		      fprintf(stderr, "Realloc error while re-assigning memory to vector\n");
		      exit(1);
		    }
		}		   		
	    }
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

    } // End While
  
      // save last sequence
# ifdef DEBUG  
  printf("Sequence %d of size %ld is %s \n", n_seq, sq_len, sq_buffer);
# endif
  all_sq[n_seq - 1] = (char*) malloc ((sq_len + 1)*sizeof(char));
  strcpy(all_sq[n_seq - 1], sq_buffer);
  if ( n_seq % MAX_SQ == 0)
    {
      all_sq_sz += MAX_SQ;
      all_sq = (char **) realloc (all_sq, (all_sq_sz*sizeof(char*)));
    }		  
  fclose(infp);

  // process all sequences
  gettimeofday(&t1, NULL);
  process_all_sq (all_sq, n_seq, k_mers, histogram);
  gettimeofday(&t2, NULL);
  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms}
  printf("Processing time: %5.3f ms\n", elapsedTime);
  
  // create an output file
  FILE *outfp = fopen(out_file, "w");
  char fq;
  for (i = 0; i < max_ent; i++)
    {
      if((fq = histogram[i])!=0)
	fprintf(outfp,"0x%.8lX\t%hu\n", i, fq);
    }
  fclose(outfp);
  free(histogram);
  return 0;
}

void process_all_sq (char** all, size_t sq_num, int k_mers, unsigned short* histogram)
{
  int i, j, sq_len;
  unsigned long in;
  for(i = 0; i < sq_num; i++)
    {
      sq_len = strlen(all[i]);
      char sub_sq[k_mers + 1]; // including '\0' char
      for(j = 0; j <= sq_len - k_mers; j++)
	{
	  memcpy(sub_sq, &all[i][j], k_mers);
	  sub_sq[k_mers] = '\0';
	  in = get_index(sub_sq, k_mers);
	  histogram[in]++;
#     ifdef DEBUG
	  printf("sub sq %s , index = 0x%.8lX \n",sub_sq, in);
#     endif
	} 
    } 
}


unsigned long get_index(char* sq, size_t sz)
{
  unsigned long i, index = 0;
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
