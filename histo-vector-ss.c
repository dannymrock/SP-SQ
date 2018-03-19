/**
 *   \file histo-vector-ss.c
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
#include <sys/time.h>


#define MAX_SQ 5000
#define MAX_LINE 100

//#define DEBUG

void process_supersq (char* sq, size_t sq_len, int num_items, char* histogram);
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
  char size_buf[MAX_LINE];
  //size_t sq_len, ln_len, i;
  size_t i;
  char* super_seq;
  struct timeval t1, t2;
  double elapsedTime;
  
  strcpy(in_file, argv[1]);
  num_items = strtol(argv[2], NULL, 10);
  strcpy(out_file, argv[3]);
  
  // create vector
  // using (8-bits) char type values to keep the frequency of the histogram
  size_t max_ent = pow(4, num_items);
  char* histogram = (char *) calloc (max_ent, sizeof(char));

  /* Open file to load a sequence */
  FILE *infp = fopen(in_file, "r");

  // Load first line with vector size
  fgets(size_buf, MAX_LINE, infp);
  int size = strtol(size_buf, NULL, 10);

# ifdef DEBUG
  printf("Size: %d\n", size);
# endif

  super_seq = (char *) malloc (size + 1);
  // Load second line 
  fgets(super_seq, size + 1, infp);
# ifdef DEBUG
  printf("Super_seq: %s\n", super_seq);
# endif
  fclose(infp);

  gettimeofday(&t1, NULL);
 
  process_supersq(super_seq, size, num_items, histogram);

  gettimeofday(&t2, NULL);
  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

  
  printf("Processing time: %5.3f ms\n", elapsedTime);

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


void process_supersq (char* sq, size_t sq_len, int num_items, char* histogram)
{
  int i;
  char sub_sq[num_items + 1]; // including '\0' char
  for(i = 0; i <= sq_len - num_items; i++)
    {
      memcpy(sub_sq, &sq[i], num_items);

      //End of the subsequence string
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
