/**
 *   \file histo-hash.c
 *   \brief Creates a histogram from a "fna" or "fasta" file.
 *
 *  Detailed description
 *  This program creates an histogram form a "fna" or "fasta" file 
 *  using a hashmap to represent the histogram array.  
 *  Using Pete Warden simple hashmap implementation
 *     - http://petewarden.typepad.com/
 *     - https://github.com/petewarden/c_hashmap
 *
 *   Compile: gcc -Wall -o histo-hash histo-hash1.c hashmap.o -lm
 *   Use:  ./histo-hash Bancomini.dat 31 out.dat
 *  
 *   \Author: Danny Múnera - Parallel Programing UdeA 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include "hashmap.h"

#define MAX_SQ 5000
#define MAX_LINE 100

#define KEY_MAX_LENGTH (256)
#define KEY_COUNT (1024*1024)

typedef struct mapent_s
{
    char key_string[KEY_MAX_LENGTH];
    int number;
} mapent_t;

// Global variable Hashmap 
map_t mymap;

//#define DEBUG

void process_all_sq (char** all, size_t sq_num, int k_mers);
int printent(void* fd, void * data);
  
int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      fprintf(stderr, "ERROR - usage: histo <file> k_mers <outfile>\n");
      exit(1);
    }
  char in_file[20];
  char out_file[20];
  int k_mers;
  char sq_buffer[MAX_SQ], temp_buf[MAX_LINE];
  size_t sq_len, ln_len;
  struct timeval t1, t2;
  double elapsedTime;  

  mymap = hashmap_new();
  
  strcpy(in_file, argv[1]);
  k_mers = strtol(argv[2], NULL, 10);
  strcpy(out_file, argv[3]);
  
  // Data structure for sequences
  int all_sq_sz = MAX_SQ;
  char** all_sq = (char **) malloc(all_sq_sz * sizeof(char*));
  
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
	      strcpy(all_sq[n_seq - 1], sq_buffer);
	      if ( n_seq % MAX_SQ == 0)
		{
		  all_sq_sz += MAX_SQ;
		  all_sq = realloc(all_sq,(all_sq_sz*sizeof(char*)));
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
  process_all_sq (all_sq, n_seq, k_mers);
  gettimeofday(&t2, NULL);
  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms}
  printf("Processing time: %5.3f ms\n", elapsedTime);
  
  // create an output file
  FILE *outfp = fopen(out_file, "w");
  hashmap_iterate(mymap, &printent, outfp);
  fclose(outfp);
  // Destroy the map 
  hashmap_free(mymap);
  return 0;
}

void process_all_sq (char** all, size_t sq_num, int k_mers)
{
  int i, j, sq_len;
  for(i = 0; i < sq_num; i++)
    {
      sq_len = strlen(all[i]);
      char sub_sq[k_mers + 1]; // including '\0' char
      for(j = 0; j <= sq_len - k_mers; j++)
	{
	  memcpy(sub_sq, &all[i][j], k_mers);
	  sub_sq[k_mers] = '\0';

	  mapent_t* value; // = malloc(sizeof(data_struct_t));
	  if (hashmap_get(mymap, sub_sq, (void**)(&value)) == MAP_MISSING)
	    {
	      //printf("Map missing \n");
	      value = malloc(sizeof(mapent_t));
	      strcpy(value->key_string, sub_sq);
	      value->number=1;
	      int error = hashmap_put(mymap, value->key_string, value);
	      assert(error==MAP_OK);
	    }
	  else
	    {
	      //printf("Map ok\n");
	      value->number++;
	    }	  
	  //size_t in = get_index(sub_sq, k_mers);
	  //histogram[in]++;
#     ifdef DEBUG
	  printf("sub sq %s , index = 0x%.8lX \n",sub_sq, in);
#     endif
	} 
    } 
}

int printent(void* fd, void* data)
{
  //printf("printing\n");
  fprintf((FILE *)fd,"%s\t%d\n", ((mapent_t*)data)->key_string, ((mapent_t*)data)->number);
  return MAP_OK;
}
