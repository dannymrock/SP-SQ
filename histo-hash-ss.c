/**
 *   \file histo-hash-ss.c
 *   \brief A Documented file.
 *
 *  Detailed description
 *  This programs creates subsequences from a given ADN sequences file
 *  Get histogram info into a hash map structure
 *  Using Pete Warden simple hashmap implementation
 *     - http://petewarden.typepad.com/
 *     - https://github.com/petewarden/c_hashmap
 *
 *   Compile: gcc -Wall -o histo-hash-ss histo-hash-ss.c hashmap.o -lm
 *   Use:  ./histo-hash-ss Bancomini.dat 31 out.dat
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

#define MAX_LINE 100

//#define DEBUG

#define KEY_MAX_LENGTH (256)
#define KEY_COUNT (1024*1024)

typedef struct mapent_s
{
    char key_string[KEY_MAX_LENGTH];
    int number;
} mapent_t;

map_t mymap;

void process_supersq (char* sq, size_t sq_len, int num_items);
size_t get_index(char* sq, size_t sz);
//Function used for hashmap iterator to write data to outfile
int printent(void* fd, void * data);
  
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
  char* super_seq;
  struct timeval t1, t2;
  double elapsedTime;

  mymap = hashmap_new();
  
  strcpy(in_file, argv[1]);
  num_items = strtol(argv[2], NULL, 10);
  strcpy(out_file, argv[3]);
  

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
 
  process_supersq(super_seq, size, num_items);

  gettimeofday(&t2, NULL);
  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

  
  printf("Processing time: %5.3f ms\n", elapsedTime);

  // create an output file
  FILE *outfp = fopen(out_file, "w");
  hashmap_iterate(mymap, &printent, outfp);
  fclose(outfp);

  // Destroy the map 
  hashmap_free(mymap);
  
  return 0;
}


void process_supersq (char* sq, size_t sq_len, int num_items)
{
  int i;
  char sub_sq[num_items + 1]; // including '\0' char
  for(i = 0; i <= sq_len - num_items; i++)
    {
      memcpy(sub_sq, &sq[i], num_items);
      //End of the subsequence string
      sub_sq[num_items] = '\0';

      mapent_t* value; // = malloc(sizeof(data_struct_t));
      if (hashmap_get(mymap, sub_sq, (void**)(&value)) == MAP_MISSING)
	{
	  //printf("Map missing \n");
	  /* TODO: Possible memory leak due to following malloc
	   *       I need to check how hashmap is working
	   *       I suppose I need to iterate over each element and to free it
	   */
	  value = malloc(sizeof(mapent_t));
	  strcpy(value->key_string, sub_sq);
	  value->number=1;
	  int error = hashmap_put(mymap, value->key_string, value);
	  assert(error==MAP_OK);
	} else
	{
	  //printf("Map ok\n");
	  value->number++;
	}
      
#     ifdef DEBUG
      printf("sub sq %s \n",sub_sq);
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

int printent(void* fd, void* data)
{
  //printf("printing\n");
  fprintf((FILE *)fd,"%s\t%d\n", ((mapent_t*)data)->key_string, ((mapent_t*)data)->number);
  return MAP_OK;
}
