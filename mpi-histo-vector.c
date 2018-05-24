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
#include <mpi.h>
#include <assert.h>

#define MAX_SQ 5000
#define MAX_LINE 1000

//#define DEBUG

void process_all_sq (char** all, size_t sq_num, int k_mers,
		     unsigned int* histogram, int low, int high);
void get_index(char* sq, size_t sz, long long * index);
void get_char(char* sq, size_t sz, long long index);
  
int main(int argc, char *argv[])
{
  char in_file[200];
  char out_file[200];
  int k_mers;
  char sq_buffer[MAX_SQ], temp_buf[MAX_LINE];
  size_t sq_len, ln_len, i;
  //struct timeval t1, t2;
  //double elapsedTime;
  int c_size, myr;

  MPI_Init(&argc, &argv);
  MPI_Comm c = MPI_COMM_WORLD;
  MPI_Comm_size(c, &c_size);
  MPI_Comm_rank(c, &myr);

  if (argc != 4)
    {
      fprintf(stderr, "ERROR - usage: histo <file> k_mers <outfile>\n");
      exit(1);
    }

  strcpy(in_file, argv[1]);
  k_mers = strtol(argv[2], NULL, 10);
  strcpy(out_file, argv[3]);
  
  // Each process creates a vector
  // using unsigned ints to keep the frequency of the histogram
  long long max_ent = pow(4, k_mers);
  // each process reports my_ent entries to the histogram 
  long long my_ent = max_ent / c_size;
  unsigned int* histogram = (unsigned int*) calloc (my_ent,
						    sizeof(unsigned int));
  assert(histogram != NULL);

#ifdef DEBUG
  printf("myr: %d c_size: %d k_mers: %d max_ent: %lld my_ent: %lld\n", myr,
	 c_size, k_mers, max_ent, my_ent);
#endif // DEBUG
 
  char** all_sq;
  int all_sq_sz;
  int n_seq;
  int* sizes; 
  if(myr == 0)
    {
      // Data structure for sequences
      all_sq_sz = MAX_SQ;
      all_sq = (char **) malloc(all_sq_sz * sizeof(char*));
      sizes = (int *) malloc(all_sq_sz * sizeof(int));
      assert(all_sq != NULL);
      assert(sizes != NULL);
      
      /* Open file to load a sequence */
      FILE *infp = fopen(in_file, "r");
      assert(infp != NULL);
      n_seq = 0;
      while (fgets(temp_buf,MAX_LINE,infp) != NULL)
	{
	  // check if line is a properties line (strating with '>')
	  if(temp_buf[0] == '>')
	    {
	      if(n_seq > 0)
		{
#ifdef DEBUG  
		  printf("Sequence %d of size %ld is %s \n", n_seq,
			 sq_len, sq_buffer);
#endif
		  sizes[n_seq - 1] = sq_len + 1;
		  all_sq[n_seq - 1] = (char*) malloc ((sq_len + 1)
						      *sizeof(char));
		  assert(all_sq[n_seq - 1] != NULL);
		  
		  strcpy(all_sq[n_seq - 1], sq_buffer);
		  if ( n_seq % MAX_SQ == 0)
		    {
		      all_sq_sz += MAX_SQ;
		      all_sq = realloc(all_sq,(all_sq_sz*sizeof(char*)));
		      sizes = realloc(sizes,(all_sq_sz*sizeof(int)));
		      assert(all_sq != NULL);
		      assert(sizes != NULL);
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
	} // End While

      // save last sequence
# ifdef DEBUG  
      printf("Sequence %d of size %ld is %s \n", n_seq, sq_len, sq_buffer);
# endif
      sizes[n_seq - 1] = sq_len + 1;
      all_sq[n_seq - 1] = (char*) malloc ((sq_len + 1)*sizeof(char));
      strcpy(all_sq[n_seq - 1], sq_buffer);
      fclose(infp);
      // Broadcast size of the vector (all_sq_sz) 
      MPI_Bcast(&all_sq_sz, 1, MPI_INT, 0, c);
      // Broadcast number of sequences (nseq)
      MPI_Bcast(&n_seq, 1, MPI_INT, 0, c);
      // Broadcast sizes array
      MPI_Bcast(sizes, n_seq, MPI_INT, 0, c);
      for(i = 0; i < n_seq; i++)
	{
	  MPI_Bcast(all_sq[i], sizes[i], MPI_CHAR, 0, c);
	}
    } //end if -  process 0 loading sequences
  else
    {
      MPI_Bcast(&all_sq_sz, 1, MPI_INT, 0, c);
      MPI_Bcast(&n_seq, 1, MPI_INT, 0, c);
      all_sq = (char **) malloc(all_sq_sz * sizeof(char*));
      assert(all_sq != NULL);
      // Broadcast sizes array
      sizes = (int*) malloc(all_sq_sz * sizeof(int));
      assert(sizes != NULL);
#ifdef DEBUG
      printf("******myr: %d all_sq_sz %d n_seq %d\n", myr, all_sq_sz, n_seq);
#endif // DEBUG

      MPI_Bcast(sizes, n_seq, MPI_INT, 0, c);
      for(i = 0; i < n_seq; i++)
	{
#ifdef DEBUG
	  printf("----------myr: %d sizes[%ld] %d\n", myr, i, sizes[i]);
#endif // DEBUG

	  all_sq[i] = (char*) malloc ((sizes[i])*sizeof(char));
	  MPI_Bcast(all_sq[i], sizes[i], MPI_CHAR, 0, c);
	}
    }


  int my_low = myr * my_ent;
  int my_high = (myr+1) * my_ent;
#ifdef DEBUG
  printf("Process %d ready to process sequences\n", myr);
#endif // DEBUG
  
  // process all sequences
  // gettimeofday(&t1, NULL);
  process_all_sq (all_sq, n_seq, k_mers, histogram, my_low, my_high);
  //gettimeofday(&t2, NULL);
  //elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms}
  //printf("Processing time: %5.3f ms\n", elapsedTime);

  //Free data structure
   for(i = 0; i < n_seq; i++)
     {
      free(all_sq[i]);
    }
   free(all_sq);
   free(sizes);
  
   // create an output file for each process
   char par_file[100];
   sprintf(par_file,"out-%d.out",myr);
   FILE *outfp = fopen(par_file, "w");
   assert(outfp != NULL);
   unsigned short fq;
   char buff[100];
   long long index;
   for (index = 0LL; index < my_ent; index++)
     {
       if((fq = histogram[index])!=0)
	 {
	   get_char(buff, k_mers, index + my_low);
	   fprintf(outfp,"%s\t%u\n", buff, fq);
	 }	  
     }
   fclose(outfp);
   free(histogram);
   MPI_Finalize();
   
   return 0;
}

void process_all_sq (char** all, size_t sq_num, int k_mers,
		     unsigned int* histogram, int low, int high)
{
  int i, j, sq_len;
  long long in;
  for(i = 0; i < sq_num; i++)
    {
      sq_len = strlen(all[i]);
      char sub_sq[k_mers + 1]; // including '\0' char
      for(j = 0; j <= sq_len - k_mers; j++)
	{
	  memcpy(sub_sq, &all[i][j], k_mers);
	  sub_sq[k_mers] = '\0';
	  get_index(sub_sq, k_mers, &in);

	  // report index only if is in my process range
	  if(in >= low)
	    if(in < high)
	      histogram[in - low]++; // = *(histogram+in) + 1;
	  
#     ifdef DEBUG
	  printf("sub sq %s , index = %lld \n",sub_sq, in);
#     endif
	}
    } 
}

void get_index(char* sq, size_t sz, long long *index)
{
  long long i;
  *index = 0LL;
  for(i = 0; i < sz; i++)
    {
      // compare string from last position to first
      //printf("char i %ld = %c\n",i, sq[sz-(i+1)] );
      switch (sq[sz-(i+1)]) {
      case 'A': 
	//index += 0;   
	break;
      case 'C':
	*index += 1LL << 2LL * i;   
	break;
      case 'G': 
	*index += 2LL << 2LL * i;      
	break;
      case 'T': 
	*index += 3LL << 2LL * i;
	break;
      default:
	break;
      }
    }
}

void get_char(char* sq, size_t sz, long long index)
{
  long long mask, masked, value;
  long long i, in;
  for(i = 0; i < sz; i++)
    {
      mask = 3LL << (i*2LL);
      masked = index & mask;
      value = masked >> (i*2);
      in = sz - 1 - i;
      switch (value) {
      case 0: 
	sq[in] = 'A';  
	break;
      case 1:
	sq[in] = 'C';  
	break;
      case 2: 
	sq[in] = 'G';
	break;
      case 3: 
	sq[in] = 'T';
	break;
      default:
	sq[in] = 'N';
	break;
      }
    }
  sq[sz] = '\0';
# ifdef DEBUG
  printf("index = %lld -> sq %s \n", index, sq);
# endif
}
