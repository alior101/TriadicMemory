/*
triadicmemory.c

Copyright (c) 2022 Peter Overmann

C-language reference implementation of the Triadic Memory algorithm published in
   https://github.com/PeterOvermann/Writings/blob/main/TriadicMemory.pdf
This source file can be compiled as a stand-alone command line program or as a library.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the “Software”), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Building the command line version: cc -Ofast triadicmemory.c -o /usr/local/bin/triadicmemory

This command line tool instantiates a new memory instance. It can be used
to store triples {x,y,z} of sparse distributed representations (SDRs), and to
recall one part of a triple by specifying the other two parts.
{x,y,_} recalls the third part, {x,_,z} recalls the second part, {_,y,z} recalls
the first part.

An SDR is given by a set of p integers in the range from 1 to n.
Typical values are n = 1000 and p = 10.

Command line arguments: triadicmemory <n> <p>

Command line usage examples:

Store {x,y,z}:
{37 195 355 371 471 603 747 914 943 963, 73 252 418 439 461 469 620 625 902 922, 60 91 94 128 249 517 703 906 962 980}

Recall x:
{_ , 73 252 418 439 461 469 620 625 902 922,  60 91 94 128 249 517 703 906 962 980}

Recall y:
{37 195 355 371 471 603 747 914 943 963, _ , 160 91 94 128 249 517 703 906 962 980}

Recall z:
{37 195 355 371 471 603 747 914 943 963, 73 252 418 439 461 469 620 625 902 922, _}

Delete {x,y,z}:
-{37 195 355 371 471 603 747 914 943 963, 73 252 418 439 461 469 620 625 902 922, 60 91 94 128 249 517 703 906 962 980}

Print random SDR:
random

Print version number:
version

Terminate process:
quit

*/

// remove the following if compiling as library
#define TRIADICMEMORY_COMMANDLINE



// pull out the following if compiling as library

// -------------------- begin triadicmemory.h --------------------

#define MEMTYPE unsigned char

typedef struct
	{
	MEMTYPE* m;
	int n, p;
	} TriadicMemory;
	
typedef struct
	{
	int *a, n, p;
	} SDR;
	
SDR *sdr_random (SDR*, int n);
SDR *sdr_new (int n);
int sdr_distance( SDR*x, SDR*y);

TriadicMemory *triadicmemory_new(int n, int p);

void triadicmemory_write  (TriadicMemory *T, SDR *x, SDR *y, SDR *z);
void triadicmemory_delete (TriadicMemory *T, SDR *x, SDR *y, SDR *z);

SDR* triadicmemory_read_x  (TriadicMemory *T, SDR *x, SDR *y, SDR *z);
SDR* triadicmemory_read_y  (TriadicMemory *T, SDR *x, SDR *y, SDR *z);
SDR* triadicmemory_read_z  (TriadicMemory *T, SDR *x, SDR *y, SDR *z);
	
// -------------------- end triadicmemory.h --------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
/*
void swapElements(int * x, int* y) 
{ 
	int temp = *x; 
	*x = *y; 
	*y = temp; 
}   
// Partition function
int partition (int  arr[], int lowIndex, int highIndex) 
{ 
	int pivotElement = arr[highIndex];
	int i = (lowIndex - 1); 
	for (int j = lowIndex; j <= highIndex- 1; j++) 
	{ 
		if (arr[j] <= pivotElement) 
		{ 
			i++; 
			swapElements(&arr[i], &arr[j]); 
		} 
	} 
	swapElements(&arr[i + 1], &arr[highIndex]); 
	return (i + 1); 
}   
// QuickSort Function
void quickSort(int arr[], int lowIndex, int highIndex) 
{ 
	if (lowIndex < highIndex) 
	{ 
		int pivot = partition(arr, lowIndex, highIndex); 
		// Separately sort elements before & after partition 
		quickSort(arr, lowIndex, pivot - 1); 
		quickSort(arr, pivot + 1, highIndex); 
	} 
}   
// Function to print array
void printArray(int arr[], int size) 
{ 
	int i; 
	for (i=0; i < size; i++) 
	printf("%d ", arr[i]); 
}   
*/

static int cmpfunc (const void * a, const void * b)
	{ return  *(int*)a - *(int*)b; }
		
static SDR* binarize (SDR *v, int targetsparsity)
	{
	int sorted[v->n], rankedmax;
	
	for ( int i=0; i < v->n; i++ )
		sorted[i] = v->a[i];
	
	qsort( sorted, v->n, sizeof(int), cmpfunc);
	
	rankedmax = sorted[ v->n - targetsparsity ];
	
	if(rankedmax == 0)
		rankedmax = 1;
		
	v->p = 0;
	for ( int i = 0; i < v->n; i++)
		if (v->a[i] >= rankedmax)
			v->a[v->p++] = i;
	
	return (v);
	}


SDR* sdr_random( SDR*s, int p)
	{
	static int *range = 0;
	int n = s->n;
	s->p = p;
	
	if (! range)
		{
		srand((unsigned int)time(NULL));				// make sure this is only called once
		range = malloc(n*sizeof(int));					// integers 0 to n-1
		for (int i = 0; i < n; i++) range[i] = i;
		}
		
	for (int k = 0; k < p; k++) // random selection of p integers in the range 0 to n-1
		{
		int r = rand() % n;
		s->a[k] = range[r];
		int tmp = range[n-1]; range[n-1] = range[r]; range[r] = tmp; // swap selected value to the end
		n--; // values swapped to the end won't get picked again
		}
		
	qsort( s->a, p, sizeof(int), cmpfunc);
	return (s);
	}
	
	SDR * 
	build_random_sdr(int n, int p)
	{
	static int *range = 0;
	SDR * s = sdr_new(n);
	s->p = p;
	
	srand((unsigned int)time(NULL));				
	range = malloc(n*sizeof(int));					// integers 0 to n-1
	for (int i = 0; i < n; i++) range[i] = i;
		
	for (int k = 0; k < p; k++) // random selection of p integers in the range 0 to n-1
		{
		int r = rand() % n;
		s->a[k] = range[r];
		int tmp = range[n-1]; range[n-1] = range[r]; range[r] = tmp; // swap selected value to the end
		n--; // values swapped to the end won't get picked again
		}
		
	qsort( s->a, p, sizeof(int), cmpfunc);
	
	return s;
	}

int sdr_distance( SDR*x, SDR*y) // Hamming distance
	{
	int i = 0, j = 0, h = x->p + y->p;
	
	while (i < x->p && j < y->p )
		{
		if (x->a[i] == y->a[j])
			{ h -= 2; i++; j++; }
		else if (x->a[i] < y->a[j]) ++i;
		else ++j;
		}
	
	return h;
	}
	
SDR *sdr_new(int n)
	{
	SDR *s = malloc(sizeof(SDR));
	s->a = malloc(n*sizeof(int));
	s->n = n;
	s->p = 0;
	return s;
	}
	
TriadicMemory *triadicmemory_new(int n, int p)
	{
	TriadicMemory *t = malloc(sizeof(TriadicMemory));
	// limitation: malloc may fail for large n, use virtual memory instead in this case
	
	t->m = (MEMTYPE*) malloc( n*n*n * sizeof(MEMTYPE));
	t->n = n;
	t->p = p;
	
	MEMTYPE *cube = t->m;
	
	for (int i = 0; i < n*n*n; i++)  *(cube ++) = 0;
	
	return t;
	}
	
void triadicmemory_write (TriadicMemory *T, SDR *x, SDR *y, SDR *z)
	{
	int N = T->n;
	
	for( int i = 0; i < x->p; i++) for( int j = 0; j < y->p; j++) for( int k = 0; k < z->p; k++)
		++ *( T->m + N*N*x->a[i] + N*y->a[j] + z->a[k] ); // counter overflow is unlikely
	}
	
void triadicmemory_delete (TriadicMemory *T, SDR *x, SDR *y, SDR *z)
	{
	int N = T->n;
	
	for( int i = 0; i < x->p; i++) for( int j = 0; j < y->p; j++) for( int k = 0; k < z->p; k++)
		if (*( T->m + N*N*x->a[i] + N*y->a[j] + z->a[k] ) > 0) // checking for counter underflow
			-- *( T->m + N*N*x->a[i] + N*y->a[j] + z->a[k] );
	}
	


SDR* triadicmemory_read_x (TriadicMemory *T, SDR *x, SDR *y, SDR *z)
	{
	int N = T->n;
	
	for( int i = 0; i < N; i++ ) x->a[i] = 0;
		for( int j = 0; j < y->p; j++)  for( int k = 0; k < z->p; k++) for( int i = 0; i < N; i++)
			x->a[i] += *( T->m + N*N*i + N*y->a[j] + z->a[k]);
						
	return( binarize(x, T->p) );
	}

SDR* triadicmemory_read_y (TriadicMemory *T, SDR *x, SDR *y, SDR *z)
	{
	int N = T->n;
	
	for( int j = 0; j < N; j++ ) y->a[j] = 0;
		for( int i = 0; i < x->p; i++) for( int j = 0; j < N; j++) for( int k = 0; k < z->p; k++)
			y->a[j] += *( T->m + N*N*x->a[i] + N*j + z->a[k]);
						
	return( binarize(y, T->p) );
	}

SDR* triadicmemory_read_z (TriadicMemory *T, SDR *x, SDR *y, SDR *z)
	{
	int N = T->n;
	
	for( int k = 0; k < N; k++ ) z->a[k] = 0;
	for( int i = 0; i < x->p; i++) for( int j = 0; j < y->p; j++) for( int k = 0; k < N; k++)
		z->a[k] += *( T->m + N*N*x->a[i] + N*y->a[j] + k);
						
	return( binarize(z, T->p));
	}
	
	
#ifdef TRIADICMEMORY_COMMANDLINE
	
#define SEPARATOR ','
#define QUERY '_'

#define VERSIONMAJOR 1
#define VERSIONMINOR 2


static char* parse (char *buf, SDR *s)
	{
	int *i;
	s->p = 0;
	
	while ( *buf != 0 && *buf != SEPARATOR && *buf != '}')
		{
		while (isspace(*buf)) buf++;
		if (! isdigit(*buf)) break;
		
		i = s->a + s->p;
		sscanf( buf, "%d", i);
		
		if ( (*i)-- > s->n || *i < 0 )
			{
			printf("position out of range: %s\n", buf);
			exit(2);
			}
		s->p ++;
		
		while (isdigit(*buf)) buf++;
		while (isspace(*buf)) buf++;
		}
		
	if (*buf == QUERY && s->p == 0)  { s->p = -1; buf++; while (isspace(*buf)) buf++;}
	
	if (*buf == SEPARATOR) buf++;
	
	return buf;
	}

static void sdr_print(SDR *s)
	{
	for (int r = 0; r < s->p; r++)
		{
		printf("%d", s->a[r] + 1);
		if (r < s->p -1) printf(" ");
		}
	printf("\n"); fflush(stdout);
	}


int test_triadic(int N, int P)
{
	int test_count = 10000;
	int i,errors = 0;
   	clock_t start,end;

	printf("Testing for N=%d and P=%d for %d SDR insertions and queryings\n",N,P,test_count);

	TriadicMemory *T = triadicmemory_new(N, P);

	SDR** xvec = (SDR**) malloc( test_count * sizeof(SDR*));
	SDR** yvec = (SDR**) malloc( test_count * sizeof(SDR*));
	SDR** zvec = (SDR**) malloc( test_count * sizeof(SDR*));
	SDR** qvec = (SDR**) malloc( test_count * sizeof(SDR*)); // query holding vector

	for (i=0;i<test_count;i++)
	{
		xvec[i] = build_random_sdr( N, P);
		yvec[i] = build_random_sdr( N, P);
		zvec[i] = build_random_sdr( N, P);
		qvec[i] = build_random_sdr( N, P);
	}
	
	start = clock() ;
	for (i=0;i<test_count;i++)
	{
		SDR * x = xvec[i];
		SDR * y = yvec[i];
		SDR * z = zvec[i];
		triadicmemory_write(T,x,y,z);
	}
	end = clock() ;
	double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
	printf("SDR Insertion time: %f - %.f per sec\n", elapsed_time, test_count/elapsed_time);

	// Query Z!
	start = clock() ;
	for (i=0;i<test_count;i++)
	{
		SDR * x = xvec[i];
		SDR * y = yvec[i];
		SDR * q = qvec[i];
		triadicmemory_read_z (T, x, y, q);
	}
	end = clock() ;
	elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
	printf("Z SDR query time: %f - %.f per sec: ", elapsed_time, test_count/elapsed_time);

	start = clock() ;
	for (i=0;i<test_count;i++)
	{
		for (int j=0;j<P;j++)
		{
			if (qvec[i]->a[j] != zvec[i]->a[j])
			{
				printf("error: at cycle %d: %d queried bit is %d vs original %d\n", i, j,qvec[i]->a[j], zvec[i]->a[j]);
				errors++;
			}
		}
	}
	end = clock() ;
	elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
	printf("%d errors. (check time: %f - %.f per sec)\n", errors, elapsed_time, test_count/elapsed_time);


	// Query Y!
	start = clock() ;
	for (i=0;i<test_count;i++)
	{
		SDR * x = xvec[i];
		SDR * q = qvec[i];
		SDR * z = zvec[i];
		triadicmemory_read_y (T, x, q, z);
	}
	end = clock() ;
	elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
	printf("Y SDR query time: %f - %.f per sec: ", elapsed_time, test_count/elapsed_time);

	start = clock() ;
	for (i=0;i<test_count;i++)
	{
		for (int j=0;j<P;j++)
		{
			if (qvec[i]->a[j] != yvec[i]->a[j])
			{
				printf("error: at cycle %d: %d queried bit is %d vs original %d\n", i, j,qvec[i]->a[j], yvec[i]->a[j]);
				errors++;
			}
		}
	}
	end = clock() ;
	elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
	printf("%d errors. (check time: %f - %.f per sec)\n", errors, elapsed_time, test_count/elapsed_time);


// query X
	// Query Z!
	start = clock() ;
	for (i=0;i<test_count;i++)
	{
		SDR * q = qvec[i];
		SDR * y = yvec[i];
		SDR * z = zvec[i];
		triadicmemory_read_x (T, q, y, z);
	}
	end = clock() ;
	elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
	printf("X SDR query time: %f - %.f per sec: ", elapsed_time, test_count/elapsed_time);

	start = clock() ;
	for (i=0;i<test_count;i++)
	{
		for (int j=0;j<P;j++)
		{
			if (qvec[i]->a[j] != xvec[i]->a[j])
			{
				printf("error: at cycle %d: %d queried bit is %d vs original %d\n", i, j,qvec[i]->a[j], xvec[i]->a[j]);
				errors++;
			}
		}
	}
	end = clock() ;
	elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
	printf("%d errors. (check time: %f - %.f per sec)\n", errors, elapsed_time, test_count/elapsed_time);

	return errors;
}


int main(int argc, char *argv[])
	{
	char *buf, inputline[10000];
	
	if (argc != 3)
		{
		printf("usage: triadicmemory n p\n");
		printf("n is the hypervector dimension, typically 1000\n");
		printf("p is the target sparse population, typically 10\n");

		}
        
	int N, P;  // SDR dimension and target sparse population, received from command line

    sscanf( argv[1], "%d", &N);
    sscanf( argv[2], "%d", &P);

	if (argc == 4 && !strcmp(argv[3], "test"))
	{
		test_triadic(1000, 10);
		exit(0);
	}
   	TriadicMemory *T = triadicmemory_new(N, P);
    	
	SDR *x = sdr_new(N);
	SDR *y = sdr_new(N);
	SDR *z = sdr_new(N);
	
	while (	fgets(inputline, sizeof(inputline), stdin) != NULL)
		{
		if (! strcmp(inputline, "quit\n"))
			exit(0);
		
		else if (! strcmp(inputline, "random\n"))
			sdr_print(sdr_random(x, P));

		else if ( strcmp(inputline, "version\n") == 0)
			printf("%d.%d\n", VERSIONMAJOR, VERSIONMINOR);
			
		else // parse input of the form { 1 2 3, 4 5 6, 7 8 9 }
			{
			int delete = 0;
			buf = inputline;
			
			if (*buf == '-')
				{ delete = 1; ++buf; }
		
			if (*buf != '{')
				{ printf("expecting '{', found %s\n ", inputline); exit(4); }
		
			buf = parse(parse(parse(buf+1, x), y), z);
		
			if( *buf != '}')
				{ printf("expecting '}', found %s\n ", inputline); exit(4); }
		
			if ( x->p >= 0 && y->p >= 0 && z->p >= 0) // write or delete x, y, z
				{
				if (delete == 0) // write
					triadicmemory_write  (T, x, y, z);
				else // delete
					triadicmemory_delete (T, x, y, z);
				}

			else if ( x->p >= 0 && y->p >= 0 && z->p == -1) // read z
				sdr_print( triadicmemory_read_z (T, x, y, z));
				
			else if ( x->p >= 0 && y->p == -1 && z->p >= 0) // read y
				sdr_print( triadicmemory_read_y (T, x, y, z));

			else if ( x->p == -1 && y->p >= 0 && z->p >= 0) // read x
				sdr_print( triadicmemory_read_x (T, x, y, z));

			else
				{ printf("invalid input\n"); exit(3); }
			}
		}
			
	return 0;
	}
	
#endif // TRIADICMEMORY_COMMANDLINE
