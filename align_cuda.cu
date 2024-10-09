/*
 * Exact genetic sequence alignment
 * (Using brute force)
 * CUDA
 */
#include <complex.h>
#include <stddef.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<sys/time.h>

/* Headers for the CUDA assignment versions */
#include<cuda.h>

/* Example of macros for error checking in CUDA */
#define CUDA_CHECK_FUNCTION( call )	{ cudaError_t check = call; if ( check != cudaSuccess ) fprintf(stderr, "CUDA Error in line: %d, %s\n", __LINE__, cudaGetErrorString(check) ); }
#define CUDA_CHECK_KERNEL( )	{ cudaError_t check = cudaGetLastError(); if ( check != cudaSuccess ) fprintf(stderr, "CUDA Kernel Error in line: %d, %s\n", __LINE__, cudaGetErrorString(check) ); }

/* Arbitrary value to indicate that no matches are found */
#define	NOT_FOUND	-1

/* Arbitrary value to restrict the checksums period */
#define CHECKSUM_MAX	65535


/* 
 * Utils: Function to get wall time
 */
double cp_Wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.0e-6 * tv.tv_usec;
}

/*
 * Utils: Random generator
 */
#include "rng.c"

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 * DO NOT USE OpenMP IN YOUR CODE
 *
 */
/* ADD KERNELS AND OTHER FUNCTIONS HERE */

//5-- Hacer la funcion kernel para generar la secuencia con el nuevo rng_next()
__global__ void new_generate_rng_sequence(rng_t random, float prob_G, float prob_C, float prob_A, char *seq, unsigned long length){
		// CADA HILO VA A LA POSICION DE SU ID GLOBAL
		unsigned long  tid_block = (threadIdx.x);
		unsigned block_size = (blockDim.x);
		unsigned long tid_global = (blockIdx.x * block_size + tid_block);

		//Si su tid esta fuera de rango entonces salirse y no hacer nada
		if( tid_global >= length) return;

		// Dar tantos saltos en como la posicion que te toca
		rng_skip(&random, tid_global);
		float prob = rng_next(&random);

		if( prob < prob_G) seq[tid_global] = 'G';
		else if ( prob < prob_C) seq[tid_global] = 'C';
		else if ( prob < prob_A) seq[tid_global] = 'A';
		else seq[tid_global] = 'T';
}

__global__ void search_pat_by_gpu(char *sequence, unsigned long seq_length, int pat_number, unsigned long *pat_length, char **pattern, unsigned long *pat_found, int *seq_matches){
		unsigned long  tid_block = (threadIdx.x);
		unsigned block_size = ( blockDim.x);
		unsigned long tid_global = ( blockIdx.x * block_size + tid_block);

		int pat;
		unsigned long lind;
		//Go through all the patterns
		if( tid_global < seq_length) {
			for(pat=0; pat<pat_number; pat++){

					//Check it dosent get more far than the sequence length
					if(tid_global > seq_length - pat_length[pat]) continue;

					//Search for the current pattern that belongs to this thread
					for(lind = 0; lind < pat_length[pat]; lind++){
							// The pattern doesnt match
							if(sequence[lind + tid_global] != pattern[pat][lind]) break;
					}
					//Check if the loop has ended with a match of the pattern
					if( lind == pat_length[pat] ){
							//put in pat found the min position found of the pattern
							//BECAUSE notfound is -1 and in unsigned is the biggest number possible
							atomicMin((unsigned long long *) &pat_found[pat], tid_global);
					}
			}
		}
		//Wait to all threads to start to cound seq_matches
		__syncthreads();
		if(pat_found[tid_global] == (unsigned long) NOT_FOUND) return;
		if(tid_global >= pat_number) return;

		//Now increment_matches
		for(lind=0; lind<pat_length[tid_global]; lind++){
				if(seq_matches[ pat_found[tid_global] + lind ] == (unsigned long) NOT_FOUND)
						seq_matches[pat_found[tid_global] + lind] = 0;
				else
						atomicAdd(&seq_matches[pat_found[tid_global] + lind], 1);
		}
}

__global__ void count_number_of_pat(int *return_pat_matches, unsigned long *pat_found, int pat_found_length, unsigned long *return_checksum_patfound, int seq_length, int *seq_matches, unsigned long* return_checksum_matches){
		unsigned long  tid_block = (threadIdx.x);
		unsigned block_size = ( blockDim.x);
		unsigned long tid_global = ( blockIdx.x * block_size + tid_block);

		if( tid_global < pat_found_length){
		
			if(pat_found[tid_global] != (unsigned long) NOT_FOUND){
				   	atomicAdd(return_pat_matches, 1);
					atomicAdd((unsigned long long*)return_checksum_patfound, pat_found[tid_global]);
			}
		}

		if(tid_global >= seq_length) return;
		if(seq_matches[tid_global] != (unsigned long) NOT_FOUND)
					atomicAdd((unsigned long long*)return_checksum_matches, seq_matches[tid_global]);

}


/*
 * Function: Increment the number of pattern matches on the sequence positions
 * 	This function can be changed and/or optimized by the students
 */
void increment_matches( int pat, unsigned long *pat_found, unsigned long *pat_length, int *seq_matches ) {
	unsigned long ind;	
	for( ind=0; ind<pat_length[pat]; ind++) {
		if ( seq_matches[ pat_found[pat] + ind ] == NOT_FOUND )
			seq_matches[ pat_found[pat] + ind ] = 0;
		else
			seq_matches[ pat_found[pat] + ind ] ++;
	}
}
/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

/*
 * Function: Allocate new patttern
 */
char *pattern_allocate( rng_t *random, unsigned long pat_rng_length_mean, unsigned long pat_rng_length_dev, unsigned long seq_length, unsigned long *new_length ) {

	/* Random length */
	unsigned long length = (unsigned long)rng_next_normal( random, (double)pat_rng_length_mean, (double)pat_rng_length_dev );
	if ( length > seq_length ) length = seq_length;
	if ( length <= 0 ) length = 1;

	/* Allocate pattern */
	char *pattern = (char *)malloc( sizeof(char) * length );
	if ( pattern == NULL ) {
		fprintf(stderr,"\n-- Error allocating a pattern of size: %lu\n", length );
		exit( EXIT_FAILURE );
	}

	/* Return results */
	*new_length = length;
	return pattern;
}

/*
 * Function: Fill random sequence or pattern
 */
void generate_rng_sequence( rng_t *random, float prob_G, float prob_C, float prob_A, char *seq, unsigned long length) {
	unsigned long ind; 
	for( ind=0; ind<length; ind++ ) {
		double prob = rng_next( random );
		if( prob < prob_G ) seq[ind] = 'G';
		else if( prob < prob_C ) seq[ind] = 'C';
		else if( prob < prob_A ) seq[ind] = 'A';
		else seq[ind] = 'T';
	}
}

/*
 * Function: Copy a sample of the sequence
 */
void copy_sample_sequence( rng_t *random, char *sequence, unsigned long seq_length, unsigned long pat_samp_loc_mean, unsigned long pat_samp_loc_dev, char *pattern, unsigned long length) {
	/* Choose location */
	unsigned long  location = (unsigned long)rng_next_normal( random, (double)pat_samp_loc_mean, (double)pat_samp_loc_dev );
	if ( location > seq_length - length ) location = seq_length - length;
	if ( location <= 0 ) location = 0;

	/* Copy sample */
	unsigned long ind; 
	for( ind=0; ind<length; ind++ )
		pattern[ind] = sequence[ind+location];
}

/*
 * Function: Regenerate a sample of the sequence
 */
void generate_sample_sequence( rng_t *random, rng_t random_seq, float prob_G, float prob_C, float prob_A, unsigned long seq_length, unsigned long pat_samp_loc_mean, unsigned long pat_samp_loc_dev, char *pattern, unsigned long length ) {
	/* Choose location */
	unsigned long  location = (unsigned long)rng_next_normal( random, (double)pat_samp_loc_mean, (double)pat_samp_loc_dev );
	if ( location > seq_length - length ) location = seq_length - length;
	if ( location <= 0 ) location = 0;

	/* Regenerate sample */
	rng_t local_random = random_seq;
	rng_skip( &local_random, location );
	generate_rng_sequence( &local_random, prob_G, prob_C, prob_A, pattern, length);
}


/*
 * Function: Print usage line in stderr
 */
void show_usage( char *program_name ) {
	fprintf(stderr,"Usage: %s ", program_name );
	fprintf(stderr,"<seq_length> <prob_G> <prob_C> <prob_A> <pat_rng_num> <pat_rng_length_mean> <pat_rng_length_dev> <pat_samples_num> <pat_samp_length_mean> <pat_samp_length_dev> <pat_samp_loc_mean> <pat_samp_loc_dev> <pat_samp_mix:B[efore]|A[fter]|M[ixed]> <long_seed>\n");
	fprintf(stderr,"\n");
}



/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[]) {
	/* 0. Default output and error without buffering, forces to write immediately */
	setbuf(stdout, NULL);
	setbuf(stderr, NULL);

	/* 1. Read scenary arguments */
	/* 1.1. Check minimum number of arguments */
	if (argc < 15) {
		fprintf(stderr, "\n-- Error: Not enough arguments when reading configuration from the command line\n\n");
		show_usage( argv[0] );
		exit( EXIT_FAILURE );
	}

	/* 1.2. Read argument values */
	unsigned long seq_length = atol( argv[1] );
	float prob_G = atof( argv[2] );
	float prob_C = atof( argv[3] );
	float prob_A = atof( argv[4] );
	if ( prob_G + prob_C + prob_A > 1 ) {
		fprintf(stderr, "\n-- Error: The sum of G,C,A,T nucleotid probabilities cannot be higher than 1\n\n");
		show_usage( argv[0] );
		exit( EXIT_FAILURE );
	}
	prob_C += prob_G;
	prob_A += prob_C;

	int pat_rng_num = atoi( argv[5] );
	unsigned long pat_rng_length_mean = atol( argv[6] );
	unsigned long pat_rng_length_dev = atol( argv[7] );
	
	int pat_samp_num = atoi( argv[8] );
	unsigned long pat_samp_length_mean = atol( argv[9] );
	unsigned long pat_samp_length_dev = atol( argv[10] );
	unsigned long pat_samp_loc_mean = atol( argv[11] );
	unsigned long pat_samp_loc_dev = atol( argv[12] );

	char pat_samp_mix = argv[13][0];
	if ( pat_samp_mix != 'B' && pat_samp_mix != 'A' && pat_samp_mix != 'M' ) {
		fprintf(stderr, "\n-- Error: Incorrect first character of pat_samp_mix: %c\n\n", pat_samp_mix);
		show_usage( argv[0] );
		exit( EXIT_FAILURE );
	}

	unsigned long seed = atol( argv[14] );

#ifdef DEBUG
	/* DEBUG: Print arguments */
	printf("\nArguments: seq_length=%lu\n", seq_length );
	printf("Arguments: Accumulated probabilitiy G=%f, C=%f, A=%f, T=1\n", prob_G, prob_C, prob_A );
	printf("Arguments: Random patterns number=%d, length_mean=%lu, length_dev=%lu\n", pat_rng_num, pat_rng_length_mean, pat_rng_length_dev );
	printf("Arguments: Sample patterns number=%d, length_mean=%lu, length_dev=%lu, loc_mean=%lu, loc_dev=%lu\n", pat_samp_num, pat_samp_length_mean, pat_samp_length_dev, pat_samp_loc_mean, pat_samp_loc_dev );
	printf("Arguments: Type of mix: %c, Random seed: %lu\n", pat_samp_mix, seed );
	printf("\n");
#endif // DEBUG

        CUDA_CHECK_FUNCTION( cudaSetDevice(0) );

	/* 2. Initialize data structures */
	/* 2.1. Skip allocate and fill sequence */
	rng_t random = rng_new( seed );
	rng_skip( &random, seq_length );

	/* 2.2. Allocate and fill patterns */
	/* 2.2.1 Allocate main structures */
	int pat_number = pat_rng_num + pat_samp_num;
	unsigned long *pat_length = (unsigned long *)malloc( sizeof(unsigned long) * pat_number );
	char **pattern = (char **)malloc( sizeof(char*) * pat_number );
	if ( pattern == NULL || pat_length == NULL ) {
		fprintf(stderr,"\n-- Error allocating the basic patterns structures for size: %d\n", pat_number );
		exit( EXIT_FAILURE );
	}

	/* 2.2.2 Allocate and initialize ancillary structure for pattern types */
	int ind;
	unsigned long lind;
	#define PAT_TYPE_NONE	0
	#define PAT_TYPE_RNG	1
	#define PAT_TYPE_SAMP	2
	char *pat_type = (char *)malloc( sizeof(char) * pat_number );
	if ( pat_type == NULL ) {
		fprintf(stderr,"\n-- Error allocating ancillary structure for pattern of size: %d\n", pat_number );
		exit( EXIT_FAILURE );
	}
	for( ind=0; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_NONE;

	/* 2.2.3 Fill up pattern types using the chosen mode */
	switch( pat_samp_mix ) {
	case 'A':
		for( ind=0; ind<pat_rng_num; ind++ ) pat_type[ind] = PAT_TYPE_RNG;
		for( ; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_SAMP;
		break;
	case 'B':
		for( ind=0; ind<pat_samp_num; ind++ ) pat_type[ind] = PAT_TYPE_SAMP;
		for( ; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_RNG;
		break;
	default:
		if ( pat_rng_num == 0 ) {
			for( ind=0; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_SAMP;
		}
		else if ( pat_samp_num == 0 ) {
			for( ind=0; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_RNG;
		}
		else if ( pat_rng_num < pat_samp_num ) {
			int interval = pat_number / pat_rng_num;
			for( ind=0; ind<pat_number; ind++ ) 
				if ( (ind+1) % interval == 0 ) pat_type[ind] = PAT_TYPE_RNG;
				else pat_type[ind] = PAT_TYPE_SAMP;
		}
		else {
			int interval = pat_number / pat_samp_num;
			for( ind=0; ind<pat_number; ind++ ) 
				if ( (ind+1) % interval == 0 ) pat_type[ind] = PAT_TYPE_SAMP;
				else pat_type[ind] = PAT_TYPE_RNG;
		}
	}

	/* 2.2.4 Generate the patterns */
	for( ind=0; ind<pat_number; ind++ ) {
		if ( pat_type[ind] == PAT_TYPE_RNG ) {
			pattern[ind] = pattern_allocate( &random, pat_rng_length_mean, pat_rng_length_dev, seq_length, &pat_length[ind] );
			generate_rng_sequence( &random, prob_G, prob_C, prob_A, pattern[ind], pat_length[ind] );
		}
		else if ( pat_type[ind] == PAT_TYPE_SAMP ) {
			pattern[ind] = pattern_allocate( &random, pat_samp_length_mean, pat_samp_length_dev, seq_length, &pat_length[ind] );
#define REGENERATE_SAMPLE_PATTERNS
#ifdef REGENERATE_SAMPLE_PATTERNS
			rng_t random_seq_orig = rng_new( seed );
			generate_sample_sequence( &random, random_seq_orig, prob_G, prob_C, prob_A, seq_length, pat_samp_loc_mean, pat_samp_loc_dev, pattern[ind], pat_length[ind] );
#else
			copy_sample_sequence( &random, sequence, seq_length, pat_samp_loc_mean, pat_samp_loc_dev, pattern[ind], pat_length[ind] );
#endif
		}
		else {
			fprintf(stderr,"\n-- Error internal: Paranoic check! A pattern without type at position %d\n", ind );
			exit( EXIT_FAILURE );
		}
	}
	free( pat_type );

	/* Allocate and move the patterns to the GPU */
	unsigned long *d_pat_length;
	char **d_pattern;
	CUDA_CHECK_FUNCTION( cudaMalloc( &d_pat_length, sizeof(unsigned long) * pat_number ) );
	CUDA_CHECK_FUNCTION( cudaMalloc( &d_pattern, sizeof(char *) * pat_number ) );

	char **d_pattern_in_host = (char **)malloc( sizeof(char*) * pat_number );
	if ( d_pattern_in_host == NULL ) {
		fprintf(stderr,"\n-- Error allocating the patterns structures replicated in the host for size: %d\n", pat_number );
		exit( EXIT_FAILURE );
	}
	for( ind=0; ind<pat_number; ind++ ) {
		CUDA_CHECK_FUNCTION( cudaMalloc( &(d_pattern_in_host[ind]), sizeof(char *) * pat_length[ind] ) );
        	CUDA_CHECK_FUNCTION( cudaMemcpy( d_pattern_in_host[ind], pattern[ind], pat_length[ind] * sizeof(char), cudaMemcpyHostToDevice ) );
	}
	CUDA_CHECK_FUNCTION( cudaMemcpy( d_pattern, d_pattern_in_host, pat_number * sizeof(char *), cudaMemcpyHostToDevice ) );

	/* Avoid the usage of arguments to take strategic decisions
	 * In a real case the user only has the patterns and sequence data to analize
	 */
	argc = 0;
	argv = NULL;
	pat_rng_num = 0;
	pat_rng_length_mean = 0;
	pat_rng_length_dev = 0;
	pat_samp_num = 0;
	pat_samp_length_mean = 0;
	pat_samp_length_dev = 0;
	pat_samp_loc_mean = 0;
	pat_samp_loc_dev = 0;
	pat_samp_mix = '0';

	/* 2.3. Other result data and structures */
	int pat_matches = 0;

	/* 2.3.1. Other results related to patterns */
	unsigned long *pat_found;
	pat_found = (unsigned long *)malloc( sizeof(unsigned long) * pat_number );
	if ( pat_found == NULL ) {
		fprintf(stderr,"\n-- Error allocating aux pattern structure for size: %d\n", pat_number );
		exit( EXIT_FAILURE );
	}
	
	/* 3. Start global timer */
    CUDA_CHECK_FUNCTION( cudaDeviceSynchronize() );
	double ttotal = cp_Wtime();

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 * DO NOT USE OpenMP IN YOUR CODE
 *
 */

	/* 2.1. Allocate and fill sequence */

	//TODO quitar cuando ya no se necesite
	char *sequence ;

	//2-- Reservar espacio suficiente en la GPU
	char *device_sequence;
	CUDA_CHECK_FUNCTION(cudaMalloc(&device_sequence, sizeof(char) * seq_length););
	
	random = rng_new( seed );
	//generate_rng_sequence( &random, prob_G, prob_C, prob_A, sequence, seq_length);

	//3--Generar la sequencia con una funcion Kernel
	// El numero de bloques y su tamaño ( que cada hilo inicialice un de las posiciones de la secuencia)
	// numero de hilos = tamaño de la secuencia
	// numero de grids = numero de threads / 1024
	int blockSize = 1024; //1024 hilos por bloque
	int nunmBlocks;
	if (seq_length % blockSize == 0 ) nunmBlocks = seq_length / blockSize;
	else nunmBlocks = (seq_length / blockSize ) + 1;

	//Usando la memoria reservada para cuda
	new_generate_rng_sequence<<<nunmBlocks, blockSize>>>(random,  prob_G,  prob_C,  prob_A,  device_sequence, seq_length);

#ifdef DEBUG
	/* DEBUG: Print sequence and patterns */
	printf("-----------------\n");
	printf("Sequence: ");
	for( lind=0; lind<seq_length; lind++ ) 
		printf( "%c", sequence[lind] );
	printf("\n-----------------\n");
	printf("Patterns: %d ( rng: %d, samples: %d )\n", pat_number, pat_rng_num, pat_samp_num );
	int debug_pat;
	for( debug_pat=0; debug_pat<pat_number; debug_pat++ ) {
		printf( "Pat[%d]: ", debug_pat );
		for( lind=0; lind<pat_length[debug_pat]; lind++ ) 
			printf( "%c", pattern[debug_pat][lind] );
		printf("\n");
	}
	printf("-----------------\n\n");
#endif // DEBUG

	/* 2.3.2. Other results related to the main sequence */
	int *seq_matches;
	seq_matches = (int *)malloc( sizeof(int) * seq_length );
	if ( seq_matches == NULL ) {
		fprintf(stderr,"\n-- Error allocating aux sequence structures for size: %lu\n", seq_length );
		exit( EXIT_FAILURE );
	}

	/* 4. Initialize ancillary structures */
	for( ind=0; ind<pat_number; ind++) {
		pat_found[ind] = (unsigned long)NOT_FOUND;
	}

	for( lind=0; lind<seq_length; lind++) {
		seq_matches[lind] = 0;
	}

	/* 5. Search for each pattern */

	/* copy to gpu mem pat_found and pat_matches */
	//unsigned long *device_pat_matches;
	int * device_pat_matches;
	int *device_seq_matches;
	unsigned long *device_pat_found;
	unsigned long *device_checksum_found;
	unsigned long *device_checksum_matches;
	unsigned long checksum_found = 0;
	unsigned long checksum_matches =0;

	CUDA_CHECK_FUNCTION(cudaMalloc( &device_seq_matches, sizeof(int)*seq_length));
	CUDA_CHECK_FUNCTION(cudaMemcpy(device_seq_matches, seq_matches, sizeof(int)*seq_length, cudaMemcpyHostToDevice));
	
	CUDA_CHECK_FUNCTION(cudaMalloc( &device_pat_matches, sizeof(int)));
	CUDA_CHECK_FUNCTION(cudaMemcpy(device_pat_matches, &pat_matches, sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CHECK_FUNCTION(cudaMalloc(&device_pat_found, sizeof(unsigned long)*pat_number));
	CUDA_CHECK_FUNCTION(cudaMemcpy(device_pat_found ,pat_found, sizeof(unsigned long)*pat_number,cudaMemcpyHostToDevice));

	CUDA_CHECK_FUNCTION(cudaMalloc(&device_checksum_found, sizeof(unsigned long)));
	CUDA_CHECK_FUNCTION(cudaMemcpy(device_checksum_found, &checksum_found, sizeof(unsigned long), cudaMemcpyHostToDevice));

	CUDA_CHECK_FUNCTION(cudaMalloc(&device_checksum_matches, sizeof(unsigned long)));
	CUDA_CHECK_FUNCTION(cudaMemcpy(device_checksum_matches, &checksum_matches, sizeof(unsigned long), cudaMemcpyHostToDevice));

	// copy the length of each pattern
	CUDA_CHECK_FUNCTION(cudaMemcpy(d_pat_length, pat_length, sizeof(unsigned long)*pat_number, cudaMemcpyHostToDevice));	

	search_pat_by_gpu<<<nunmBlocks, blockSize>>>(device_sequence, seq_length, pat_number, d_pat_length, d_pattern, device_pat_found, device_seq_matches);
	CUDA_CHECK_KERNEL();

	count_number_of_pat<<<nunmBlocks, blockSize>>>(device_pat_matches, device_pat_found, pat_number, device_checksum_found, seq_length, device_seq_matches, device_checksum_matches);
	CUDA_CHECK_KERNEL();

	//copy the data to the host
	CUDA_CHECK_FUNCTION(cudaMemcpy(&pat_matches, device_pat_matches, sizeof(int), cudaMemcpyDeviceToHost));
	CUDA_CHECK_FUNCTION(cudaMemcpy(&checksum_found, device_checksum_found, sizeof(unsigned long) ,cudaMemcpyDeviceToHost));
	CUDA_CHECK_FUNCTION(cudaMemcpy(&checksum_matches, device_checksum_matches, sizeof(unsigned long) ,cudaMemcpyDeviceToHost));
	
	/* 7. Check sums */
	// calc now the rest operation
	checksum_found = checksum_found % CHECKSUM_MAX;
	checksum_matches= checksum_matches % CHECKSUM_MAX;

#ifdef DEBUG
	/* DEBUG: Write results */
	printf("-----------------\n");
	printf("Found start:");
	for( debug_pat=0; debug_pat<pat_number; debug_pat++ ) {
		printf( " %lu", pat_found[debug_pat] );
	}
	printf("\n");
	printf("-----------------\n");
	printf("Matches:");
	for( lind=0; lind<seq_length; lind++ ) 
		printf( " %d", seq_matches[lind] );
	printf("\n");
	printf("-----------------\n");
#endif // DEBUG

	/* Free local resources */	
	free( sequence );
	free( seq_matches );

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

	/* 8. Stop global timer */
    CUDA_CHECK_FUNCTION( cudaDeviceSynchronize() );
	ttotal = cp_Wtime() - ttotal;

	/* 9. Output for leaderboard */
	printf("\n");
	/* 9.1. Total computation time */
	printf("Time: %lf\n", ttotal );
	
	/* 9.2. Results: Statistics */
	printf("Result: %d, %lu, %lu\n\n", 
			pat_matches,
			checksum_found,
			checksum_matches );
		
	/* 10. Free resources */	
	int i;
	for( i=0; i<pat_number; i++ ) free( pattern[i] );
	free( pattern );
	free( pat_length );
	free( pat_found );

	/* 11. End */
	return 0;
}
