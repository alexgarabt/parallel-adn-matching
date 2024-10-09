/*
 * Exact genetic sequence alignment
 * (Using brute force)
 *
 * MPI version
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<sys/time.h>
#include<mpi.h>


/* Arbitrary value to indicate that no matches are found */
#define	NOT_FOUND	-1

/* Arbitrary value to restrict the checksums period */
#define CHECKSUM_MAX	65535

/* TAGS for comunnications */
#define TAG_NO_COMPLETE 69
#define END_SIG -2

/* Macro for checking mpi funcs and print if there is a error */
#define MPI_CHECK_FUNCTION( call ) {int ierr = call; int resultlen; if (ierr != MPI_SUCCESS){ char error_msg[MPI_MAX_ERROR_STRING]; MPI_Error_string(ierr, error_msg, &resultlen); fprintf(stderr, "MPI fallo en linea (%d):\n%s\n\n", __LINE__, error_msg); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); fflush(stdout);}}

/*
 * Function: Increment the number of pattern matches on the sequence positions
 * 	This function can be changed and/or optimized by the students
 */
/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */
void increment_matches( int pat, unsigned long *pat_found, unsigned long *pat_length, int *seq_matches ) {
	int ind;	
	for( ind=0; ind<pat_length[pat]; ind++) {
		if ( seq_matches[ pat_found[pat] + ind ] == (unsigned long)NOT_FOUND )
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
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
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
	/* 1.0. Init MPI before processing arguments */
	MPI_Init( &argc, &argv );
	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	/* 1.1. Check minimum number of arguments */
	if (argc < 15) {
		fprintf(stderr, "\n-- Error: Not enough arguments when reading configuration from the command line\n\n");
		show_usage( argv[0] );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* 1.2. Read argument values */
	unsigned long seq_length = atol( argv[1] );
	float prob_G = atof( argv[2] );
	float prob_C = atof( argv[3] );
	float prob_A = atof( argv[4] );
	if ( prob_G + prob_C + prob_A > 1 ) {
		fprintf(stderr, "\n-- Error: The sum of G,C,A,T nucleotid probabilities cannot be higher than 1\n\n");
		show_usage( argv[0] );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
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
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	unsigned long seed = atol( argv[14] );

#ifdef DEBUG
	/* DEBUG: Print arguments */
	if ( rank == 0 ) {
		printf("\nArguments: seq_length=%lu\n", seq_length );
		printf("Arguments: Accumulated probabilitiy G=%f, C=%f, A=%f, T=1\n", prob_G, prob_C, prob_A );
		printf("Arguments: Random patterns number=%d, length_mean=%lu, length_dev=%lu\n", pat_rng_num, pat_rng_length_mean, pat_rng_length_dev );
		printf("Arguments: Sample patterns number=%d, length_mean=%lu, length_dev=%lu, loc_mean=%lu, loc_dev=%lu\n", pat_samp_num, pat_samp_length_mean, pat_samp_length_dev, pat_samp_loc_mean, pat_samp_loc_dev );
		printf("Arguments: Type of mix: %c, Random seed: %lu\n", pat_samp_mix, seed );
		printf("\n");
	}
#endif // DEBUG

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
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* 2.2.2 Allocate and initialize ancillary structure for pattern types */
	unsigned long ind;
	#define PAT_TYPE_NONE	0
	#define PAT_TYPE_RNG	1
	#define PAT_TYPE_SAMP	2
	char *pat_type = (char *)malloc( sizeof(char) * pat_number );
	if ( pat_type == NULL ) {
		fprintf(stderr,"\n-- Error allocating ancillary structure for pattern of size: %d\n", pat_number );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
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
			fprintf(stderr,"\n-- Error internal: Paranoic check! A pattern without type at position %lu\n", ind );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}
	}
	free( pat_type );

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
	pat_samp_mix = '0';

	/* 2.3.1. Other results related to patterns */
	unsigned long *pat_found;
	pat_found = (unsigned long *)malloc( sizeof(unsigned long) * pat_number );
	if ( pat_found == NULL ) {
		fprintf(stderr,"\n-- Error allocating aux pattern structure for size: %d\n", pat_number );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	
	/* 3. Start global timer */
	MPI_Barrier( MPI_COMM_WORLD );
	double ttotal = cp_Wtime();

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */
	/* 2.1. Allocate and fill sequence */

	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	unsigned long my_size = seq_length/nprocs;
	unsigned long resto = seq_length%nprocs;
	unsigned long my_begin;// = my_size*rank;

	if(rank < resto) my_size++;

	if (rank < resto)
            my_begin = rank * my_size;
        else
            my_begin = rank * my_size + resto;


	char *sequence = (char *)malloc( sizeof(char) * my_size );
	if ( sequence == NULL ) {
		fprintf(stderr,"\n-- Error allocating the sequence for size: %lu\n", my_size );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	random = rng_new( seed );
	rng_t random_local = random;
	rng_skip(&random_local, my_begin);	
	generate_rng_sequence( &random_local, prob_G, prob_C, prob_A, sequence, my_size);

#ifdef DEBUG
	/* DEBUG: Print sequence and patterns */
	printf("-----------------\n");
	printf("[%d]Sequence: ",rank);
	for( ind=0; ind<my_size; ind++ ) 
		printf( "%c", sequence[ind] );
	printf("\n-----------------\n");
	printf("Patterns: %d ( rng: %d, samples: %d )\n", pat_number, pat_rng_num, pat_samp_num );
	unsigned long debug_pat;
	for( debug_pat=0; debug_pat<pat_number; debug_pat++ ) {
		printf( "Pat[%lu]: ", debug_pat );
		for( ind=0; ind<pat_length[debug_pat]; ind++ ) 
			printf( "%c", pattern[debug_pat][ind] );
		printf("\n");
	}
	printf("-----------------\n\n");
#endif // DEBUG

	/* 2.3.2. Other results related to the main sequence */
	

	/* 4. Initialize ancillary structures */
	for( ind=0; ind<pat_number; ind++) {
		pat_found[ind] = seq_length;
	}

	/* 5. Search for each pattern */
	//Legacy var
	int *seq_matches;
	int pat_matches = 0;

	unsigned long start;
	unsigned long pat;
	//MPI_Request request;
	MPI_Status status;


	// [0]=my_begin+start | [1]= tamaño ya leido del patorn | [2]= numero de patrón
	unsigned long *envio_array = (unsigned long*) malloc(sizeof(unsigned long)*3);
	if ( envio_array == NULL ) {
			fprintf(stderr,"Error al crear envio_array en proc:%d\n", rank );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	
	unsigned long *recibo_array = (unsigned long*) malloc(sizeof(unsigned long)*3);
	if ( recibo_array == NULL ) {
			fprintf(stderr,"Error al crear recibo_array en proc:%d\n", rank );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* info for rank and my_size distrib
	printf("rank[%d] my_size=%lu, my_begin=%lu\n", rank,  my_size, my_begin);
	fflush(stdout);
	*/

	MPI_CHECK_FUNCTION(MPI_Barrier(MPI_COMM_WORLD));

	for( pat=0; pat < pat_number; pat++ ) {

		/* 5.1. For each posible starting position */
		for( start=0; start < my_size;start++) {

			/* 5.1.1. For each pattern element */
			for( ind=0; ind<pat_length[pat]; ind++) {
				/* Stop this test when different nucleotids are found */
				if( (start + ind >= my_size) && (rank != nprocs-1)){

					envio_array[0]=my_begin + start;	
					envio_array[1]=ind; //Last position checked
					envio_array[2]=pat;

					//Envio el parton para que continue el siguente que le toca
					MPI_CHECK_FUNCTION(MPI_Send(envio_array, 3, MPI_UNSIGNED_LONG, rank+1, TAG_NO_COMPLETE, MPI_COMM_WORLD));
						//printf("\n ENVIA(%d) rank=%d a rank+1=%d pat=%lu, START=%lu, LENGHT=%lu\n", (1), rank, rank+1, envio_array[2], envio_array[0], envio_array[1]);
						//fflush(stdout);
					break;
				}
				else if ( sequence[start + ind] != pattern[pat][ind] ) break;
			}
			/* 5.1.2. Check if the loop ended with a match */
			if ( ind == pat_length[pat] ) {
				pat_found[pat] = my_begin + start;
				//printf("RANK(%d) encuentra patrón=%lu\n",rank, pat);
				break;
			}
		}
	}

	/**
	MPI_CHECK_FUNCTION(MPI_Barrier(MPI_COMM_WORLD));
	printf("==================\n rank=%d LLEGA pat=%lu, sig=%d\n", rank, pat, rank+1);
	fflush(stdout);
	MPI_CHECK_FUNCTION(MPI_Barrier(MPI_COMM_WORLD));
	*/
	
	//NOW PROCESS CHECK FOR THE REAMAING PATTERN THAT ITS PREVIOUS DOSENT HAVE FOUND
	unsigned long start_r, pat_r, length_r;
	while(1 && (rank != 0)){

		//First check if the previous process has finished & have nothing more to do
		MPI_CHECK_FUNCTION(MPI_Recv(recibo_array, 3, MPI_UNSIGNED_LONG, rank-1, TAG_NO_COMPLETE, MPI_COMM_WORLD, &status));

		//else continue because the rank-1 send work to do
		start_r = recibo_array[0];   //Were the patter was start to be detected
		length_r = recibo_array[1];  //Lenght already checked
		pat_r = recibo_array[2];     //Pattern searching

 		//	printf("\n RECIBE rank=%d desde rank-1=%d pat=%lu con TAG=%d, START=%lu, LENGHT=%lu\n" ,rank, status.MPI_SOURCE, pat_r,  status.MPI_TAG, start_r, length_r);
		//	fflush(stdout);
		
		//THEN its end sigal the first so no more thing to do
		if(pat_r == (unsigned long) END_SIG){
		//	printf("\n GET FIN rank=%d desde rank-1=%d \n",  rank, status.MPI_SOURCE);
		//	fflush(stdout);
			break;
		}
		
		//search for the pattern in my space
		//we only check like start = 0 / because of it will be a continuos pattern that start ranks before this
		for( ind=0; ind+length_r<pat_length[pat_r]; ind++) {
				//The pattern gets out of this region so send to the next rank
				if( (ind >= my_size) && (rank != nprocs-1)){
					envio_array[0]= start_r;	
					envio_array[1]=ind+length_r; //Length already checked
					envio_array[2]=pat_r;

					//Send to the next rank
					MPI_CHECK_FUNCTION(MPI_Send(envio_array, 3, MPI_UNSIGNED_LONG, rank+1, TAG_NO_COMPLETE, MPI_COMM_WORLD));
					//printf("\n ENVIA(%d) rank=%d a rank+1=%d pat=%lu, START=%lu, LENGHT=%lu\n", (2), rank, rank+1, envio_array[2], envio_array[0], envio_array[1]);
					//fflush(stdout);
					break;
				}
				else if ( sequence[ind] != pattern[pat_r][ind+length_r] ) break; 
			}
			// else we have a match
			if ( ind+length_r == pat_length[pat_r] ) {
				pat_found[pat_r] = start_r; 
				//printf("RANK(%d) encuentra patrón=%lu\n",rank, pat_r);
			}
	}

	if(rank != nprocs-1){
		envio_array[0] =(unsigned long) END_SIG;
		envio_array[1] =(unsigned long) END_SIG;
		envio_array[2] = (unsigned long)END_SIG;
		MPI_CHECK_FUNCTION(MPI_Send(envio_array, 3, MPI_UNSIGNED_LONG, rank+1, TAG_NO_COMPLETE, MPI_COMM_WORLD));
		//printf("\n ENVIA FIN rank=%d a rank+1=%d \n",  rank, rank+1);
		//fflush(stdout);
	}

	//Wait until all the process finish
	//MPI_CHECK_FUNCTION(MPI_Barrier(MPI_COMM_WORLD));
	/*
	MPI_CHECK_FUNCTION(MPI_Barrier(MPI_COMM_WORLD));
	printf("==================\n rank=%d LLEGA pat=%lu, sig=%d\n", rank, pat, rank+1);
	fflush(stdout);
	*/
		
		
	unsigned long *recive_pat_found;
	//Define a recv buffer for for the master
	if(rank==0){
		recive_pat_found = (unsigned long *)malloc( sizeof(unsigned long) * pat_number );
		if ( recive_pat_found == NULL ) {
			fprintf(stderr,"\n-- Error allocating aux pattern structure for size: %d\n", pat_number );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}
	}
	
	MPI_CHECK_FUNCTION(MPI_Reduce(pat_found,  recive_pat_found, pat_number, MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD));
	
	//debug	 info
	
	
		
	


	/* 7. Check sums */
	unsigned long checksum_matches = 0;
	unsigned long checksum_found = 0;

	// TODO Increment_matches 

	// rank 0 calculate the pat_matches
	if( rank == 0) {
		for(ind=0;ind<pat_number;ind++)
			if(recive_pat_found[ind] != seq_length){
				pat_matches++;
				checksum_found = (checksum_found + recive_pat_found[ind]) % CHECKSUM_MAX;
				checksum_matches = ( checksum_matches + pat_length[ind] ) % CHECKSUM_MAX;
			} 

	}
	
#ifdef DEBUG
	/* DEBUG: Write results */
	printf("-----------------\n");
	printf("[%d]Found start:",rank);
	for( debug_pat=0; debug_pat<pat_number; debug_pat++ ) {
		printf( " %lu", pat_found[debug_pat] );
	}
	printf("\n");
	printf("-----------------\n");
	printf("Matches:");
	/*for( ind=0; ind<my_size; ind++ ) 
		printf( " %d", seq_matches[ind] );*/
	printf("\n");
	printf("-----------------\n");
#endif // DEBUG

	/* Free local resources */	
	free( sequence );
	//free( seq_matches );

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

		/* 8. Stop global time */
	MPI_Barrier( MPI_COMM_WORLD );
	ttotal = cp_Wtime() - ttotal;

	/* 9. Output for leaderboard */
	if ( rank == 0 ) {
		printf("\n");
		/* 9.1. Total computation time */
		printf("Time: %lf\n", ttotal );

		/* 9.2. Results: Statistics */
		printf("Result: %d, %lu, %lu\n\n", 
				pat_matches,
				checksum_found,
				checksum_matches );

		
	}
	/*
	fflush(stdout);
	MPI_CHECK_FUNCTION(MPI_Barrier(MPI_COMM_WORLD));
	unsigned long ind1;
	printf("rank[%d] pat_found:", rank);
	for( ind1=0; ind1<pat_number; ind1++ ) {
		printf( " %lu", pat_found[ind1]);
	}
	printf("\n");
	fflush(stdout);

	MPI_CHECK_FUNCTION(MPI_Barrier(MPI_COMM_WORLD));
	if(rank==0){
	//debug info
		printf("recive_pat_found:");
		for( ind1=0; ind1<pat_number; ind1++ ) {
			printf( " %lu", recive_pat_found[ind1] );
		}
	printf("\n");
	fflush(stdout);
	}
	*/


				
	/* 10. Free resources */	
	int i;
	for( i=0; i<pat_number; i++ ) free( pattern[i] );
	free( pattern );
	free( pat_length );
	free( pat_found );

	/* 11. End */
	MPI_Finalize();
	return 0;
}
