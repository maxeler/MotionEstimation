/*********************************************************************
 * Maxeler Technologies : Motion Estimation                          *
 *                                                                   *
 * Version: 1.0                                                      *
 * Date:    12th August 2014                                         *
 *                                                                   *
 * CPU code source file                                              *
 *                                                                   *
 *********************************************************************/

/*
 * motion_estimation.c
 *
 *  last edit: 13 March 2015
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Maxfiles.h>
#include <time.h>

/*
 * prints the proper usage of the program
 */

void print_cmd_help(){
	printf("Usage: motionEstimation <x> <y>\n");
	printf("       x - number of blocks in width\n");
	printf("       y - number of blocks in height\n");
}



/*
 * Randomly assigns values to source and reference input arrays.
 * ref2d is a duplicate of reference array organized as a 2D array, used in verifying the DFE results
 * src2d is a duplicate of source array organized as a 2D array, used in verifying the DFE results
 */

void create_test_data(int width_in_blocks, int height_in_blocks, int halo_blocks_x, int halo_blocks_y, int **ref2d, int **src2d, int16_t* reference, int16_t* source){
	srand(time(NULL));
	for (int i = 0; i < height_in_blocks + 2*halo_blocks_y; i++){
		for (int j = 0 ; j < width_in_blocks + 2*halo_blocks_x;  j++){
			for (int k = 0; k < ME_blockSize; k++){
				for (int l = 0; l < ME_blockSize; l ++){
					int index = ((i * (width_in_blocks + 2*halo_blocks_x) + j) * ME_blockSize + k) * ME_blockSize + l;
					reference[index] = rand() % 1024;
					ref2d[i*ME_blockSize+k][j*ME_blockSize+l] = reference[index];
					if(j < width_in_blocks && i < height_in_blocks){
						int index_src = ((i * (width_in_blocks) + j) * ME_blockSize + k) * ME_blockSize + l;
						source[index_src] = rand() % 1024;
						src2d[i*ME_blockSize+k][j*ME_blockSize+l] = source[index_src];
					}
				}
			}
		}
	}
}



/*
 * Rounds up the input "numToRound" to a multiple of input "multiple"
 * e.g. roundUp(64,20) returns 80
 */

int roundUp(int numToRound, int multiple)
{
	if(multiple == 0)
	{
		return numToRound;
	}

	int remainder = numToRound % multiple;
	if (remainder == 0)
		return numToRound;
	return numToRound + multiple - remainder;
}



/*
 * Returns the size of each element of the output
 * each element consists of 32 bit motion vectors for each sub-block of the source block
 */

int getOutputTypeSizeInBytes(int block_size){
	int total_bits = 0;

	for(int i = 4; i <= block_size; i*=2){   // iterate through different sizes of sub blocks
		int numBlocks = block_size * block_size / (i*i);
		total_bits += numBlocks * 32;
	}
	return roundUp(total_bits, 128)/8;      // The output size should be a multiple of 128 bits
}



/*
 * Returns the search window coordinates for all the source blocks.
 * The using routine can call get_window_coords and then iterate through the coords arrays on all indexes
 *    to get the top left corner coordinates of all the reference blocks in the search windows
 */

void get_window_coords(int **coords, int width_in_blocks, int height_in_blocks, int halo_blocks_x, int halo_blocks_y){
	int idx = 0;
	for (int row = 0; row < height_in_blocks; row ++){
		for (int col = 0 ; col < width_in_blocks; col++){

			/* find the top left corner of the search block */
			int top = (row+halo_blocks_y) * ME_blockSize - ME_windowY / 2;
			int left = (col+halo_blocks_x) * ME_blockSize - ME_windowX / 2;

			/* iterate through all pixels in the search window */
			for (int y = top; y < top + ME_windowY; y++){
				for(int x = left; x < left + ME_windowX; x++){
					coords[0][idx] = x;
					coords[1][idx] = y;
					idx++;
				}
			}

		}
	}
}


/* copies a block of source 2d array into a 2d buffer with the size of a block
 * The copy takes place based on an index referring to all the blocks organized for iterating through all the blocks
 */

void get_source_block(int index, int width_in_blocks, int **buffer, int **src){
	int x = (index % width_in_blocks)* ME_blockSize;
	int y = (index / width_in_blocks)* ME_blockSize;

	for (int i = 0; i < ME_blockSize; i++){
		for (int j = 0; j < ME_blockSize; j++){
			buffer[i][j] = src[i + y][j + x];
		}
	}
}


/* copies a block of reference 2d array into a 2d buffer with the size of a block
 * The coordinates on the top left corner of the block should be passed to the function
 */

void get_reference_block(int x, int y, int **buffer, int **ref){

	for (int i = 0; i < ME_blockSize; i++){
		for (int j = 0; j < ME_blockSize; j++){
			buffer[i][j] = ref[i + y][j + x];
		}
	}
}



/*
 * Calculates the SAD between a source block and a reference block. The SAD is reported for all the sub-blocks within the blocks
 */

void calculate_sad(int **src, int ** ref, int **sad){
	int block_sizes = log(ME_blockSize /4)/log(2) +1;
	int max_Sub_Blocks = ME_blockSize * ME_blockSize / 16;

	/* calculates the SAD for sub_blocks of size 4*4 */
	for (int i = 0; i < max_Sub_Blocks; i++) {  /* iterates for all 4*4 sub-blocks within a block */
		sad[0][i] = 0;
		int x = (i % (ME_blockSize / 4)) * 4;
		int y = (i / (ME_blockSize / 4)) * 4;
		/* for each sub-block, compares each pixel and adds the different to the sad */
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				sad[0][i] += abs(src[y + j][x + k] - ref[y + j][x + k]);
			}
		}
	}

	/* calculates the sad for the larger sub-blocks by adding sad of smaller sub-blocks within the block */
	for (int i = 1; i < block_sizes; i++) {
		int subBlocks = (ME_blockSize * ME_blockSize) / (1 << (2 * i + 4));
		int widthInBlocks = ME_blockSize / (1 << (i + 2));
		for (int j = 0; j < subBlocks; j++) {
			int row = 2 * (j / widthInBlocks);
			int col = 2 * (j % widthInBlocks);
			int top    = sad[i - 1][row       * 2 * widthInBlocks + col] + sad[i - 1][row       * 2 * widthInBlocks + col + 1];
			int bottom = sad[i - 1][(row + 1) * 2 * widthInBlocks + col] + sad[i - 1][(row + 1) * 2 * widthInBlocks + col + 1];
			sad[i][j] = top + bottom;
		}
	}


}



/*
 * This Function produces the perfect motion vector for a source block
 * input "cycle" determines the iteration index for the source block currently being processed
 * input "coords" is an array that holds all the coordinates for the reference blocks in the search window, for understanding of organization
 *    of data in coords, look at get_window_coords function
 */

int*** get_golden_output(int ** coords, int cycle, int width_in_blocks, int **src, int **ref){
	int blocks_per_window = ME_windowX * ME_windowY;
	int block_sizes = log(ME_blockSize /4)/log(2) +1;
	int max_Sub_Blocks = ME_blockSize * ME_blockSize / 16;

	int *** output = calloc(ME_blockSize, sizeof(*output));
	int ** src_block = calloc(ME_blockSize, sizeof(*src_block));
	int ** ref_block = calloc(ME_blockSize, sizeof(*ref_block));
	int ** sad = calloc(ME_blockSize, sizeof(*sad));

	for(int i = 0; i < ME_blockSize; i++){
		output[i] = calloc(max_Sub_Blocks, sizeof(**output));
		for(int j = 0; j < max_Sub_Blocks; j++){
			output[i][j] = calloc(3, sizeof(***output));
		}

		src_block[i] = calloc(ME_blockSize, sizeof(**src_block));
		ref_block[i] = calloc(ME_blockSize, sizeof(**ref_block));
		sad[i] = calloc(max_Sub_Blocks, sizeof(**sad));

	}

	/* gets the current source block being analysed */
	get_source_block(cycle, width_in_blocks, src_block, src);


	for( int i = 0; i < blocks_per_window; i++){ /* iterates through all the reference blocks within the current search window */

		/* Using coords input, gets a reference block from the search window */
		get_reference_block(coords[0][cycle*blocks_per_window + i], coords[1][cycle * blocks_per_window + i], ref_block, ref);

		calculate_sad(src_block, ref_block, sad);

		/* calculates the motion vector from the current source block to the reference block */
		int x = i % ME_windowX - (ME_windowX/2);
		int y = i / ME_windowX - (ME_windowY/2);


		/* compares the cost of all the motion vectors and stores the minimum in the output[][][], cost is vector_cost + SAD */
		for(int j = 0; j < block_sizes; j++){
			int sub_blocks = (ME_blockSize * ME_blockSize) / (1 << (2 * j + 4));
			for (int k = 0; k < sub_blocks; k++) {
				if (i == 0) {
					output[j][k][0] = x;
					output[j][k][1] = y;
					output[j][k][2] = sad[j][k];
				} else {
					if (sad[j][k] < output[j][k][2]) {
						output[j][k][0] = x;
						output[j][k][1] = y;
						output[j][k][2] = sad[j][k];
					}
				}
			}
		}

	}

	for(int i = 0; i < ME_blockSize; i++){
		free(src_block[i]);
		free(ref_block[i]);
		free(sad[i]);

	}
	free(src_block);
	free(ref_block);
	free(sad);

	return output;

}

int verify_output(int **src, int **ref, int16_t *mv, int width_in_blocks, int height_in_blocks, int halo_blocks_x, int halo_blocks_y, int totalBufferOutputs, int output_size){

	int mv_index = 0;

	/* Determines the number of bits per output, and calculates the padding required to make it a multiple of 128, needed for alignment in i/o streams */
	int total_bits = 0;
	for(int i = 4; i <= ME_blockSize; i*=2){
		int numBlocks = ME_blockSize * ME_blockSize / (i*i);
		total_bits += numBlocks * 32;
	}
	int padding =( 128 - (total_bits % 128)) /16; //2 bytes per data


	/* Create a 2D array holding the coordinates of items in the output stream corresponding to the 2D matrix organization */
	int **coords = calloc(2, sizeof(*coords));
	for(int i = 0; i < 2; i ++){
		coords[i] = calloc(totalBufferOutputs*ME_numPipes, sizeof(int));
	}
	get_window_coords(coords, width_in_blocks, height_in_blocks, halo_blocks_x, halo_blocks_y);


	int return_code = 1;
	for (int i = 0; i < output_size  && return_code; i++){ /* iterates for each source block */

		/* calculates the motion vector with minimum cost for this source block */
		int ***golden_output = get_golden_output(coords, i, width_in_blocks, src, ref);

		/* compares all the outputs from the DFE to the expected outputs */
		int index = 0;
		for (int j = 4; j <= ME_blockSize  && return_code; j *= 2, index++) {
			int numBlocks = ME_blockSize * ME_blockSize / (j*j);
			for (int k = 0; k < numBlocks*2; k++) {
				if (mv[mv_index] != golden_output[index][k / 2][k % 2]) {
					printf("Output does not match on cycle %d, block size %d x %d, element %d\n", i, j, j, k);
					printf("Expected: (%d, %d), Got: (%d, %d)\n", golden_output[index][k / 2][0], golden_output[index][k / 2][1],
							mv[2*(mv_index/2)], mv[2*(mv_index/2)+1]);
					return_code = 0;
					break;
				}
				mv_index ++;
			}
		}
		free(golden_output);

		mv_index += padding;
	}

	for(int i = 0; i < 2; i ++){
		free(coords[i]);
	}
	free(coords);

	return return_code;
}


int main(int argc, char** args){

	if(argc != 3){
		print_cmd_help();
		return 1;
	}


	/* width and height are passed as program arguments and are in units of blocks */
	int width_in_blocks = atoi(args[1]);
	int height_in_blocks = atoi(args[2]);

	/* determine the number of halo blocks for the boundary calculations, blocks to fit half a window are extended in each direction */
	int halo_blocks_x = (int) ceil( (float)((ME_windowX/2) / ME_blockSize) );
	int halo_blocks_y = (int) ceil( (float)((ME_windowY/2) / ME_blockSize) );

	/* halos should be multiplied by 2 as they are added in both of the opposing directions */
	int width_in_blocks_with_halos = width_in_blocks + 2 * halo_blocks_x;
	int height_in_blocks_with_halos = height_in_blocks + 2 * halo_blocks_y;


	 /* The number of blocks that must be filled in the buffer for the buffer to hold all the search windows of the first row of blocks */
	int fillBlocks = width_in_blocks_with_halos * (2 * halo_blocks_y + 1);
	int fillCycles = fillBlocks * ME_blockSize * ME_blockSize / ME_numPipes;  /* the cycles that should be run till that number of blocks are inputted in the buffer */

	/* total number of elements that are read from the buffer */
	int totalBufferOutputs = (width_in_blocks * height_in_blocks * ME_windowX * ME_windowY) / ME_numPipes;

	/* the program total number of cycles consists of the cycles spent on filling the buffer and reading all the necessary information from it */
	int num_cycles = fillCycles + totalBufferOutputs;

	/* the number of source blocks that should be processed until all the reference blocks are read from the CPU Input.
	 * refer to the CoarseMeBufferContrl state machine to understand the formula
	 */
	int toComplete = width_in_blocks * (height_in_blocks - 1);

	size_t total_size_ref = (width_in_blocks_with_halos) * (height_in_blocks_with_halos) * ME_blockSize * ME_blockSize;
	size_t total_size_src = width_in_blocks * height_in_blocks * ME_blockSize * ME_blockSize;

	/* creating the source and reference arrays, the tiled version is used for passing the data to DFE, 2d version is used for verifying the output */
	int16_t *ref_tiled = calloc(total_size_ref, sizeof(*ref_tiled));
	int16_t *src_tiled = calloc(total_size_src, sizeof(*src_tiled));
	int **ref = calloc(height_in_blocks_with_halos * ME_blockSize, sizeof(*ref));
	int **src = calloc(height_in_blocks * ME_blockSize, sizeof(*src));
	for(int i = 0; i < height_in_blocks_with_halos * ME_blockSize; i++){
		ref[i] = calloc(width_in_blocks_with_halos*ME_blockSize, sizeof(**ref));
		if (i < height_in_blocks * ME_blockSize)
			src[i] = calloc(width_in_blocks*ME_blockSize, sizeof(**src));  //!!  + halo_blocks
	}

	int output_size = width_in_blocks * height_in_blocks;
	int output_typeSize = getOutputTypeSizeInBytes(ME_blockSize);

	int16_t *mv_out = calloc(output_size, output_typeSize);

	create_test_data(width_in_blocks, height_in_blocks, halo_blocks_x, halo_blocks_y, ref, src, ref_tiled, src_tiled);



	max_file_t *maxfile = ME_init();

	max_engine_t *engine = max_load(maxfile, "*");

	ME_actions_t actions;

	actions.ticks_CoarseMeKernel = num_cycles;
	actions.inscalar_CoarseMeKernel_blocksToComplete = toComplete;
	actions.inscalar_CoarseMeKernel_fillBlocks = fillBlocks;
	actions.inscalar_CoarseMeKernel_heightInBlocks = height_in_blocks_with_halos;
	actions.inscalar_CoarseMeKernel_heightPlusHaloInBlocks = height_in_blocks_with_halos;
	actions.inscalar_CoarseMeKernel_numBlocks = width_in_blocks * height_in_blocks;
	actions.inscalar_CoarseMeKernel_widthInBlocks = width_in_blocks;
	actions.inscalar_CoarseMeKernel_widthPlusHaloInBlocks = width_in_blocks_with_halos;
	actions.instream_reference = ref_tiled;
	actions.instream_size_reference = sizeof(*ref_tiled) * total_size_ref;
	actions.instream_source = src_tiled;
	actions.instream_size_source = sizeof(*src_tiled) * total_size_src;
	actions.outstream_mv = mv_out;
	actions.outstream_size_mv = (output_typeSize * output_size);

	printf("running DFE...\n");

	ME_run(engine, &actions);
	printf("done.\n");

	printf("Verifying output...\n");

	int return_code = 0;
	if(verify_output(src, ref, mv_out, width_in_blocks, height_in_blocks, halo_blocks_x, halo_blocks_y, totalBufferOutputs, output_size)){
		printf("*****PASS*****\n");
		return_code = 0;
	}
	else{
		printf("*****FAIL*****\n");
		return_code = 1;
	}


	free(ref);
	free(src);
	free(ref_tiled);
	free(src_tiled);
	ME_free();


	return return_code;
}
