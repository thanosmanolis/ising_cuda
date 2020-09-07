/*
**********************************************
*         - Ising model using CUDA -         *
*    GPU with multiple moments per thread    *
*          Using shared GPU memory			 *
**********************************************
*/

#include "../inc/cuda.h"

//! Define threads per block as a power of 2
//! TILE_DIM must be an integral multiple of BLOCK_ROWS
#define TILE_DIM 32
#define BLOCK_ROWS 4
#define DEPTH 2

void ising(int *G, double *w, int k, int n)
{
    //! Store weights array to GPU
	double *gpu_w;
	cudaMalloc(&gpu_w, 25*sizeof(double));
	cudaMemcpy(gpu_w, w, 25*sizeof(double), cudaMemcpyHostToDevice);

	//! Store G array to GPU
	int *gpu_G;
	cudaMalloc(&gpu_G, n*n*sizeof(int));
	cudaMemcpy(gpu_G, G, n*n*sizeof(int), cudaMemcpyHostToDevice);

	//! GPU array to store the updated values
	int *gpu_G_new;
	cudaMalloc(&gpu_G_new, n*n*sizeof(int));

	//! Temp pointer to swap gpu_G and gpu_G_new
	int *temp;

	//! Blocks & Threads
	dim3 threads( TILE_DIM, BLOCK_ROWS );
	dim3 blocks( n/threads.x, n/threads.x );

	//! Flag to see if changes were made (also store it to gpu to pass it to the kernel)
	int changes_made;
	int *gpu_changes_made;
	cudaMalloc(&gpu_changes_made, (size_t)sizeof(int));

	//! Implement the process for k iterations
	for(int i = 0; i < k; i++)
	{
		//! Initialize changes_made as zero
		changes_made = 0;
		cudaMemcpy(gpu_changes_made, &changes_made, (size_t)sizeof(int), cudaMemcpyHostToDevice);

        kernel<<< blocks , threads >>>(n, gpu_w, gpu_G, gpu_G_new, gpu_changes_made);

        //! Synchronize threads before swapping pointers
		cudaDeviceSynchronize();

		//! Swap pointers for next iteration
		temp = gpu_G;
		gpu_G = gpu_G_new;
		gpu_G_new = temp;

		//! Terminate if no changes were made
		cudaMemcpy(&changes_made, gpu_changes_made,  (size_t)sizeof(int), cudaMemcpyDeviceToHost);
		if(changes_made == 0)
			break;
	}

    //! Copy GPU final data to CPU memory
	cudaMemcpy(G, gpu_G, n*n*sizeof(int), cudaMemcpyDeviceToHost);

	//! Free allocated GPU memory
    cudaFree(gpu_w);
	cudaFree(gpu_G);
	cudaFree(gpu_G_new);
}

__global__ void kernel(int n,  double* gpu_w, int* gpu_G, int* gpu_G_new, int* flag_changes_made)
{
	//! Array in shared memory to store the examined moments
	//! (with their neighbors), every each iteration (part of gpu_G)
	__shared__ int sh_G[(TILE_DIM+2*DEPTH) * (BLOCK_ROWS+2*DEPTH)];

	//! Number of rows and cols of sh_G
	int sh_cols = (blockDim.x + 2*DEPTH);

	//! Array in shared memory to store the weights matrix (gpu_w)
	__shared__ double sh_w[25];

	//! Fill sh_w
	for(int i=0; i<25; i++)
		sh_w[i] = gpu_w[i];

	//! Moment's coordinates
	int mom_row = blockDim.y*blockIdx.y + threadIdx.y;
	int mom_col = blockDim.x*blockIdx.x + threadIdx.x;

	//! Moment's coordinates in shared array
	int sh_row = threadIdx.y + DEPTH;
	int sh_col = threadIdx.x + DEPTH;

	//The step of each thread
	int stepRow = blockDim.y *gridDim.y;
	int stepCol = blockDim.x *gridDim.x;

	//! The indices of the examined neighbors
	int idx_row, idx_col;

	//! Coordinates of the neighbors in shared array
	int neigh_row, neigh_col;

	//! Accessing the spins in the global lattice and "transfer" them in the shared matrix.
	for(int i=mom_row; i<n+DEPTH ;i+=stepRow)
	{
		for(int j=mom_col; j<n+DEPTH;j+=stepCol)
		{
			//! Every thread read its own element in shared memory
			sh_G[sh_row*sh_cols+sh_col] = gpu_G[((i + n)%n)*n + ( (j + n)%n )];

			//! Add left and right neighbors
			if(threadIdx.x < DEPTH)
			{
				neigh_row = sh_row;
				idx_row = (i + n)%n;

				for(int p=0; p<2; p++)
				{
					int adder = (p-1)*DEPTH + p*blockDim.x;
					neigh_col = sh_col + adder;
					idx_col = (j + adder + n)%n;
					sh_G[neigh_row*sh_cols + neigh_col] = gpu_G[idx_row*n + idx_col];
				}
			}

			//! Add top and bottom neighbors
			if(threadIdx.y < DEPTH)
			{
				neigh_col = sh_col;
				idx_col = (j + n)%n;

				for(int p=0; p<2; p++)
				{
					int adder = (p-1)*DEPTH + p*blockDim.y;
					neigh_row = sh_row + adder;
					idx_row = (i + adder + n)%n;
					sh_G[neigh_row*sh_cols + neigh_col] = gpu_G[idx_row*n + idx_col];
				}
			}

			//! Add corner neighbors
			if( (threadIdx.x < DEPTH) && (threadIdx.y<DEPTH) )
			{
				for(int p=0; p<4; p++)
				{
					int adder_row = (p%2 - 1)*DEPTH + (p%2)*blockDim.y;
					neigh_row = sh_row + adder_row;
					idx_row = (i + adder_row + n)%n;

					int adder_col = ((p+3)%(p+1)/2 - 1)*DEPTH + ((p+3)%(p+1)/2)*blockDim.x;
					neigh_col = sh_col + adder_col;
					idx_col = (j + adder_col + n)%n;

					sh_G[neigh_row*sh_cols + neigh_col] = gpu_G[idx_row*n + idx_col];
				}
			}

			//! Synchronize to make sure all threads have added what they were supposed to
			__syncthreads();

			if((i<n)&&(j<n))
			{
				//! Variable to store the value of each moment
				double sum_value = 0;

				//! Iterate through the moment's neighbors (m->X, l->Y axis)
				for(int m=-2; m<3; m++)
			      	for(int l=-2; l<3; l++)
			        {
			            //! Only edit the neighbors of the examined element
			            if((m == 0) && (l == 0))
			                continue;

			            //! Calculate the new value
			            sum_value += sh_w[(2+m)*5 + (2+l)] * sh_G[(m+sh_row)*sh_cols + (l+sh_col)];
			        }

				//! If positive -> 1
				//! If negative -> -1
				if(sum_value > 1e-3)
				{
					gpu_G_new[i*n + j] = 1;
					if(sh_G[sh_row*sh_cols + sh_col] == -1)
						*flag_changes_made = 1;
				}
				else if(sum_value < -1e-3)
				{
					gpu_G_new[i*n + j] = -1;
					if(sh_G[sh_row*sh_cols + sh_col] == -1)
						*flag_changes_made = 1;
				}
				else
					gpu_G_new[i*n + j] = sh_G[sh_row*sh_cols + sh_col];
			}

			//! Synchronize to make sure no thread adds next iteration's values
			__syncthreads();
		}
	}
}
