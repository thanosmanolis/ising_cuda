/*
****************************************
*      - Ising model using CUDA -      *
*    GPU with one thread per moment    *
****************************************
*/

#include "../inc/cuda.h"

//! Define threads per block as a power of 2
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
	dim3 blocks( (n+threads.x-1)/threads.x, (n+threads.x-1)/threads.x );

	// printf("            blocks: %d\n", blocks.x*blocks.y);
	// printf(" threads per block: %d\n", threads.x*threads.y);
	// printf("moments per thread: %d\n", threads.x/threads.y);
	// printf("  threads in total: %d\n", blocks.x*blocks.y*threads.x*threads.y);

	//! Implement the process for k iterations
	for(int i = 0; i < k; i++)
	{
        kernel<<< blocks , threads >>>(n, gpu_w, gpu_G, gpu_G_new);

        //! Synchronize threads before swapping pointers
		cudaDeviceSynchronize();

		//! Swap pointers for next iteration
		temp = gpu_G;
		gpu_G = gpu_G_new;
		gpu_G_new = temp;
	}

    //! Copy GPU final data to CPU memory
	cudaMemcpy(G, gpu_G, n*n*sizeof(int), cudaMemcpyDeviceToHost);

	//! Free allocated GPU memory
    cudaFree(gpu_w);
	cudaFree(gpu_G);
	cudaFree(gpu_G_new);
}

__global__ void kernel(int n,  double* gpu_w, int* gpu_G, int* gpu_G_new)
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

	//Accessing the spins in the global lattice and "transfer" them in the shared matrix.
	for(int i=mom_row; i<n+DEPTH ;i+=stepRow)
	{
		for(int j=mom_col; j<n+DEPTH;j+=stepCol)
		{
			//Every thread read its own element in shared memory
			sh_G[sh_row*sh_cols+sh_col] = gpu_G[((i + n)%n)*n + ( (j + n)%n )];

			//! Add neighbors if the examined moment is not at a corner
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

			//! Add neighbors if the examined moment is not at a corner
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

			if( (threadIdx.x < DEPTH) && (threadIdx.y<DEPTH) )
			{
				//1st corner (4 spots up and left)
				int sharedDiagAccessorX= (sh_row -DEPTH)*sh_cols +(sh_col-DEPTH);
				int GDiagAccessorX=( ( (i-DEPTH)  + n) % n)*n+( ( (j-DEPTH)  + n) % n);
				sh_G[sharedDiagAccessorX]=gpu_G[GDiagAccessorX];

				//2nd diagonial (4 spots down and left)
				sharedDiagAccessorX= (sh_row+blockDim.y)*sh_cols +(sh_col-DEPTH);
				GDiagAccessorX=( ( (i+blockDim.y)  + n) % n)*n+( ( (j-DEPTH)  + n) % n);
				sh_G[sharedDiagAccessorX]=gpu_G[GDiagAccessorX];

				//3rd corner (4 spots down and right)
				sharedDiagAccessorX= (sh_row+blockDim.y)*sh_cols +(sh_col+blockDim.x);
				GDiagAccessorX=( ( (i+blockDim.y)  + n) % n)*n+( ( (j+blockDim.x)  + n) % n);
				sh_G[sharedDiagAccessorX]=gpu_G[GDiagAccessorX];

				//4rd diagonial (4 spots up and right)
				sharedDiagAccessorX= (sh_row -DEPTH)*sh_cols+(sh_col+blockDim.x);
				GDiagAccessorX=( ( (i-DEPTH)  + n) % n)*n+( ( (j+blockDim.x)  + n) % n);
				sh_G[sharedDiagAccessorX]=gpu_G[GDiagAccessorX];
			}

			//Here we synchronize the block threads in order Shared G values are
			//updated for each thread
			__syncthreads();

			if((i<n)&&(j<n))
			{
				//! Variable to store the value of each moment
				double sum_value = 0;

				//! Iterate through the moment's neighbors (k->X, l->Y axis)
				for(int k=sh_row-2; k<sh_row+3; k++)
			      	for(int l=sh_col-2; l<sh_col+3; l++)
			        {
			            //! Only edit the neighbors of the examined element
			            if((k == sh_row) && (l == sh_col))
			                continue;

			            //! Calculate the new value
			            sum_value += sh_w[(2+k-sh_row)*5 + (2+l-sh_col)] * sh_G[k*sh_cols + l];
			        }

				//! If positive -> 1
				//! If negative -> -1
				if(sum_value > 1e-3)
					gpu_G_new[i*n + j] = 1;
				else if(sum_value < -1e-3)
					gpu_G_new[i*n + j] = -1;
				else
					gpu_G_new[i*n + j] = sh_G[sh_row*sh_cols + sh_col];
			}

			__syncthreads();
		}
	}
}
