/*
****************************************
*      - Ising model using CUDA -      *
*    GPU with one thread per moment    *
****************************************
*/

#include "../inc/cuda.h"

#define BLOCKSIZE 256

__global__ void kernel(int n,  double* gpu_w, int* gpu_G, int* gpu_G_new)
{
	//! Thread ID (from 0 to n*n)
	int thread_id = blockIdx.x*blockDim.x + threadIdx.x;

	//! Variable to store the value of each moment
	double sum_value;

	//! Moment's coordinates
	int mom_X = thread_id%n;
	int mom_Y = thread_id/n;

	//! The indices of the examined neighbors
	int idx_X, idx_Y;

	if( thread_id < n*n )
	{
		//! Iterate through the moment's neighbors (k->X, l->Y axis)
	    for(int k=0; k<5; k++)
	        for(int l=0; l<5; l++)
	        {
	            //! Only edit the neighbors of the examined element
	            if((k == 2) && (l == 2))
	                continue;

	            //! Find the index of the examined neigbor
	            //! If the element is at a special position (i.e. a corner)
	            //! continue to the other side of the matrix
	            idx_X = (mom_X + (k-2) + n) % n;
	            idx_Y = (mom_Y + (l-2) + n) % n;

	            //! Calculate the new value
	            sum_value += gpu_w[l*5 + k] * gpu_G[idx_Y*n + idx_X];
	        }

	    //! If positive -> 1
	    //! If negative -> -1
	    if(sum_value > 1e-3)
	        gpu_G_new[mom_Y * n + mom_X] = 1;
	    else if(sum_value < -1e-3)
	        gpu_G_new[mom_Y * n + mom_X] = -1;
	    else
	        gpu_G_new[mom_Y * n + mom_X] = gpu_G[mom_Y * n + mom_X];
	}
}

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
	
	//! Grid size comes from a combination of n and BLOCKSIZE
	int threads = BLOCKSIZE;
	int blocks;
	if(n*n%threads == 0)
		blocks = n*n/threads;
	else
		blocks = n*n/threads + 1;

	printf(" threads per block: %d\n", threads);
	printf("            blocks: %d\n", blocks);
	printf("  threads in total: %d\n", blocks*threads);

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
