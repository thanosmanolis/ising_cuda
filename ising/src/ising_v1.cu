/*
****************************************
*      - Ising model using CUDA -      *
*    GPU with one thread per moment    *
****************************************
*/

#include "../inc/cuda.h"

//! Define threads per block as a power of 2
#define BLOCKSIZE 256

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
	int blocks = (n*n + threads - 1)/threads;

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

__global__ void kernel(int n, double* gpu_w, int* gpu_G, int* gpu_G_new, int* flag_changes_made)
{
	//! Thread ID (from 0 to n*n)
	int thread_id = blockIdx.x*blockDim.x + threadIdx.x;

	//! Variable to store the value of each moment
	double sum_value;

	//! Moment's coordinates
	int mom_row = thread_id/n;
	int mom_col = thread_id%n;

	//! The indices of the examined neighbors
	int idx_row, idx_col;

	if( thread_id < n*n )
	{
		//! Iterate through the moment's neighbors (m->X, l->Y axis)
	    for(int m=-2; m<3; m++)
			for(int l=-2; l<3; l++)
	        {
	            //! Only edit the neighbors of the examined element
	            if((m == 0) && (l == 0))
					continue;

	            //! Find the index of the examined neigbor
	            //! If the element is at a special position (i.e. a corner)
	            //! continue to the other side of the matrix
	            idx_row = (mom_row + m + n) % n;
	            idx_col = (mom_col + l + n) % n;

	            //! Calculate the new value
				sum_value += gpu_w[(2+m)*5 + (2+l)] * gpu_G[idx_row*n + idx_col];
	        }

	    //! If positive -> 1
	    //! If negative -> -1
	    if(sum_value > 1e-3)
		{
	        gpu_G_new[mom_row*n + mom_col] = 1;
			if(gpu_G[mom_row*n + mom_col] == -1)
				*flag_changes_made = 1;
		}
		else if(sum_value < -1e-3)
		{
		    gpu_G_new[mom_row*n + mom_col] = -1;
			if(gpu_G[mom_row*n + mom_col] == -1)
				*flag_changes_made = 1;
		}
	    else
	        gpu_G_new[mom_row*n + mom_col] = gpu_G[mom_row*n + mom_col];
	}
}
