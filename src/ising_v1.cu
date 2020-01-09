/*
****************************************
*      - Ising model using CUDA -      *
*    GPU with one thread per moment    *
****************************************
*/

#include "../inc/ising.h"

struct timeval startwtime, endwtime;
double p_time;

#define BLOCKSIZE 16

__global__ void kernel(int n,  double* gpu_w, int* gpu_G, int* gpu_G_new)
{
	//! Moment's coordinates
	int mom_X = threadIdx.x + blockIdx.x * blockDim.x;
	int mom_Y = threadIdx.y + blockIdx.y * blockDim.y;

    //! The indices of the examined neighbors
	int idx_X, idx_Y;

    //! Variable to store the value of each moment
    double sum_value = 0;

	if( (mom_X < n) && (mom_Y < n) )
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

	//! Find a number of threads (per block) that is a divisor of n
	//! Else, threads (per block) take the default value
	int threads = 0;
	for(int i=10; i<=32; i++)
		if(n%i == 0)
			threads = i;
	if(threads == 0)
		threads = BLOCKSIZE;

	//! Grid size comes from a combination of n and threads
	int blocks;
	if(n%threads == 0)
		blocks = n/threads;
	else
		blocks = n/threads + 1;

	printf("threads:%d\n", threads);
	printf("blocks:%d\n", blocks);
	printf("%d\n", blocks*blocks*threads*threads - n*n);

	//! Define block and grid
	dim3 dimBlock( threads, threads );
	dim3 dimGrid ( blocks, blocks );

	//! Implement the process for k iterations
	for(int i = 0; i < k; i++)
	{
        kernel<<< dimGrid , dimBlock >>>(n, gpu_w, gpu_G, gpu_G_new);

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

int main(int argc, char *argv[])
{
	int n, k;

    if(argc > 1)
    {
        n = atoi(argv[1]);  // # number of elements (n*n)
        k = atoi(argv[2]);  // # iterations
    }else
    {
		n = 517;	// default value for n
        k = 1;		// default value for k
    }

	//! Array that will keep the init binary file info
	int *G = (int*)malloc(n*n * sizeof(int));

	//! Weights
    double weights[] = {0.004, 0.016, 0.026, 0.016, 0.004,
                		0.016, 0.071, 0.117, 0.071, 0.016,
            			0.026, 0.117, 0    , 0.117, 0.026,
            			0.016, 0.071, 0.117, 0.071, 0.016,
            			0.004, 0.016, 0.026, 0.016, 0.004};

	//! Open binary file and write contents to G array
    FILE *fptr = fopen("conf-files/conf-init.bin","rb");
    if (fptr == NULL)
	{
        printf("Error opening file");
        exit(1);
    }
    fread(G, sizeof(int), n*n, fptr);
	fclose(fptr);

	//! ========= START POINT =========
    gettimeofday (&startwtime, NULL);

    //! Implement ising procedure
    ising(G, weights, k, n);

	//! ========= END POINT =========
    gettimeofday (&endwtime, NULL);
    p_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
  		      + endwtime.tv_sec - startwtime.tv_sec);

	//! Name of conf file depending on k value
	char filename[25];
	snprintf(filename, sizeof(filename), "conf-files/conf-%d.bin", k);

	//! Compare updated data with the correct data (for k = 1, 4, 11)
	int *data = (int*)malloc(n*n * sizeof(int));
	int isWrong = 0;

	fptr = fopen(filename,"rb");
	fread(data, sizeof(int), n*n, fptr);
	fclose(fptr);
	for(int v = 0; v < n*n; v++)
		if(data[v] != G[v])
			isWrong = 1;

	//! Check if any comparison failed
	if (!isWrong)
		printf("[k=%d]" GREEN " CORRECT\n" RESET_COLOR, k);
	else
		printf("[k=%d]" RED " WRONG\n" RESET_COLOR, k);

	printf(RED "Real Time: %f\n", p_time);

	//! Free allocated GPU memory
    free(G);
	free(data);

    return 0;
}
