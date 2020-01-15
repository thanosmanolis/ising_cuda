/*
************************************
*    - Ising model using CUDA -    *
*           Main function          *
************************************
*/

#include "ising.h"

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

	//! If n!=517, fill the matrix with random starting values
	//! Else, fill the matrix using the conf-init bin file
	if(n == 517 && (k==1 || k==4 || k==11))
	{
		//! Open binary file and write contents to G array
	    FILE *fptr = fopen("conf-files/conf-init.bin","rb");
	    if (fptr == NULL)
		{
	        printf("Error opening file");
	        exit(1);
	    }
	    fread(G, sizeof(int), n*n, fptr);
		fclose(fptr);
	}else
	{
		// ! Fill G with random numbers (-1 or 1)
		for(int i=0; i<n*n; i++)
		{
			G[i] = (int)rand();
			if(G[i]%2)
				G[i] = 1;
			else
				G[i] = -1;
		}
	}

	//! ========= START POINT =========
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start);

    //! Implement ising procedure
    ising(G, weights, k, n);

	//! ========= END POINT =========
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	float time = 0;
	cudaEventElapsedTime(&time, start, stop);

	//! Print runtime to file or terminal
	if(argc > 3)
	{
		//! Write output to file
		FILE *f = fopen(argv[3], "a");
	    if (f == NULL)
	    {
	        printf("Error opening file!\n");
	        exit(1);
	    }
		fprintf(f, "[%d,%d]: %f sec\n", n, k, time/1000);
	    fclose(f);
	}else
		printf(YELLOW "[%d,%d]: %f sec\n" RESET_COLOR, n, k, time/1000);

	//! Check if result is correct
	if(n == 517 && (k==1 || k==4 || k==11))
	{
		//! Name of conf file depending on k value
		char filename[25];
		snprintf(filename, sizeof(filename), "conf-files/conf-%d.bin", k);

		//! Compare updated data with the correct data (for k = 1, 4, 11)
		int *data = (int*)malloc(n*n * sizeof(int));
		int isWrong = 0;

		FILE *fptr_2 = fopen(filename,"rb");
		fread(data, sizeof(int), n*n, fptr_2);
		fclose(fptr_2);
		for(int i = 0; i < n*n; i++)
			if(data[i] != G[i])
				isWrong = 1;

		//! Check if any comparison failed
		if (!isWrong)
			printf(GREEN "CORRECT\n" RESET_COLOR);
		else
			printf(RED "WRONG\n" RESET_COLOR);

		//! Free allocated memory
		free(data);
	}

	//! Free allocated memory
    free(G);

    return 0;
}
