/*!
  \file   tester.c
  \brief  Validate ising model implementation.

  \author Thanos Manolis
  \date   2020-01-04
*/

#include "ising.h"

int main()
{
	//! Matrix dimensions (n) and number of iterations (k)
	int n = 517;
	int k = 1;

	//! Array that will keep the init binary file info
	int *G = malloc(n*n * sizeof(int));

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

    //! Implement ising procedure
    ising(G, weights, k, n);

	//! Name of conf file depending on k value
	char filename[25];
	snprintf(filename, sizeof(filename), "conf-files/conf-%d.bin", k);

	//! Compare updated data with the correct data (for k = 1, 4, 11)
	int *data = malloc(n*n * sizeof(int));
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

    return 0;
}
