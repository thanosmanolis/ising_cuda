/*
********************************
*    Ising model Sequential    *
********************************
*/

#include "../inc/ising.h"

void ising(int *G, double *w, int k, int n)
{
	//! Array to store the updated values
	int *G_new = malloc(n*n * sizeof(int));

	//! Temp pointer to swap G and G_new
	int *temp;

	//! Variable to store the value of each moment
	double sum_value = 0;

	//! The indices of the examined neighbors
	int idx_X, idx_Y;

	//! Implement the process for k iterations
	for(int i = 0; i < k; i++)
	{
		//! Access every single moment of the matrix
		for(int x=0; x<n; x++)
			for(int y=0; y<n; y++)
			{
				sum_value = 0;

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
						idx_X = (x + (k-2) + n) % n;
						idx_Y = (y + (l-2) + n) % n;

						//! Calculate the new value
						sum_value += w[l*5 + k] * G[idx_Y*n + idx_X];
					}

				//! If positive -> 1
				//! If negative -> -1
				if(sum_value > 1e-3)
					G_new[y * n + x] = 1;
				else if(sum_value < -1e-3)
					G_new[y * n + x] = -1;
				else
					G_new[y * n + x] = G[y * n + x];
			}

		//! Swap pointers for next iteration
		temp = G;
		G = G_new;
		G_new = temp;
	}

	//! At the last iteration, if the k is odd,
	//! G points to G_new and G_new points to G
	if(k%2 != 0)
		memcpy(G_new, G, n*n*sizeof(int));
}
