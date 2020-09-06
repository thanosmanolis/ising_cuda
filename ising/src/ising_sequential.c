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
	int idx_row, idx_col;

	//! Flag to see if changes were made
	int changes_made;

	//! Implement the process for k iterations
	for(int iter=0; iter<k; iter++)
	{
		changes_made = 0;

		//! Access every single moment of the matrix
		for(int i=0; i<n; i++)
			for(int j=0; j<n; j++)
			{
				sum_value = 0;

				//! Iterate through the moment's neighbors (m->row, l->column)
				for(int m=-2; m<3; m++)
					for(int l=-2; l<3; l++)
					{
						//! Only edit the neighbors of the examined element
						if((m == 0) && (l == 0))
							continue;

						//! Find the index of the examined neigbor
			            //! If the element is at a special position (i.e. a corner)
			            //! continue to the other side of the matrix
						idx_row = (i + m + n) % n;
			            idx_col = (j + l + n) % n;

						//! Calculate the new value
						sum_value += w[(2+m)*5 + (2+l)] * G[idx_row*n + idx_col];
					}

				//! If positive -> 1
				//! If negative -> -1
				if(sum_value > 1e-3)
				{
					G_new[i*n + j] = 1;
					changes_made = 1;
				}
				else if(sum_value < -1e-3)
				{
					G_new[i*n + j] = -1;
					changes_made = 1;
				}
				else
					G_new[i*n + j] = G[i*n + j];
			}

		//! Swap pointers for next iteration
		temp = G;
		G = G_new;
		G_new = temp;

		//! Terminate if no changes were made
		if(!changes_made)
		{
			k = iter+1;
			break;
		}
	}

	//! At the last iteration, if the k is odd,
	//! G points to G_new and G_new points to G
	if(k%2 != 0)
		memcpy(G_new, G, n*n*sizeof(int));
}
