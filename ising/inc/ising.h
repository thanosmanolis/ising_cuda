#ifndef ISING_H
#define ISING_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <sys/time.h>
#include <sys/times.h>

//! Clearing the shell using escape sequences
#define clear() printf("\033[H\033[J")
#define RED "\033[0;31m"
#define GREEN "\033[0;32m"
#define GREEN_BOLD "\033[1;32m"
#define YELLOW "\033[0;33m"
#define RESET_COLOR "\033[0m"

/*
*************************************************
*    @file   ising.h                            *
*    @author amanolis <amanolis@ece.auth.gr>    *
*    @date   Sat Jan 4 16:03:30 2020            *
*    @brief  Main utility functions             *
*************************************************
*/

/*
********************************************************************
*    Ising model evolution                                         *
*                                                                  *
*    - param G    Spins on the square lattice         [n-by-n]     *
*    - param w    Weight matrix                       [5-by-5]     *
*    - param k    Number of iterations                [scalar]     *
*    - param n    Number of lattice points per dim    [scalar]     *
*                                                                  *
*    NOTE: Both matrices G and w are stored in row-major format    *
********************************************************************
*/

void ising(int *G, double *w, int k, int n);

#endif
