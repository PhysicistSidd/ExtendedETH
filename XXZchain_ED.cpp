#include <iostream>
#include <bitset>
#include <string>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <chrono>
#include <math.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

const unsigned int N = 5;
const unsigned long dimension = ( ((long) 1) << N);

const double j1 = 1.0;
const double j2 = 0.85;
const double j1z = -0.05;
const double j2z = -0.05;
const double hz1 = 0.5;
const double hz2 = -1.2;

/*
 * Fast, Stable implementation of the choose function
 * (p choose q)
 */
unsigned int Choose(unsigned int p, unsigned int q){
    if (q > p) return 0;
    if (q * 2 > p) q = p-q;
    if (q == 0) return 1;

    int result = n;
    for( int i = 2; i <= q; ++i ) {
        result *= (p-i+1);
        result /= i;
    }
    return result;
}

/**
 * @brief Computes the Hilbert space dimension for each number sector
 * 
 * @param NCn Stores Hilbert space dimensions
 */
void computeDimension(unsigned int * NCn){
    for ( int n = 0; n <= N; n++){ //N+1 total sectors from zero spins up to all spins up
        *(dn + n) = Choose(N,n); //The Hilbert space dimension is just the choose function
    }
}

/**
 * @brief Computes the integers which represent the index in the full Hilbert space for states in the sector specific hilbert space
 * 
 * @param n the number of spins pointing up
 * @param dn dimension in n spin sector
 * @param stateIndices stores the indices for each state
 * 
 * Example: for a 2 site system (in binary) the index for the state with only the right spin pointing up is '01' (or 1). 
 * However, its index in the sector is simple 0. The array that store the indices for the n=1 sector would be
 * stateIndices[1] = { 1, 2}, thus the zeroth element is 1 and the first element is 2 ('10').
 * 
 * Another example would be the 3 spin sector of 5 site chain: 
 * stateIndices[3] = { 7, 11, 13, 14, 19, 21, 22, 25, 26, 28}
 * 
 */
void computeStateIndices(unsigned int n, unsigned int dn, unsigned int * stateIndices){
    unsigned int index = 0; //full space index
    unsigned int nextIndex; //needed for updating
    stateIndices = new unsigned int[dn]; 

    for ( unsigned int m = 0; m < n; m++ ){
        index += (1 << m); //initialize to all up spins on the right end
    }

    for ( unsigned int sectorIndex = 0; sectorIndex < dn; sectorIndex++){
        stateIndices[sectorIndex] = index; //store full space index to sector index

        /*
         * Incomprehensible but wicked fast implementation for bitwise updating of indices. See: https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
         */
        unsigned int t = index | (index - 1); // bitwise 'or' isolates most significant bit in t=1
        nextIndex = (t+1) | (((~t & -~t) - 1) >> (__builtin_ctz(index) + 1)); // right side of the 'or' updates least significant bits

        index = nextIndex; //update
    }
}

/**
 * @brief 
 * 
 * @param n 
 * @param dn 
 * @param Hamiltonian 
 * @param stateIndices 
 */
void computeHamiltonianElements(unsigned int n, unsigned int dn, MatrixXcd * Hamiltonian, unsigned int * stateIndices){
    for (int in = 0; in < dn; in++){
        for (int out = 0; out < dn; out++){

        }
    }
}

int main(){
    unsigned int n; //sector number
    unsigned int * dn = new unsigned int[N+1]; //dimension of said sector
    unsigned int * stateIndices[N+1];

    MatrixXcd * Hamiltonians = new MatrixXcd[N+1];
    MatrixXcd * EnergyStates = new MatrixXcd[N+1];
    VectorXd * EnergyLevels = new VectorXcd[N+1];

    computeDimensions(dn);

    for ( n = 0 ; n < N+1 ; n++ ){
        computeStateIndices(n, dn[n], stateIndices[n]);
        computeHamiltonianElements(n, dn[n], (Hamiltonians + n), stateIndices[n]);
    }
}