#ifndef PARAMSGEN_H
#define PARAMSGEN_H

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdint.h>

/* ---------- Constants ---------- */

/*
 * We will base the cryptosystem on the cyclotomic field
 * 				Z[x]/\Phi_N(x)
 * For our purposes, we will set 
 * 					N = p^e	
 * where p is prime.
 */

const int p = 5;
const int p_prime = 5;
const int e = 4;
const int q = 625;									// q  = 5^4
const int Nprime = 625;								// N' = N/p
const int N = Nprime*p;								// N = p^e

/*
 * The dimension of the GSW instance is given by the degree of 
 * the cyclotomic polynomial, which is the Euler totient function
 * of N that is given by
 *					(p-1)p^(e-1)
 */
const int phi_N = N - Nprime;

/*
 *	The dimension of the output of FFT (r2c) is N/2+1. See
 * 
 *	http://www.fftw.org/doc/What-FFTW-Really-Computes.html
 *	http://www.fftw.org/doc/The-1d-Real_002ddata-DFT.html
 *
 *	for more details.
 */
const int FFTdim = (N/2)+1;

const int n = 500;									// The dimension of the LWE instance

/*
 * The cyclotomic polynomial on N is the only irreducible polynomial 
 * with integer coefficients that divides
 *						x^N - 1
 * but does not divide any polynomial of the form x^M - 1 for M<N. 
 * In this case it can be written as
 * 				\sum_{i=0}^{p-1} x^{i*N/p}
 */

const int d_Q = 32;
const long int Q = (long int) 1 << d_Q;				// Q  = 2^32



/* ------- Type Definitions ------- */

//typedef int32_t ZmodQ;							// Defining the type ZmodQ as an integer of 32 bits. ALERT: This type must change whenever Q changes.
struct uZmodQ 										// u is for unsigned
{
	uint64_t value : d_Q;
};					

struct ZmodQ 
{
	int64_t value : d_Q;
};


typedef ZmodQ Ring_ModQ[phi_N];
typedef fftw_complex Ring_FFT[FFTdim];

/* 
	-------- Gadget Matrix --------
 
	Attention:
	The variable v is called u on the paper 
*/
const int d_g = 4;
const int e_g = 8;

const int d_g2 = 2*d_g;

const ZmodQ v = { int (Q/(p*p_prime)) };       		// Q/p^2 (floor)

//const ZmodQ v = {171798691};       				// Q/p^2 (floor)
//const ZmodQ v_inverse = {204522251};				// v^{-1} mod Q

const ZmodQ vB_g[d_g] = {v, {v.value << e_g}, {v.value << (2*e_g)}, {v.value << (3*e_g)}};
const int g_bits[d_g] = {8, 8, 8, 8};				// Used for the descomposition of the accumulator
const int g_bits_32[d_g] = {d_Q - g_bits[0], d_Q - g_bits[1], d_Q - g_bits[2], d_Q - g_bits[3]};

/* ----- Key-Switch Parameters ----- */

const int B_ks = 25;
const int d_ks = 7;
const ZmodQ KS_table[d_ks] = 
{	
	{1},
	{25},
	{25*25},
	{25*25*25},
	{25*25*25*25},
	{25*25*25*25*25},
	{25*25*25*25*25*25}
};													// The first 7 powers of 25

/* ----- Base-Switch Parameters ----- */

const int B_r = 25;
const int d_r = 2;
const int Br_table[d_r] = {1,25};



/*
 * The increase routine is implemented in parallel using 6 different 
 * sections, this fixes the number of threads to use to 6
 */
const int NumThreads = 8;



#endif
