#include <iostream>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <cstdlib>
#include "../time_profiler.cpp"
#include <string.h>



int main ()
{
	//int p = 5;
	int N = 512;										// This will be the length of the fourier transform
	//int phi_N = 20;									// The degree of the polynomials we are working with
	int dim_fft_esp = N/2+1;							// The length of the vectors in the FFT space is divided by 2
	int poly[N];										// Polynomial
	long long t1, t2;									// Time beacons


	double* re;											// real array for FFT 
	fftw_complex* copx;									// complex array for FFT
	fftw_plan forw, back;								// FFT plans

	//int exponent = 22;								// This is just self-explanatory...


	/* 
		///// --------	Setup phase	--------	\\\\\
	*/
	t1 = profiler::rdtsc();
	/*
		Filling the polynomial
		Now set as poly(x) = x
	*/
	for (int i = 0; i < N; i++)
	{
		poly[i] = rand();
	}

	/*
		Initializing arrays and plans
	*/
	re = (double*) fftw_malloc(sizeof(double) * N);
 	copx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (dim_fft_esp));
	
	forw = fftw_plan_dft_r2c_1d(N, re, copx, FFTW_PATIENT);
	back = fftw_plan_dft_c2r_1d(N, copx, re, FFTW_PATIENT);
	t2 = profiler::rdtsc();
	std::cout << "Setup time: " << t2 - t1 << "\n";


	/*
		///// -------- 	Executing forward	---------	\\\\\
	*/

	

	t1 = profiler::rdtsc();

	for (int c=0; c<100; c++)
	{
		for (int i = 0; i < N; i++)					// Filling the real array with the coefficients of the polynomial
		{
			re[i] = poly[i];
		}

		/*std::cout << "The input real array is: ";		// Printing the current status...
		for (int i = 0; i < N; i++)						
		{
			std::cout << re[i];
		}
		std::cout << "\n";
		*/


//		std::cout << "computing...\n";
		fftw_execute(forw);		
	}

	for (int i = 0; i < N; i++)					// Filling the real array with the coefficients of the polynomial
	{
		re[i] = poly[i];
	}

	/*std::cout << "The input real array is: ";		// Printing the current status...
	for (int i = 0; i < N; i++)						
	{
		std::cout << re[i];
	}
	std::cout << "\n";
	*/


//	std::cout << "computing...\n";
	fftw_execute(forw);

	t2 = profiler::rdtsc();
	std::cout << "Execution forward time: " << t2 - t1 << "\n";



	/*
	/*
		///// ------- 	Exponentiation	--------	\\\\\\\
	*+/
	t1 = profiler::rdtsc();
	for (int i = 0; i < dim_fft_esp; i++)
	{
		double complex current_entry = (double complex) copx[i];
		for (int j = 0; j < exponent; j++)
		{
			copx[i] = ((double complex) copx[i]) * (current_entry);
		}
	}
	t2 = profiler::rdtsc();
	std::cout << "Exponentiation time: " << t2 - t1 << "\n";

	/*
		///// -------- Execution backwards 	--------	\\\\\\
	*+/


	for (int i = 0; i < N/p; i++)
	{
		copx[i*p] = (double complex) 1.0;
	}

	fftw_execute(back);

	for (int i = 0; i < phi_N; i++)					// Copying back the firts entries of the result... fingers crossed...
	{
		poly[i] = (long int) round(re[i]/N);
		std::cout << poly[i];	
	}

	std::cout << "\n";

	for (int i = 0; i < phi_N; i++)
	{
		if (poly[i]!=0)
		{
			std::cout << i << "\n";
		}
	}

	*/

}