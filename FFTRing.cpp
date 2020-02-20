/*
 *	///////////// ----- FFT - Ring Transformation ----- \\\\\\\\\\\\\\\\\
 *	
 *	The GSW cryptosystem is defined over the ring 
 *							R = Z[x]/\Phi_N(x)
 *	where \Phi_N(x) is the cyclotomic polynomial associated to N. In order
 *	to do fast multiplication we use the Fast Fourier Transform (FFT). 
 *	The FFT is a fast algorithm to compute the Fourier Transform, which is 
 *	a function 
 *					FT:C(x)/x^N-1 ----> C(x)/x^N-1
 *	defined as
 *			FT(a(x)) = (a(\zeta_0), \ldots, a(\zeta_{N-1}))
 *	where \zeta_0, \ldots, \zeta_{N-1} are the N-roots of unity in the 
 *	complex numbers. One last step is to reduce the resulting polynomial
 *	modulo \Phi_N(x). We know that when N is a prime power (N = p^e), the 
 *	corresponding cyclotomic polynomial is
 *		\Phi_N(x) = 1 + x^{N/p} + x^{2N/p} + \ldots + x^{\varphi(N)}.
 *	Therefore, in R we have that
 *		x^{\varphi(N)} = - 1 - x^{N/p} - \ldots - x^{(N-2)N/p}.
 */

#include <iostream>
#include <complex.h>
#include <fftw3.h>
#include "FFTRing.h"
#include "params.h"

namespace FFTRing
{
	/*
	 *	In order to run FFT in parallel it is necessary to declare a 
	 *	different plan and arrays for each thread (alternatively it is 
	 *	possible to declare only a collection of FFT arrays and use a 
	 *	special plan executer).
	 */
	double *re[d_g2];
	fftw_complex *co[d_g2];
	fftw_plan r2c[d_g2], c2r[d_g2];

	void Setup ()
	{
		for (int index=0; index<d_g2; index++)					// Initializing the arrays and plans previously declared
		{
			re[index] = (double*) fftw_malloc(sizeof(double) * N);
			co[index] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * FFTdim);
			r2c[index] = fftw_plan_dft_r2c_1d(N, re[index], co[index], FFTW_PATIENT);
			c2r[index] = fftw_plan_dft_c2r_1d(N, co[index], re[index], FFTW_PATIENT);
		}

	}

	void Ring2FFT (Ring_FFT *coVector, const Ring_ModQ *reVector, const int index)
	{
		/*
		 *	Remember that the size of the polynomial is \phi(N), since it 
		 *	is  in the cyclotomic field. The rest of the entries are set to 
		 *	zero.
		 */

		for (int i=0; i<phi_N; i++)							// Filling the real array with the coefficients of the real polynomial.
		{
			re[index][i] = (double) ((*reVector)[i]).value;
		}

		for (int i=phi_N; i<N; i++)
		{
			re[index][i] = 0.0;
		}

		fftw_execute(r2c[index]);

		for (int i=0; i< FFTdim; i++)
		{
			(*coVector)[i] = co[index][i];
		}
	}

	void FFT2Ring (Ring_ModQ *reVector, const Ring_FFT *coVector, const int index)
	{
		/*
		 *	coVector is the vector the containing evaluation of (half of) 
		 *	the N-roots of unity in C (the set of complex numbers). After 
		 *	filling reVector with the first phi_N coefficients we reduce the 
		 *	rest modulo \Phi_N(x).
		 */

		//std::cerr << "\nraw input :";
		for (int i=0; i<FFTdim; i++)						// Filling the complex array with the entries of the complex vector.
		{
			co[index][i] = (double complex) (*coVector)[i]/N;
		//	std::cerr << creal(co[index][i])*N << " + i" << cimag(co[index][i])*N << "::";
		}
		//std::cerr << "\n\n";



		fftw_execute(c2r[index]);

		//std::cerr << "\nraw output :";
		for (int i=0; i<phi_N; i++)							// Initializing the polynomial with the first phi_N terms
		{
			((*reVector)[i]).value = (long int) round(re[index][i]);
			//std::cerr << std::fixed << round(re[index][i]) << " ";
		}
		//std::cerr << "\n\n";

		for (int i=0; i<(Nprime); i++)						// Reduction modulo \Phi_N(x) of the last terms
		{
			for (int k=0; k<p-1; k++)
			{
				((*reVector)[Nprime*k+i]).value -= (long int) round(re[index][phi_N+i]); 
			}
		}
	}


	void Destructor ()
	{
		for (int index=0; index<d_g2; index++)					// Initializing the arrays and plans previously declared
		{
			fftw_free(re[index]);
			fftw_free(co[index]);
			fftw_destroy_plan(r2c[index]);
			fftw_destroy_plan(c2r[index]);
		}
	}		
}
