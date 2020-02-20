/*
 *	/////////////// ----- Homomorphic Accumulator ----- \\\\\\\\\\\\\\\\\
 *	
 *	The Homomorphic Accumulator is a bootstrapping technique to refresh the
 *	LWE ciphertexts. It is based on the Gentry-Sahai-Waters (GSW) 
 *	cryptosystem.
 *	
 */

#include <iostream>
#include <cstdlib>
#include <omp.h>
#include "FFTRing.h"
#include "HomACC.h"
#include "tools.h"
 
#include <complex.h>
#include <fftw3.h>


namespace HomACC 
{

	/*
	 *	Each state is recovered by a characteristic function: a function 
	 *	that returns 1 if the element belongs to a set and 0 otherwise (in 
	 *	this context is 1 if a variable has a certain value and 0 
	 *	otherwise). Each one of the stater requires a different test 
	 *	function and, therefore, a different test vector.
	 *	
	 *	A Ring_ModQ element is an polynomial on the ring
	 *					R = Z[x]/\Phi_N(x)
	 *	Therefore it can be seen as a vector of \varphi(n) (phi_N) elements.
	 */


	Ring_ModQ characFuncVector_real[p];						// Array of vectors that whose entry value is p-1 on a specific set and -1 everywhere else
	ZmodQ v_inverse;

	/*
	 *	\\\\\\\\\\\\\\\\\\\\ ----- Setup ----- ////////////////////////
	 *
	 *	The Setup step declares the necessary elements for the rest of the 
	 *	functions.
	 *		-Set a time-depending seed for the rand() function
	 *		-Call the Setup function for FFTRing 
	 *		-Compute the test vectors (in FFT form)
	 */

	void Setup ()
	{
		srand(time(NULL));
		v_inverse = tools::inv_modQ(v);			// Write inv_modQ

		
		FFTRing::Setup();									// Call the Setup for the FFT-Ring translator. The arrays and plans are declared at this step
		std::cerr << "FFT Setup ready... \n";

		

		for (int i=0; i<p-1; i++)
		{
			for (int j=0; j<i*(Nprime); j++)				
			{
				characFuncVector_real[i][j].value = -1;
			}
			for (int j=i*(Nprime); j<(i+1)*(Nprime); j++)
			{
				characFuncVector_real[i][j].value = p-1;
			}
			for (int j=(i+1)*(Nprime); j<phi_N; j++)
			{
				characFuncVector_real[i][j].value = -1;
			}
		}
		
		for (int j=0; j<phi_N; j++)
		{
			characFuncVector_real[p-1][j].value = -1;
		}

		std::cerr << "Characteristic function vectors initialized... \n";

	}

	/*
	 *	\\\\\\\\\\\\\\\\\\\\\ ----- Y_m_coeff ----- /////////////////////
	 *
	 *	The polynomial Y=X^{N/q} is a qth root of unity in R=Z[x]/Phi_N(x). 
	 *	Given a number m\in Z_q, this function computes an integer array 
	 *	with the coefficients of the polynomial. It is used in the Encrypt
	 *	and InitializeACC functions.
	 */
	void Y_m_coeff (int m, int* Y_m)
	{
		m = ((m%q)+q)%q;									// Sanitize
		for (int i=0; i<phi_N; i++)							// Fill the vector with zeros
		{
			Y_m[i] = 0;
		}
		if (m*(N/q)<phi_N)									// m<phi_N
		{
			Y_m[m*(N/q)] = 1;
		}
		else 												// phi_N <= m < N
		{
			int xs = m*(N/q) - phi_N;
			for (int i=0; i<p-1; i++)
			{
				Y_m[i*Nprime+xs] = -1;
			}
		}
	}

	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\ ----- Encrypt ----- //////////////////////
	 *
	 */
	void Encrypt (ct_FFT ct, Ring_FFT sk_FFT, int m)
	{
		Ring_FFT a_i;										
		ct_ModQ ct_Q;										// GSW ciphertext
		m = (((m % q) + q) % q);							// The message is first reduced mod q (in the set {0, 1, ..., q-1})

		int* Y_m;
		Y_m = (int*) malloc(sizeof(int)*phi_N);
		Y_m_coeff (m, Y_m);									// Precomputing Y^m 

		for (int i=0; i<d_g2; i++)
		{
			for (int j=0; j<phi_N; j++)
			/*
			 *	Fill the corresponding entry on the first column with 
			 *	a random vector. 
			 *	TODO: For cryptographic purposes use a cryptographic PRG
			 *	instead of rand().
			 */
			{
				ct_Q[i][0][j].value = rand(); 
			}
			FFTRing::Ring2FFT(&a_i, &ct_Q[i][0]);			// Push the vector to the FFT space

			for (int j=0; j<FFTdim; j++)					// Multiply the vector (as polynomials) by the corresponding secret key
			{
				a_i[j] = ((double complex) a_i[j]) * ((double complex) sk_FFT[j]);
			}
			FFTRing::FFT2Ring(&ct_Q[i][1], &a_i);			// Writing the result on the second column

			for (int j=0; j<phi_N; j++)						// Adding a small error on the vector
			{
				ct_Q[i][1][j].value += Sample(Chi1);
			}
		}
		
		for (int i=0; i<d_g; i++)
		/*
		 *	Add a multiple of the gadget matrix. The gadget matrix is 
		 *	defined as 
		 *		G = (I B_gI ... B_g^{d_g-1}I)
		 *	Where I is the 2x2 identity matrix (in the ring Z[x]/Phi_N(x)).
		 *	The corresponding multiple is 
		 *		uY^m
		 */
		{
			for (int j=0; j<phi_N; j++)
			{
				ct_Q[2*i  ][0][j].value += Y_m[j]*vB_g[i].value;
				ct_Q[2*i+1][1][j].value += Y_m[j]*vB_g[i].value;
			}
		}

		for (int i=0; i<d_g2; i++)							//	Push the result to the FFT space writing it on the FFT 	ciphertext ct
		{
			for (int j=0; j<2; j++)
			{
				FFTRing::Ring2FFT(&ct[i][j], &ct_Q[i][j]);
			}
		}
	}


	/*
	 *	\\\\\\\\\\\\\\\\\---- Eval-Key Generation ----///////////////////
	 *
	 *	An FHEW evaluation key consists of two components: a boostrapping 
	 *	key and LWE key-switch key.
	 *	
	 *		* A bootstrapping key is a 3-dimensional array of pointers to 
	 *			FHEW-ciphertexts in FFT format.
	 *		* An LWE key-switch key is a 3-dimensional array of pointers
	 *			to LWE-ciphertexts over the Q-space (Z_Q).
	 */

	void KeyGen (EvalKey *EK, const LWE::SecretKey LWEsk)
	{
		/*
		 *	The fist part creates a key with wich the LWEsk will be 
		 *	encrypted (Such a key is called z in the paper)
		 */
		LWE::SecretKey_phiN HomACCsk;
		LWE::KeyGen_QphiN (HomACCsk);						// Generate a local LWE key of length phi_N 

		/*
		 *	The second part creates and writes the Key-Switch key, writing 
		 *	it in FFT format
		 */

		LWE::SwitchingKeyGen (EK->KSkey, LWEsk, HomACCsk);	// Using the key just generated, generate and write the KeySwitch key

		Ring_FFT HomACCskFFT;
		FFTRing::Ring2FFT(&HomACCskFFT, &HomACCsk);

		std::cout << "Key-Switch Key generated ...\n";
		
		/*
		 *	The third part creates and writes the bootstrapping key
		 */
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<B_r; j++)
			{
				for (int k=0; k<d_r; k++)
				{
					EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT));
					Encrypt((*EK->BSkey[i][j][k]), HomACCskFFT, LWEsk[i] * j * Br_table[k]);
				}
			}
		}
		std::cout << "Bootstrapping Key generated ...\n";
	}


	/*
	 *	\\\\\\\\\\\\\---- Adaptative Eval-Key Generation ----//////////////
	 *
	 *	The previous function generates an evaluation key from a LWE key 
	 *	that is generated inside the function itself. Such key cannot be 
	 *	accessed from any other point of the program. 
	 *	
	 *	This function is a variant where the key is generated outside. 
	 *	It can be used for experimental or debugging purposes.
	 */

	void KeyGen_ad (EvalKey *EK, const LWE::SecretKey LWEsk, const LWE::SecretKey_phiN *HomACCsk)
	{
		/*
		 *	This part creates and writes the Key-Switch key, writing 
		 *	it in FFT format
		 */

		LWE::SwitchingKeyGen (EK->KSkey, LWEsk, *HomACCsk);	// Using the key HomACCsk, generate and write the KeySwitch key

		std::cerr << "switching key generated \n";

		Ring_FFT HomACCskFFT;
		FFTRing::Ring2FFT(&HomACCskFFT, HomACCsk);

		std::cout << "Key-Switch Key generated ...\n";
		
		/*
		 *	The third part creates and writes the bootstrapping key
		 */
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<B_r; j++)
			{
				for (int k=0; k<d_r; k++)
				{
					EK->BSkey[i][j][k] = (ct_FFT*) fftw_malloc(sizeof(ct_FFT));
					Encrypt((*EK->BSkey[i][j][k]), HomACCskFFT, LWEsk[i] * j * Br_table[k]);
				}
			}
		}
		std::cout << "Bootstrapping Key generated ...\n";
	}




	/*
	 *	\\\\\\\\\\\\\\\\\ ---- ACCIncr_by_thread ---- ///////////////////
	 *
	 *	On this method we describe the individual work that each thread does
	 *	to carry out the routine of adding information to the accumulator in 
	 *	parallel. Basically works on each column of the accumulator, by first
	 *	decomposing the entries of the column and then multiplying by the 
	 *	corresponding row of C. The method will be called by 
	 *	ACCIncr_parallel_trigger. 
	 */


	void ACCIncr_by_thread (ct_FFT ACC, ct_FFT C, int i)
	{
		/*
		 *	About the following variables:
		 *		* realACC is used to pull the accumulator to the real space. 
		 *			Remember that the ACC is passed to this function as an 
		 *			element of the fourier space.
		 *		* decRealACC is the decomposed form of the realACC. 
		 *		* decFFTACC is used to carry the FFT forms of the decomposed 
		 *			ACC.
		 */

		Ring_ModQ realACC, decRealACC[d_g];
		Ring_FFT decFFTACC[2][d_g];
		
		for (int j = 0; j < 2; j++) 
		{

			FFTRing::FFT2Ring(&realACC, &ACC[i][j], i);		// Pull the accumulator to the real space
    
			for (int k = 0; k < phi_N; k++)
			{
				/*
				 *	Each entry of the vector is a double (see FFTRing.cpp). 
				 *	After the following loop the integer t will be 
				 *	decomposed into 
				 *		t = r_0(2^0) + r_1(2^11) + r_2(2^22).
				 *
				 *	The integers r_0, r_1 and r_2 will form the entries of 
				 *	the decomposed form of ACC.
				 */
				ZmodQ t = {(realACC[k]).value * (v_inverse).value};       
				for (int l = 0; l < d_g; l++) 
				{
					/*
					 *	After performing the first two bit-shifts the 
					 *	integer r (which is a 32-bit integer) ends up taking 
					 *	the integer value of the corresponding chunk of 
					 *	integers extracted from t.
					 */
					ZmodQ r 					= {t.value 				<< g_bits_32[l]};      
					r.value 					= r.value 				>> g_bits_32[l];
					t.value 					= (t.value - r.value) 	>> g_bits[l];
					(decRealACC[l][k]).value	= r.value;
				}
			}

			for (int l = 0; l < d_g; l++)
			{												// Push the accumulator back to the FFT space.
				FFTRing::Ring2FFT(&decFFTACC[j][l], &decRealACC[l], i);
			}
		}

		/*
		 *	At this point we have the decomposed forms of the accumulator in 
		 *	FFT format. On the following loop we are performing N2 matrix 
		 *	multiplications (or, well, N2 dot products since this is only 
		 *	what a single thread does.)
		 */

		for (int j = 0; j < 2; j++) 
		{
			for (int k = 0; k < FFTdim; k++) 
			{
				ACC[i][j][k] = (double complex) 0.0;
				for (int l = 0; l < d_g; l++) 
				{
					ACC[i][j][k] += ((double complex) decFFTACC[0][l][k]) * ((double complex) C[2*l][j][k]);
					ACC[i][j][k] += ((double complex) decFFTACC[1][l][k]) * ((double complex) C[2*l+1][j][k]);
				}
			}
		}
	}

	/*
     *	\\\\\\\\\\\\\\\\ ---- ACCIncr_parallel_trigger ---- ////////////////
     *
	 *	On this function we trigger the calls to the AddToACC_on_thread 
	 *	function in parallel.
	 *
	 *	Use of OpenMP: After several experiments we concluded that the best 
	 *	mode to divide the work is using omp parallel sections. A single 
	 *	thread will work on each section, sharing the variable ACC, but 
	 *	using an initialized private copy of C (for more information see the 
	 *	OpenMP documentation for firstprivate). At the end of the method, 
	 *	the accumulator will have all the needed information from C.
	 */
	void ACCIncr_parallel_trigger(ct_FFT ACC, ct_FFT C) 
	{  

		omp_set_num_threads(NumThreads);                           // Fix the number of thread to use to 6, since we are dividing into 6 sections. Change to whatever number you want for experimental purposes... or not.

		#pragma omp parallel sections firstprivate(C)     // Apparently this is going to make the thing faster...
		{

			#pragma omp section 							// Section 1
			{
				ACCIncr_by_thread(ACC, C, 0);
			}												// End of Section 1

			#pragma omp section 							// Section 2
			{
				ACCIncr_by_thread(ACC, C, 1);
			}												// End of Section 2

			#pragma omp section 							// Section 3
			{
				ACCIncr_by_thread(ACC, C, 2);
			}												// End of Section 3

			#pragma omp section 							// Section 4
			{
				ACCIncr_by_thread(ACC, C, 3);
			}												// End of Section 4 

			#pragma omp section 							// Section 5
			{
				ACCIncr_by_thread(ACC, C, 4);
			}												// End of Section 5

			#pragma omp section 							// Section 6
			{
				ACCIncr_by_thread(ACC, C, 5);
			}												// End of Section 6

			#pragma omp section 							// Section 5
			{
				ACCIncr_by_thread(ACC, C, 6);
			}

			#pragma omp section 							// Section 5
			{
				ACCIncr_by_thread(ACC, C, 7);
			}

		}													// End of pragma sections
	}

	/*
	 *	\\\\\\\\\\\\\\\\\\\\\ ---- InitializeACC ---- /////////////////////
	 *
	 *	Set a ciphertext to X^m * G (noiseless encryption of m)
     */
	
	void InitializeACC(ct_FFT ACC, int m) 
	{

		//Ring_FFT a_i;
		ct_ModQ ct_Q;
		m = ((( (m) % q) + q) % q);							// The message is first reduced mod q

		int* Y_m;
		Y_m = (int*) malloc(sizeof(int)*phi_N);
		Y_m_coeff (m, Y_m);									// Precomputing Y^m


		for (int i = 0; i < d_g2; ++i)						// Fill up the matrix with zeros
		{
			for (int j = 0; j < 2; ++j) 
			{
				for (int k = 0; k < phi_N; ++k) 
				{
					ct_Q[i][j][k].value = 0;
				}
			}
		}

		for (int i=0; i<d_g; i++)
		{
			for (int j=0; j<phi_N; j++)
			{
				ct_Q[2*i  ][0][j].value += Y_m[j]*(vB_g[i]).value;
				ct_Q[2*i+1][1][j].value += Y_m[j]*(vB_g[i]).value;
			}
		}

		for (int i=0; i<d_g2; ++i) 
		{
			for (int j=0; j<2; ++j) 						//	Push the result to the FFT space writing it on the FFT 	ciphertext ct
			{
				FFTRing::Ring2FFT(&ACC[i][j], &ct_Q[i][j]);
			}
		}
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\ ---- vector2matrix ---- ////////////////////
	 *	
	 *	Given a vector a, the function f_a:R\rightarrow R defined by
	 *						f_a(b) = ab
	 *	is a linear mapping. Therefore a has an associated matrix a_arrow
	 *	that uniquely defines f_a.
	 *
	 *	This function computes the associated matrix corresponding to the 
	 *	the function f_a.
	 */

	void vector2matrix (Ring_ModQ a, ZmodQ ** a_arrow)
	{
		for (int i=0; i<phi_N; i++)							// Filling the first column with the coefficients of the vector a
		{
			a_arrow[i][0] = a[i];
		}	
		for (int j=1; j<phi_N; j++)							// Computing the rest of the columns recursively
		{
			(a_arrow[0][j]).value = 0;							// Initializing the column by shifting the previous one
			for (int i=1; i<phi_N; i++)						
			{
				a_arrow[i][j] = a_arrow[i-1][j-1];
			}

			for (int i=0; i<p-1; i++)						// Substracting the corresponding multiple of x^{\phi_N}
			{
				a_arrow[i*(Nprime)][j].value -= a_arrow[phi_N-1][j-1].value;
			}
		}
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\ ---- vec_x_mat ---- /////////////////////
	 *
	 *	This method allows us to compute the matrix multiplication 
	 *						t * M
	 *	where t is a row vector (of phi_N entries) and M is a square matrix 
	 *	of dimension phi_N by phi_N.
	 */

	void vec_x_mat (Ring_ModQ t, ZmodQ ** M, Ring_ModQ res)
	{
		for (int i=0; i<phi_N; i++)
		{
			(res[i]).value = 0;
			for (int j=0; j<phi_N; j++)
			{
				res[i].value += t[j].value * (M[j][i]).value;
			}
		}
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\ ---- innerProd ---- ////////////////////
	 *
	 *	A short method to compute the inner product of two elements in the 
	 *	ring
	 */

	ZmodQ innerProd (Ring_ModQ a, Ring_ModQ b)
	{
		ZmodQ innp = {0};
		for (int i=0; i<phi_N; i++)
		{
			innp.value += (a[i].value)*(b[i].value);
		}
		return innp;
	}



	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\ ---- innerProd ---- ////////////////////
	 *
	 *	A short method to compute the inner product of two elements in the 
	 *	ring
	 */

	int innerProdInt (int a[n], int b[n])
	{
		int innp = 0;
		for (int i=0; i<n; i++)
		{
			innp += (a[i])*(b[i]);
		}
		return innp;
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\\ ---- IncrLoop ---- //////////////////////
	 *
	 *	This function is a routine to exclusively call all the 
	 *	corresponding increasing steps. 
	 */

 	void IncrLoop (LWE::CipherText* ctRaw, ct_FFT ACC, const EvalKey* EK)
	{
		for (int i = 0; i < n; ++i) 
		{
			int ai = (q - (ctRaw->a)[i] % q) % q;

			for (int k = 0; k<d_r; ++k, ai /= B_r) 
			{
				int a0 = ai % B_r;

				if (a0) 
				{
					ACCIncr_parallel_trigger(ACC, *(EK->BSkey[i][a0][k]));
				}
			}
		}

	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\ ---- MembershipTest ---- ////////////////////
	 *
	 *	This function allows us to evaluate the membership of a non-zero 
	 *	coordinate in an specific set of Z_{phi_N} 
	 */

	void MembershipTest (LWE::CipherTextQ_phiN *ct_Q_phiN, ct_FFT ACC, int l)
	{
		l = ((l % p) + p) % p;								// Sanitizing l

		/*
		 *	About the following variables:
		 *		- a is used to store the ring-value of entries of the
		 *			accumulator
		 *		- t is the corresponding test vector (computed from the 
		 *			vectors previously declared)
		 *		- txa is the result of the product t a^arrow (as matrices)
		 */

		Ring_ModQ a, t, txa;										
		ZmodQ ** a_arrow, innp;


		//compute the corresponding test vector here
		for (int i=0; i<phi_N; i++)
		{
			t[i] = characFuncVector_real[l][i];
		}

		a_arrow = (ZmodQ**)malloc(sizeof(ZmodQ*) * phi_N);
		for (int i=0; i<phi_N; i++)
		{
			a_arrow[i] = (ZmodQ*) malloc(sizeof(ZmodQ) * phi_N);
		}


		FFTRing::FFT2Ring(&a, &ACC[1][0], 0);				// Check the zero here
		vector2matrix(a, a_arrow);
		vec_x_mat(t, a_arrow, txa); 

		for (int i=0; i<phi_N; i++)
		{
			ct_Q_phiN->a[i] = txa[i];
		}

		FFTRing::FFT2Ring(&a, &ACC[1][1]);
		innp = innerProd(t, a);

		(ct_Q_phiN->b).value = innp.value + v.value;

	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\\ ---- Collapse ---- //////////////////////
	 *
	 */	

	void Collapse (ct_FFT ACC, EvalKey* EK, LWE::CipherText * ctFresh, int l)
	{
		if (l<0 or l>=p)
		{
			std::cerr << "Variable out of rank.\n";
			ctFresh = NULL;
			return;
		}

		LWE::CipherTextQ ctQ;
		LWE::CipherTextQ_phiN ctQ_phiN;

		MembershipTest(&ctQ_phiN, ACC, l);
		LWE::KeySwitch(&ctQ, EK->KSkey, ctQ_phiN);
		LWE::ModSwitch(ctFresh, ctQ);
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\\ ---- MSBTest ---- //////////////////////
	 *
	 *	TODO: Test this function.
	 */

	void MSBTest (ct_FFT ACC, EvalKey* EK, LWE::CipherText * ctFresh)
	{
		if (p != 2)
		{
			std::cerr << "Ciphertext modulus not appropriate for MSB test\n";
			ctFresh = NULL;
			return;	
		}

		LWE::CipherTextQ ctQ;
		LWE::CipherTextQ_phiN ctQ_phiN;

		MembershipTest(&ctQ_phiN, ACC, 1);
		LWE::KeySwitch(&ctQ, EK->KSkey, ctQ_phiN);
		LWE::ModSwitch(ctFresh, ctQ);	

	}



	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\\ ---- Accumulator ---- //////////////////////
	 *
	 */

	void Accumulator (LWE::CipherText * ctRaw, ct_FFT ACC, EvalKey * EK)
	{
		int vv = (ctRaw->b + q/(2*p)) % q;

		InitializeACC(ACC, vv);    		

		IncrLoop(ctRaw, ACC, EK);
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\\ ---- LWErefresh ---- //////////////////////
	 *
	 */

 	void LWErefresh (LWE::CipherText * ctRaw, LWE::CipherText ** ctFresh, EvalKey * EK, int * testArray)
	{
		ct_FFT ACC;
		
		Accumulator(ctRaw, ACC, EK);

		#pragma omp parallel sections firstprivate(ACC, EK)
		{
			#pragma omp section
			{
				Collapse(ACC, EK, ctFresh[0], testArray[0]);
			}
			#pragma omp section
			{
				Collapse(ACC, EK, ctFresh[1], testArray[1]);
			}			
			#pragma omp section
			{
				Collapse(ACC, EK, ctFresh[2], testArray[2]);
			}
			#pragma omp section
			{
				Collapse(ACC, EK, ctFresh[3], testArray[3]);
			}
			#pragma omp section
			{
				Collapse(ACC, EK, ctFresh[4], testArray[4]);
			}
		}
	}



	/*
	 *	\\\\\\\\\\\\\\\\\\\\\ ---- TESTING PART ---- ////////////////////
	 *
	 */

	int f_i (int m, int i)
	{
		if (m == i)
		{
			return 1;
		}
		return 0;
	}

	void correctnessTest (int m, int result, int i, int j)
	{

		if (f_i(m, i) != result)
		{
			std::cerr << "fatal error at the iteration " << j << "\n";

			std::cerr << "message = " << m << "\n";
			std::cerr << "i = " << i << "\n";
			std::cerr << "result = " << result << "\n";


			exit(0);
		}
		std::cerr << "Passed...\n";

	}

	void printError (int m, int i, LWE::CipherText * ct_to_test, LWE::SecretKey sk)
	{
		int er = ((ct_to_test->b - innerProdInt(sk, ct_to_test->a) - (f_i(m, i)*q/p))%q + q)%q;

		if (er > q/2)
		{
			er -= q;
		}

		std::cerr << "The error is " << er << "\n";
	}

 	void LWErefresh_correctnessTest (LWE::CipherText* ctRaw, LWE::CipherText ctFresh[p], EvalKey* EK, int * testArray, LWE::SecretKey sk, int m, int j)
	{

		ct_FFT ACC;
		
		int result, vv = (ctRaw->b + q/(2*p)) % q;

		InitializeACC(ACC, vv);    		

		IncrLoop(ctRaw, ACC, EK);
		
		#pragma omp parallel sections firstprivate(ACC, sk)
		{
			#pragma omp section
			{
				Collapse(ACC, EK, &(ctFresh[0]), testArray[0]);
				result = LWE::Decrypt_p (sk, ctFresh[0]);
				correctnessTest (m, result, 0, j);
				printError(m, 0, &(ctFresh[0]), sk);
			}
			#pragma omp section
			{
				Collapse(ACC, EK, &(ctFresh[1]), testArray[1]);
				result = LWE::Decrypt_p (sk, ctFresh[1]);
				correctnessTest (m, result, 1, j);
				printError(m, 1, &(ctFresh[1]), sk);
			}
			#pragma omp section
			{
				Collapse(ACC, EK, &(ctFresh[2]), testArray[2]);
				result = LWE::Decrypt_p (sk, ctFresh[2]);
				correctnessTest (m, result, 2, j);
				printError(m, 2, &(ctFresh[2]), sk);
			}
			#pragma omp section
			{
				Collapse(ACC, EK, &(ctFresh[3]), testArray[3]);
				result = LWE::Decrypt_p (sk, ctFresh[3]);
				correctnessTest (m, result, 3, j);
				printError(m, 3, &(ctFresh[3]), sk);
			}
			#pragma omp section
			{
				Collapse(ACC, EK, &(ctFresh[4]), testArray[4]);
				result = LWE::Decrypt_p (sk, ctFresh[4]);
				correctnessTest (m, result, 4, j);
				printError(m, 4, &(ctFresh[4]), sk);
			}
		}
	}	

}	// end of namespace HomACC 