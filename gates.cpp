#include <iostream>
#include <cstdlib>
#include <omp.h>
#include "LWE.h"
#include "HomACC.h"
 
#include <complex.h>
#include <fftw3.h>



namespace Gates
{


	void sum_raw (LWE::CipherText * ct_r, LWE::CipherText * ct_0, LWE::CipherText * ct_1)
	{
		for (int i = 0; i < n; ++i)
		{
			ct_r->a[i] = (ct_0->a[i] + ct_1->a[i]) % q;
		}
		ct_r->b =  (ct_0->b + ct_1->b) % q;
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\\ ---- Compare ---- //////////////////////
	 *
	 *	TODO: Test this function.
	 */

	void Compare (LWE::CipherText * ct_r, LWE::CipherText * ct_0, LWE::CipherText * ct_1, HomACC::EvalKey* EK)
	{
		LWE::CipherText * ct_diff = new LWE::CipherText;
		HomACC::ct_FFT ACC;

		for (int i = 0; i < n; ++i)
		{
			ct_diff->a[i] = (ct_0->a[i] - ct_1->a[i]) % q;
		}
		ct_diff->b = (ct_0->b - ct_1->b) % q;

		HomACC::Accumulator(ct_diff, ACC, EK);

		// Collapse here to MSB
		HomACC::MSBTest(ACC, EK, ct_r);
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\\ ---- Ordering ---- //////////////////////
	 *
	 *	TODO: Test this function.
	 */


	 /*
	void Ordering ()
	{

	}
	*/


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\\ ---- FullAdder_1 ---- //////////////////////
	 *
	 *	TODO: Test this function.
	 */

	void FullAdder_1 (LWE::CipherText * ct_b1, LWE::CipherText * ct_b2, LWE::CipherText * ct_ci, LWE::CipherText * ct_s, LWE::CipherText * ct_co, HomACC::EvalKey * EK)
	{
		if (p < 4)
		{
			std::cerr << "Message modulus not appropriate for the 1-bit full-adder\n";
			ct_co = NULL;
			ct_s = NULL;
			return;
		}

		LWE::CipherText ct_sum, ct_sum_1, ct_sum_2, ct_sum_3;
		HomACC::ct_FFT ACC;

		sum_raw (&ct_sum, ct_b1, ct_b2);
		sum_raw (&ct_sum, &ct_sum, ct_ci);

		HomACC::Accumulator(&ct_sum, ACC, EK);

		//Not sure of this... I think I need different test vectors
		HomACC::Collapse (ACC, EK, &ct_sum_1, 1);				
		HomACC::Collapse (ACC, EK, &ct_sum_2, 2);
		HomACC::Collapse (ACC, EK, &ct_sum_3, 3);				

		sum_raw (ct_s, &ct_sum_1, &ct_sum_3);
		sum_raw (ct_co, &ct_sum_2, &ct_sum_3);
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\\ ---- FullAdder_k ---- //////////////////////
	 *
	 *	TODO:
	 *		- Implement LWE::CipherText_k
	 *		- Test this function.
	 */

	 /*
	void FullAdder_k (LWE::CipherText * ct_b1, LWE::CipherText * ct_b2, LWE::CipherText * ct_ci, LWE::CipherText * ct_s, LWE::CipherText * ct_co, HomACC::EvalKey * EK, int k)
	{
		if (p < (int)1 << k)
		{
			std::cerr << "Message modulus not appropriate for the " << k <<"-bit full-adder\n";
			ct_co = NULL;
			ct_s = NULL;
			return;
		}

		LWE::CipherText_k * ct_sum, * ct_sum_1, * ct_sum_2, * ct_sum_3;
		ct_FFT ACC;

		sum_raw (ct_sum, ct_b1, ct_b2);
		sum_raw (ct_sum, ct_sum, ct_ci);

		HomACC::Accumulator(ct_sum, ACC, EK);

		//Not sure of this... I think I need different test vectors
		HomACC::Collapse (ACC, EK, ct_sum_1, 1);				
		HomACC::Collapse (ACC, EK, ct_sum_2, 2);
		HomACC::Collapse (ACC, EK, ct_sum_3, 3);				

		sum_raw (ct_s, ct_sum_1, ct_sum_3);
		sum_raw (ct_co, ct_sum_2, ct_sum_3);
	}
	*/

}
