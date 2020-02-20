#ifndef HOMACC_H
#define HOMACC_H

#include "LWE.h"

namespace HomACC {


	/* 
	 *	------ Type definitions ------ 
	 *	Types needed:
	 *		Ring_ModQ (defined in params.h)
	 *		Ring_FFT (defined in params.h)
	 *		
	 */

	typedef Ring_FFT ct_FFT[d_g2][2];					// Ciphertext in FFT form
	typedef ct_FFT* BootstrappingKey[n][B_r][d_r];
	typedef struct 
	{
		BootstrappingKey BSkey;
		LWE::SwitchingKey KSkey;
	} EvalKey;

	typedef Ring_ModQ ct_ModQ[d_g2][2];					// Ciphertext in coefficient form
	typedef Ring_ModQ dct_ModQ[d_g2][d_g2];				// Decomposed Ciphertext in coeff form
	typedef Ring_FFT  dct_FFT[d_g2][d_g2];				// Decomposed Ciphertext in FFT form

	/* ---------- Functions ---------- */

	void Setup ();

	void KeyGen (EvalKey *EK, const LWE::SecretKey LWEsk);

	void Accumulator (LWE::CipherText * ctRaw, ct_FFT ACC, EvalKey * EK);
	void MSBTest (ct_FFT ACC, EvalKey* EK, LWE::CipherText * ctFresh);
	void Collapse (ct_FFT ACC, EvalKey* EK, LWE::CipherText * ctFresh, int l);



	void Encrypt (ct_FFT ct, Ring_FFT sk_FFT, int m);
	void LWErefresh (LWE::CipherText* lweCTraw, LWE::CipherText* lweCTfresh, const EvalKey* EK, int t=-1);
	void LWErefresh_correctnessTest (LWE::CipherText* ctRaw, LWE::CipherText ctFresh[p], EvalKey* EK, int * testArray, LWE::SecretKey sk, int m, int j);

}

#endif