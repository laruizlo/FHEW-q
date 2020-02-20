#ifndef GATES_H
#define GATES_H

#include "LWE.h"
#include "HomACC.h"

namespace Gates
{
	void Compare(LWE::CipherText ct_r, LWE::CipherText_4 ct_0, LWE::CipherText_4 ct_1, EvalKey* EK);

	//void Ordering();
	
	void FullAdder_1(LWE::CipherText * ct_b1, LWE::CipherText * ct_b2, LWE::CipherText * ct_ci, LWE::CipherText * ct_s, LWE::CipherText * ct_co, HomACC::EvalKey * EK);

	void FullAdder_k(LWE::CipherText * ct_b1, LWE::CipherText * ct_b2, LWE::CipherText * ct_ci, LWE::CipherText * ct_s, LWE::CipherText * ct_co, HomACC::EvalKey * EK);

}

#endif 