#include <iostream>
#include <stdio.h>
#include <math.h>
#include "params.h"
#include "distrib.cpp"
#include "FFTRing.cpp"
#include "LWE.cpp"
#include "HomACC.cpp"


#include <sys/time.h>
#include <time.h>

/*
void noiseVecExtractor(ZmodQ * noiseVec, Ring_ModQ az, Ring_ModQ b_Acc, int * Y_m, int itn)
{
//	std::cerr << "\n noiseVecExtractor called!!\n";
	for (int i = 0; i < phi_N; ++i)
	{
		noiseVec[i] = b_Acc[i] - az[i] - (v*Y_m[i]);
//		std::cerr << "(" << b_Acc[i] << " , " << noiseVec[i] << ") , ";
	}
//	std::cerr << " ---end\n";
}
*/


long long int innerProd (Ring_ModQ a, Ring_ModQ b)
{
	long long int innp = 0;
	for (int i=0; i<phi_N; i++)
	{
		innp += ((long long int)a[i])*((long long int)b[i]);
	}
	return innp;
}


/*
void IncrLoop (LWE::CipherText* lweCTraw, HomACC::ct_FFT ACC, const HomACC::EvalKey* EK)
{
//	ZmodQ * noiseVec = (ZmodQ*) malloc(sizeof(ZmodQ)*phi_N);
//	long long int noise;
//	int * Y_m = (int*) malloc(sizeof(int)*phi_N);
//	int val = m, itn = 0;
//	Ring_ModQ z, az_temp, b_Acc;
//	Ring_FFT az_FFTtemp, z_FFT;
//	for (int i = 0; i < phi_N; ++i)
//	{
//		z[i] = HomACCsk[i];
//	}
//	FFTRing::Ring2FFT(&z_FFT, &z);


	for (int i = 0; i < n; ++i) 
	{
		int ai = (q - (lweCTraw->a)[i] % q) % q;

		for (int k = 0; k<d_r; ++k, ai /= B_r) 
		{
			int a0 = ai % B_r;

			if (a0) 
			{
//				itn++;
				/*
				if (itn > 10)
				{
					exit(0);
				}
*-/

				std::cerr << "loop iniciated\n iteration number " << itn << "\n";
				HomACC::ACCIncr_parallel_trigger(ACC, *(EK->BSkey[i][a0][k]));
//				std::cerr << " incr computed\n";
				//val += sk[i]*a0*Br_table[k];
//				std::cerr << " val computed = " << val << "\n";
				
			}
		}
	}
/*
	HomACC::Y_m_coeff(val, Y_m);

//				std::cerr << " y^m computed\n";

	// computing az
	for (int j = 0; j < FFTdim; ++j)
	{
		az_FFTtemp[j] = z_FFT[j]*ACC[1][0][j];
	}
//				std::cerr << " AZ iniciated\n";
//	FFTRing::FFT2Ring(&az_temp, &az_FFTtemp);
//				std::cerr << " az in FFT form\n";


//	FFTRing::FFT2Ring(&b_Acc, &ACC[1][1]);
//				std::cerr << " b extracted\n";

//	noiseVecExtractor(noiseVec, az_temp, b_Acc, Y_m, itn);

	for (int ll = 0; ll < p; ++ll)
	{
		noise = innerProd(noiseVec, HomACC::characFuncVector_real[ll]);
		std::cerr << "total noise with vector " << ll << " is "<<  noise << "\n";	
	}
*-/
	

	std::cerr << "done!\n";

}
*/



int minimum(int a, int b)
{
	if (a<b)
	{
		return a;
	}
	return b;
}

int f_i (int m, int i)
{
	if (m == i)
	{
		return 1;
	}
	return 0;
}


int main()
{
	std::cout << "Well... this is an experiment... let's see how it goes!\n";

	HomACC::Setup();										//This calls FFTRing::Setup() as well.

	std::cerr << "Generating LWE secret key \n";
	LWE::SecretKey sk;
	LWE::KeyGen(sk);
	
	//std::cerr << "Generating LWE-HomACC key \n";
	//LWE::SecretKey_phiN HomACCsk;
	//LWE::KeyGen_QphiN (HomACCsk);



	std::cerr << "Generating evaluation key \n";
	HomACC::EvalKey EK;
//	HomACC::KeyGen_ad(&EK, sk, &HomACCsk);
	HomACC::KeyGen(&EK, sk);


	int m;
	HomACC::ct_FFT ACC;
	LWE::CipherTextQ_phiN ct_QphiN;
	LWE::CipherText ct, ctFresh;
	LWE::CipherTextQ ctQ;
//	int result_QphiN, result_Q, result, noiseQ, noise, noiseCount;
	int result, noise, noiseCount=0;

	timespec time1, time2;


	clock_gettime(CLOCK_REALTIME, &time1);

	for (int j = 0; j < 100; ++j)
	{
			
		m = (rand()%p + p)%p;

		LWE::Encrypt_p (&ct, sk, m);

		int vv = (ct.b + q/(2*p)) % q;

		HomACC::InitializeACC(ACC, vv);
		//std::cerr << "Accumulator initialized!\n";

		//IncrLoop(&ct, ACC, &EK, sk, vv);
		HomACC::IncrLoop(&ct, ACC, &EK);


		//std::cerr << "Accumulator loop finished\n";

		std::cerr << "Iteration" << j+1 << ". The message is: " << m << "\n";
		for (int i=0; i<p; i++)
		{

			HomACC::MembershipTest(&ct_QphiN, ACC, i);
			//result_QphiN = LWE::Decrypt_QphiN(HomACCsk, ct_QphiN);
			//std::cerr << "the result ZQ for the " << i <<"th test is \t" << result_QphiN << "\n";

			LWE::KeySwitch(&ctQ, EK.KSkey, ct_QphiN);
			//result_Q = LWE::Decrypt_Q(sk, ctQ);

			//std::cerr << "the result Q for the " << i <<"th test is \t" << result_Q << "\n";

			LWE::ModSwitch(&ctFresh, ctQ);			

			result = LWE::Decrypt_p (sk, ctFresh);
			//std::cerr << "the final result for the " << i <<"th test is \t" << result << "\n";

			//noiseQ = ctQ.b-(m*Q/p);
			//std::cerr << "the noiseQ is \t\t" << noiseQ << "\n";
			noise = (ctFresh.b-((int)innerProd(sk,ctFresh.a)%q))%q - (result*(q/p))%q;
			noise = minimum(abs(noise), q - abs(noise));
			//std::cerr << "\t\t\t\t\touter r" << ((ctFresh.b-((int)innerProd(sk,ctFresh.a)%q))%q + q)%q << "\n";
			//std::cerr << "the noise is \t\t" << noise << "\n";

			if (noise > q/(4*p))
			{
				noiseCount++;
				std::cerr << "aaaah\n";
			}

			if (f_i(m, i) != result)
			{
				std::cerr << "fatal error at the iteration " << j << "\n noise = " << noise << "\n";

				exit(0);
			}
			std::cerr << "Passed...\n";

		}
	}

	std::cerr << "the number of errors is " << noiseCount << "\n";
	clock_gettime(CLOCK_REALTIME, &time2);

	uint64_t ttime = 1000000000*(time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec);

	std::cerr << "the total time is " << (long long unsigned int) ttime << "\n";



/*

	int u, v;
	HomACC::ct_FFT gswCTu, gswCTv;
	Ring_FFT HomACCsk_FFT;
	FFTRing::Ring2FFT(&HomACCsk_FFT, &HomACCsk);


	for (int i = 0; i < 20; ++i)
	{
		u = (rand()%q + q)%q;
		v = (rand()%q + q)%q;


		std::cerr << "the messages to encrypt are: \n \t\t\t u=" << u << "\n\t\t\t v=" << v << "\n";
		std::cerr << "the sum mod N is " << (u+v)%N << "\n";

		HomACC::Encrypt(gswCTu, HomACCsk_FFT, u);
		HomACC::Encrypt(gswCTv, HomACCsk_FFT, v);

		HomACC::ACCIncr_parallel_trigger(gswCTu, gswCTv);

		for (int l = 0; l < 5; ++l)
		{

			HomACC::MembershipTest(&ct_QphiN, gswCTu, l);

			rawResult = LWE::Decrypt_phiN(HomACCsk, ct_QphiN);

			std::cerr << "the result for " << l << " is " << rawResult << "\n";
		}
	}

*/
	return 0;
	
}



	
	
