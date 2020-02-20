#include <iostream>
#include <stdio.h>
#include <math.h>
#include "params.h"
#include "distrib.cpp"
#include "FFTRing.cpp"
#include "LWE.cpp"
#include "HomACC.cpp"
#include "gates.cpp"
#include <omp.h>


#include <sys/time.h>
#include <time.h>


LWE::SecretKey sk;
HomACC::EvalKey EK;
int keySet = 1;


void GenSetup ()
{

	HomACC::Setup();										//This calls FFTRing::Setup() as well.

	std::cerr << "Generating LWE secret key \n";
	LWE::KeyGen(sk);

	std::cerr << "Generating evaluation key \n";
	HomACC::KeyGen(&EK, sk);

	keySet = 0;
}

void test (int numItns)
{
	
	std::cout << "Well... this is an experiment... let's see how it goes!\n";

	int testArray[p] = {0,1,2,3,4};					// To be changed

	int m;
	LWE::CipherText ct, ctFresh[p];
	GenSetup();

	timespec time1, time2;
	clock_gettime(CLOCK_REALTIME, &time1);

	for (int j = 0; j < numItns; ++j)
	{
			
		m = (rand()%p + p)%p;

		LWE::Encrypt_p (&ct, sk, m);

		std::cerr << "Iteration" << j+1 << ". The message is: " << m << "\n";

		HomACC::LWErefresh_correctnessTest(&ct, ctFresh, &EK, testArray, sk, m, j);
	}

	clock_gettime(CLOCK_REALTIME, &time2);

	uint64_t ttime = 1000000000*(time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec);

	std::cerr << "the total time is " << (long long unsigned int) ttime << "\n";




	exit(0);


}



void sum (int a, int b)
{
	if (keySet)
	{
		GenSetup();		
	}


	int bitlength = (max ( log2(a), log2(b) ) )+1;

	std::cout << bitlength << "-bit integers\n";

	int a_decomp [bitlength], b_decomp [bitlength], s_decomp [bitlength+1], s=0;
	LWE::CipherText ct_a [bitlength], ct_b [bitlength], ct_s [bitlength+1], ct_ci, ct_co;

	for (int i = 0; i < bitlength; ++i)
	{
		a_decomp[i] = a % 2;
		b_decomp[i] = b % 2;

		LWE::Encrypt_p(&ct_a[i], sk, a_decomp[i]);
		LWE::Encrypt_p(&ct_b[i], sk, b_decomp[i]);

		a /= 2;
		b /= 2;
	}
	std::cerr << "Bits encrypted!\n";

	LWE::Encrypt_p(&ct_ci, sk, 0);
	for (int i = 0; i < bitlength; ++i)
	{
		Gates::FullAdder_1 (&ct_a[i], &ct_b[i], &ct_ci, &ct_s[i], &ct_co, &EK);
	}
	ct_s[bitlength] = ct_co;

	for (int i = 0; i < bitlength+1; ++i)
	{
		s_decomp[i] = LWE::Decrypt_p (sk, ct_s[i]);
		s += (int) s_decomp[i]<<i;
	}
	std::cerr << "The result is " << s << "\n";

}



int main (int argc, char **argv) 
{

	char opt = 'a';

	if (argc == 2)
	{
		opt = *(argv[1]);
	}

	startingPoint:

		if (opt == 't')
		{
			int numItns;
			std::cout << "What is the number of times that you want to repeat the experiment?\n";
			std::cin >> numItns;
			test (numItns);
		}
		else if (opt == 's')
		{
			int a, b;
			std::cout << "What is the first number? \n";
			std::cin >> a;
			std::cout << "What is the second number? \n";
			std::cin >> b;

			sum (a, b);
		}
		else if (opt == 'x')
		{
			exit(0);
		}
		std::cerr << "\nusage: " << argv[0] << "\n\n"
		<< "\t-Type t to run a test\n"
		<< "\t-Type s to test the full adder\n"
		<< "\t-Type h for a detailed description\n"
		<< "\t-Type x exit\n";

		std::cin >> opt;

	goto startingPoint;

	return 0;
	 
}



	
	
