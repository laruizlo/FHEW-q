#include <iostream>
#include <cstdlib>
#include "LWE.h"
#include "HomACC.h"
#include "distrib.h"
#include "time_profiler.h"
#include "FHEtools.h"
#include <stdio.h>




#include "FFTRing.h"
#include <complex.h>
#include <time.h>
#include <omp.h>




int main()
{
	







	std::cout << "Initializing...\n";

	srand(round(clock()));
	HomACC::Setup();

	std::cout << "Setup Ready Commander... \n";


	LWE::SecretKey LWEsk;
	HomACC::EvalKey EK;


	std::cout << "Generating LWE Key ... \n";
	LWE::KeyGen(LWEsk);
	
	std::cout << "LWE Key Generated ...\n";


	std::cout << "Generating Evaluation Key ... \n";
	HomACC::KeyGen(&EK, LWEsk);
	
	std::cout << "Evaluation Key Generated ...\n";



	//FILE* ekFile;
	/*	  	
	ekFile = fopen("ek.txt", "w");
	FHEtools::fwrite_ek(EK, ekFile);
	fclose(ekFile);
	std::cout << "Exiting!\n";  
	exit(0);
	*/

	/*
	std::cout << "Reading evaluation key ... \n\n";
	ekFile = fopen("ek.txt", "r");
	EK = *FHEtools::fread_ek(ekFile);
	fclose(ekFile);
	std::cout << "Done!\n\n";
	*/



	int b1, b2;

	LWE::CipherText ct_b1, ct_b2, ct_sum_raw, ct_sum_fresh;


	int t1, t2;

	t1 = clock();
	for (int i=0; i<10; i++)
	{
		b1 = (int) round(5*(rand()/(RAND_MAX+1.0)));
		std::cout << "b1 " << b1%5 << "\n";

		b2 = (int) round(5*(rand()/(RAND_MAX+1.0)));
		std::cout << "b2 " << b2%5 << "\n";

		LWE::Encrypt(&ct_b1, LWEsk, b1);
		LWE::Encrypt(&ct_b2, LWEsk, b2);

		std::cout << "b1? " << LWE::Decrypt(LWEsk, ct_b1) << "\n";
		std::cout << "b2? " << LWE::Decrypt(LWEsk, ct_b2) << "\n";



		std::cout << "\n sum = " << b1 + b2 << "\n";



		for (int j=0; j<n; j++)
		{
			//std::cout << ((ct_b1.a[j] + ct_b2.a[j]) + q) % q << "\n";
			ct_sum_raw.a[j] = ((ct_b1.a[j] + ct_b2.a[j])%q + q) % q;
		}
		ct_sum_raw.b = ((ct_b1.b + ct_b2.b)%q + q) % q;

		HomACC::LWErefresh(&ct_sum_raw, &ct_sum_fresh, &EK, 0);

		int b_prime = LWE::Decrypt(LWEsk, ct_sum_fresh);

		std::cout << "b_prime " <<  b_prime << "\n\n\n";
	}
	t2 = clock();
	std::cout << (float) (t2-t1)/CLOCKS_PER_SEC << "\n\n";






	std::cerr << "experiment begins \n ";

	for (int i=0; i<100000; i++)
	{

		b1 = (int) round(5*(rand()/(RAND_MAX+1.0)));

		LWE::Encrypt(&ct_b1, LWEsk, b1);

		int b1_d = LWE::Decrypt(LWEsk, ct_b1);

		if (b1_d != b1)
		{
			std::cerr << b1 << ":";
			std::cerr << b1_d << "\n";
		}
	}

	std::cout << "done\n\n";


	/*
		freeing the ciphertext pointers
	*/			
	FFTRing::Destructor();
}

