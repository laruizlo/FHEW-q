#include <iostream>
#include "distrib.cpp"
#include "LWE.cpp"

int main()
{
	int m, mp;
	LWE::CipherTextQ_phiN ct;
	LWE::SecretKey_phiN sk;

	LWE::KeyGen_phiN(sk);

	std::cout << "Key generated \n\n";

	for (int i=0; i<1000000; i++)
	{
		m = rand()%p;
		std::cerr << "m = " << m << "\n";
		LWE::Encrypt_phiN(&ct, sk, m);
 		mp = LWE::Decrypt_phiN(sk, ct);

 		if (m == mp)
 		{
			std::cerr << " = " << mp << "\n";
 		}
 		else
 		{
 			std::cerr << "Error with the message " << m << " decrypted as " << mp << " at the iteration " << i << "\n";
 			exit(0);
 		}


	}

}