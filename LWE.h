#ifndef LWE_H
#define LWE_H

#include "params.h"
#include "distrib.h"



namespace LWE 
{
  
	typedef struct 
	{                                                                 // n-dimensional LWE ciphertext modulo q
		int a[n];                                                       // array of n integers modulo q 
		int b;                                                          // integer modulo q
	} CipherText;

	typedef struct 
	{                                                                 // n-dimensional LWE ciphertext modulo Q
		ZmodQ a[n];
		ZmodQ b;
	} CipherTextQ;
  
	typedef struct 
	{                                                                 // N-dimensional LWE ciphertext modulo Q
		ZmodQ a[phi_N];
		ZmodQ b;
	} CipherTextQ_phiN;
  
	typedef int SecretKey[n];                                         // n-dimensional LWE secret key 
	typedef ZmodQ SecretKey_phiN[phi_N];                                        // N-dimensional LWE secret key (Why is this here?) 
	typedef CipherTextQ* SwitchingKey[phi_N][B_ks][d_ks];  

	/* ---------- Functions ---------- */

	void KeyGen(SecretKey sk);
	void KeyGen_QphiN(SecretKey_phiN sk);
	void Encrypt_p(CipherText* ct, const SecretKey sk, int m);
	void Encrypt_p_prime(CipherText* ct, const SecretKey sk, int m);
	int Decrypt_p(const SecretKey sk, const CipherText& ct);
	int DecryptDetail(const SecretKey sk, const CipherText& ct);      //Uncomment this line for debugging purposes
  
	/* 
	 *	Generate key material (SwitchingKey) required by KeySwitch to 
	 *	transform LWE encryptions under old_sk into LWE encryptions under 
	 *	new_sk
	 */

	void SwitchingKeyGen(SwitchingKey res, const SecretKey new_sk, const SecretKey_phiN old_sk);
	void KeySwitch(CipherTextQ* res, const SwitchingKey K, const CipherTextQ_phiN& ct);

	// Changes an LWE ciphertext modulo Q into an LWE ciphertext modulo q
	void ModSwitch(CipherText* ct, const CipherTextQ& c); 
	int round_qQ(ZmodQ v);

}




#endif
