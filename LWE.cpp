#include <cstdlib>
#include <iostream>
#include "LWE.h"
#include "params.h"
#include <cassert>


#include "time_profiler.h"

using namespace std;
namespace LWE 
{

	/* 
     *	\\\\\\\\\\\\\\\\\---- Secret-Key Generator ----///////////////////
     *
	 *	The secret key is a vector of n entries chosen from Z_q either 
	 *	uniformly at random or as a short vector. For this implementation, 
	 *	the secret key used in the inner encryption scheme is a short vector,
	 *	whereas for the outer encryption scheme is a random vector with 
	 *	slightly larger entries.
	 *
	 *	In order to guarantee that the vector is short, the following method 
	 *	samples the entries from the distribution Chi_Binary (which takes 
	 *	its values on the set {-1, 0, 1}), and limits the total size (i.e. 
	 *	the sum) of the entries to 5 and the absolute size (i.e. the sum of 
	 *	the absolute values) to 5 + n/2.
	 *	
	 *	The second method just samples the entries from the distribution 
	 *	Chi1.
	 */

	void KeyGen(SecretKey sk) 
	{
    	/*
		 *	For some reason this function is outputting the same secret key 
		 *	every time.
		 */

		int s, ss;

		KeyGenRestart:
		s=0;
		ss=0;

		for (int i = 0; i < n; ++i) 
		{
			sk[i] = Sample(Chi_Binary);
			s+= sk[i];
			ss+= abs(sk[i]);
		}
    	

		if (abs(s)>5) goto KeyGenRestart;
		if (abs(ss - n/2)>5) goto KeyGenRestart;
	}


	void KeyGen_QphiN(SecretKey_phiN sk) 
	{
		for (int i = 0; i < phi_N; i++) 
		{
			sk[i].value = Sample(Chi1);
		}
	}

	
	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\---- Encrypt ----////////////////////////
	 *
	 *	Given a secret key s \in Z_q^n and a message m \in Z_p, the  
	 *	encryption algorithm for the LWE cryptosystem is given by
	 *
	 *			a 			<-	Z_q^n 	(uniformly at random)
	 *			e 			<-	Z_q 	(Sampled from a (sub)gaussian dist)
	 *			LWE(s, m) 	= 	(a, round ((a\cdot s) + mq/t + e))
	 *
	 *	
	 *
	 *
	 *
	 */

	void Encrypt_p(CipherText* ct, const SecretKey sk, int m) 
	{
		double b;
		m = m % (p);										// Sanitizing the message 

		b =  m * (double)q / p;							// Initializing b

		int e = Sample(Chi3);
		b += (double)e;

		for (int i = 0; i < n; i++) 
		{
			ct->a[i] = rand() % q;
			b += (ct->a[i] * sk[i]) % q;
		}
		ct->b = (((int)round(b) % q) + q) % q;
	}


	void Encrypt_p_prime(CipherText* ct, const SecretKey sk, int m) 
	{
		double b;
		m = m % (p_prime);										// Sanitizing the message 

		b =  m * (double)q / p_prime;							// Initializing b

		int e = Sample(Chi3);
		b += (double)e;

		for (int i = 0; i < n; i++) 
		{
			ct->a[i] = rand() % q;
			b += (ct->a[i] * sk[i]) % q;
		}
		ct->b = (((int)round(b) % q) + q) % q;
	}



	/*
	 *	\\\\\\\\\\\\\\\\\\\\---- Noiseless-Encrypt ----//////////////////////
	 *
	 *	
	 *	Use for debugging purposes 
	 */

	void NoiselessEncrypt(CipherText* ct, const SecretKey sk, int m) 
	{
		double b;
		m = m % (2*p);										// Sanitizing the message 

		b =  m * (double)q / (2*p);							// Initializing b

		for (int i = 0; i < n; i++) 
		{
			ct->a[i] = rand() % q;
			b += (ct->a[i] * sk[i]) % q;
		}
		ct->b = (((int)round(b) % q) + q) % q;
	}

	/*
	 *	\\\\\\\\\\\\\\\\\\\\---- Encrypt_phiN ----//////////////////////
	 *
	 *
	 *
	 *
	 */


	void Encrypt_QphiN(CipherTextQ_phiN* ct, const SecretKey_phiN sk, int m) 
	{
		ZmodQ b;											// Signed reduction modulo 2^32 is automatic
		m = m % p;											// Sanitizing the message 

		b.value =  m * (Q / p);									// Initializing b, lets see if the error matters 

		int e = Sample(Chi3);								// Sampling the error

		b.value += (double)e;										// Adding the error to b

		for (int i = 0; i < phi_N; i++) 
		{
			(ct->a[i]).value = rand();
			b.value += (ct->a[i]).value * (sk[i]).value;
		}
		ct->b = b;
	}	

	/*
	 *	\\\\\\\\\\\\\\\\\\\\\\---- Decrypt_p ----////////////////////////
	 *
	 *
	 *
	 *
	 */

	int Decrypt_p (const SecretKey sk, const CipherText& ct) 
	{
		int r = ct.b;

		for (int i = 0; i < n; i++)
		{ 
			r -= ct.a[i] * sk[i];
		}
		//std::cerr << "\t\t\t\tinner r " << (r%q + q)%q << "\n";
		r = ((r % q) + q + q/(2*p)) % q;

		return p*r/q;    
	}
	

	/*
	 *	\\\\\\\\\\\\\\\\\\\\\---- Decrypt_p_prime ----/////////////////////
	 *
	 *
	 *
	 *
	 */

	int Decrypt_p_prime(const SecretKey sk, const CipherText& ct) 
	{
		int r = ct.b;

		for (int i = 0; i < n; i++)
		{ 
			r -= ct.a[i] * sk[i];
		}
		r = ((r % q) + q) % q;

		return p_prime*r/q;    
	}


	int DecryptDetail(const SecretKey sk, const CipherText& ct) 
	{
		int r = ct.b;

		for (int i = 0; i < n; ++i) 
		{
			r -= ct.a[i] * sk[i];
		}

		r = ((r % q) + q + q/8) % q;
		int m = 4*r/q;

		cerr << "\t Value " << r - q/8 << "\t Decoded as " << m << " * " << q/4 << "\t + " << r - m * q/4 - q/8<< endl;
		cout << r - m * q/4 - q/8 << ", "; 

		return m;
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\---- Decrypt_QphiN ----//////////////////////
	 *
	 *	
	 *
	 *
	 */

	int Decrypt_QphiN(const SecretKey_phiN sk, const CipherTextQ_phiN& ct) 
	{
		uZmodQ r = {ct.b.value};

		for (int i = 0; i < phi_N; i++)
		{ 
			r.value -= ct.a[i].value * (sk[i]).value;
		}

		int m = round(p*((double)r.value/Q));								// Message

		return (m%p);
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\---- Decrypt_Q ----//////////////////////
	 *
	 *	
	 *
	 *
	 */

	int Decrypt_Q(const SecretKey sk, const CipherTextQ ctQ)
	{
		uZmodQ r = {ctQ.b.value};

		for (int i = 0; i < n; i++)
		{ 
			r.value -= ctQ.a[i].value * sk[i];
		}

		int m = round(p*((double)r.value/Q));								// Message

		return (m%p);	
	}

	/*
	 *	\\\\\\\\\\\\\\\\\\\\---- round_qQ ----//////////////////////
	 *
	 *	
	 *
	 *
	 */	

	int round_qQ(ZmodQ v) 
	{
		return (int)(floor(.5 + (double) v.value * (double) q / (double) Q) + q) % q;
	}


	/*
	 *	\\\\\\\\\\\\\\\\\\\\---- ModSwitch ----//////////////////////
	 *
	 *	
	 *
	 *
	 */

	void ModSwitch(CipherText *ct, const CipherTextQ &c) 
	{
		for (int i = 0; i < n; ++i) 
		{
			ct->a[i] = round_qQ(c.a[i]);  
		}

		ct->b = round_qQ(c.b);
	}
  
	/* 
	 * \\\\\\\\\\\\\\\\\\\\\\\ ----- SwitchingKey ----- ///////////////////
	 *
	 * The LWE key-switch key is a 3-dimensional array of pointers to LWE-
	 * ciphertexts over the Q-space (Z_Q). It takes a key (n-key) and an N-
	 * key as input, along with an array to write the result. The (i, j, k) 
	 * entry of the arrays is of the form
	 *         [(a_0,..., a_n), sum_{l=0}^n a_l s_l - z_i*j*KS^k + e] 
	 * which is an LWE encryption of z_i*j*KS.
	 * KS (KS_base) is the base on which all the possible entries are 
	 * written, KS_exp is the maximum exponent of the base, therefore it is 
	 * able to describe numbers from 0 to KS_base^KS_exp
	 * NOTE: there is no explicit reduction mod Q.
	 */

	void SwitchingKeyGen(SwitchingKey res, const SecretKey new_sk, const SecretKey_phiN old_sk) 
	{

		for (int i = 0; i < phi_N; i++) 
		{
			for (int j = 0; j < B_ks; j++) 
			{
				for (int k = 0; k < d_ks; k++) 
				{
					CipherTextQ* ct = new CipherTextQ;					// Declare the format of the corresponding entry.
					(ct->b).value = (-(old_sk[i]).value) * j * (KS_table[k].value)  +  Sample(Chi2);	// b is now equal to z_i*j*KS^k + e
					for (int l = 0; l < n; ++l) 
					{
						(ct->a[l]).value = rand();								// Sample a random vector entry by entry
						(ct->b).value += (ct->a[l]).value * new_sk[l];					// add the interior product to b entry by entry
					}
					res[i][j][k] = ct;									// Write the result
				}
			}
		}
	}




	void KeySwitch(CipherTextQ* res, const SwitchingKey K, const CipherTextQ_phiN& ct) 
	{

		for (int k = 0; k < n; k++) 
		{
			(res->a[k]).value = 0;
		}
		res->b = ct.b;
    
		for (int i = 0; i < phi_N; i++) 
		{
			uZmodQ a = {-ct.a[i].value};
			for (int j = 0; j < d_ks; j++, a.value /= B_ks) 
			{
				uZmodQ a0 = {a.value % B_ks};
				for (int k = 0; k < n; k++)
				{
					(res->a[k]).value -= ((K[i][a0.value][j])->a[k]).value;
				}
				(res->b).value -= ((K[i][a0.value][j])->b).value;
			}
		} 
 
	}

}
