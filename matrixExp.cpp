#include <iostream>
#include <stdio.h>
#include <math.h>
#include "params.h"
#include "distrib.cpp"
#include "FFTRing.cpp"
#include "LWE.cpp"
#include "HomACC.cpp"

/*
void mat_mult(int ** m_1, int ** m_2)
{

}
*/

void printVector (Ring_ModQ vec, char name)
{

	std::cerr << name << " = (" << vec[0];					// printing
	for (int i = 1; i < phi_N; ++i)
	{
		std::cerr << "," << vec[i];
	}
	std::cerr << ")\n";

}

void Poly_Mult(Ring_ModQ *a, Ring_ModQ *b, Ring_ModQ *c)
{
	std::cerr << "multiplying\n";


	Ring_FFT a_FFT, b_FFT, c_FFT;

	FFTRing::Ring2FFT(&a_FFT, a, 0);
	FFTRing::Ring2FFT(&b_FFT, b, 0);

	for (int i = 0; i < phi_N; ++i)
	{
		c_FFT[i] = ((double complex)a_FFT[i]) * ((double complex)b_FFT[i]);
	}

	FFTRing::FFT2Ring(c, &c_FFT);


}


void Matrix_Mult(int**a, int*b, int*c)
{
	for (int i = 0; i < phi_N; ++i)
	{
		c[i] = 0;
		for (int j = 0; j < phi_N; ++j)
		{
			c[i] += a[i][j]*b[j];
		}
	}
}

int main()
{

	std::cerr << "size " << sizeof(long int) << "\n";


	FFTRing::Setup();


	std::cerr << "multi experiment\n\n";

	int dim = phi_N;

	Ring_ModQ a, b, c;
	int ** a_arrow, * bprime, * cprime;


	for (int i = 0; i < phi_N; ++i)
	{
		a[i] = ((ZmodQ) rand())%(1<<11);
		b[i] = ((ZmodQ) rand())%(1<<11);
	}

	Poly_Mult(&a, &b, &c);



	a_arrow = (int**) malloc(sizeof(int*)*dim);
	bprime = (int*) malloc(sizeof(int) * dim);
	cprime = (int*) malloc(sizeof(int) * dim);

	for (int i = 0; i < dim; ++i)
	{
		a_arrow[i] = (int*) malloc(sizeof(int)*dim);
		bprime[i] = b[i];

	}

	HomACC::vector2matrix(a, a_arrow);


	printVector(a, 'a');
	printVector(b, 'b');
	printVector(c, 'c');


	Matrix_Mult(a_arrow, bprime, cprime);

	printVector(cprime, 'l');


	for (int i = 0; i < phi_N; ++i)
	{
		if (cprime[i] =! c[i])
		{
			std::cerr << "noooooo\n\n\n";
		}
	}



/*
	HomACC::vec_x_mat(a, b_arrow)
*/


}