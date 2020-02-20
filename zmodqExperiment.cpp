#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <stdint.h>


const int d_Q = 6;

struct ZmodQ
{
	int64_t	value : d_Q;
};


int main() 
{
	ZmodQ a = {1}, b = {1};


	for (int i = 0; i < 1000; ++i)
	{
		//std::cerr << "(" << (a.value) << ";" << i << ")";
		std::cerr << "(" << (b.value) << ";" << i << ")";
		a.value += 1;
		b.value = b.value << 1;
		a = b;
		std::cerr << "(" << (a.value) << ";" << i << ")";

	}	
	std::cerr << "\n";



	return 0;
}
