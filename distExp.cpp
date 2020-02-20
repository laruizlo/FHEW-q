#include <iostream>
#include <stdio.h>
#include <math.h>
#include "distrib.cpp"

int main()
{

	Distrib chi = Chi1;

	srand(time(NULL));

/*
	std::cerr << "Chi1\n";
	for (int i=0; i<300; i++)
	{
		std::cerr << Sample(Chi1) << " ";
	}
	std::cerr << "\n\n";
	
	std::cerr << "Chi2\n";
	for (int i=0; i<300; i++)
	{
		std::cerr << Sample(Chi2) << " ";
	}
	std::cerr << "\n\n";
	
	std::cerr << "Chi3\n";
	for (int i=0; i<300; i++)
	{
		std::cerr << Sample(Chi3) << " ";
	}
	std::cerr << "\n\n";
	
	std::cerr << "Chi_Binary\n";
	for (int i=0; i<300; i++)
	{
		std::cerr << Sample(Chi_Binary) << " ";
	}
	std::cerr << "\n\n";
*/	



	int r, l=0, max=0;
	

	for (int i=0; i<100000; i++)
	{
		r = Sample(chi);

		if (abs(r) > l)
		{
			std::cerr << r << " ";
			l = abs(r);
		}
	}

	int lg = (int)log10(l);

	int pow10 = 1;
	int ref = 10;
	if (lg == 0) ref =1;
	for (int i=0; i<lg; i++)
	{
		pow10 *= 10;
	}
	std::cerr << "\nl = " << l << "\n";
	l /= round(pow10/ref);
	std::cerr << "l = " << l << "\n";
	l += chi.std_dev*10/pow10 + 4;

	int k[l];
	for (int i=0; i<l; i++)
	{
		k[i]=0;
	}

	std::cout << "Chi!\nstd dev = ";
	std::cout << chi.std_dev << "\n";

	for (int i=0; i<100000000; i++)
	{
		r = Sample(chi);

		for (int j=0; j<l; j++)
		{
			if (abs(round(ref*r/pow10)) == j)
			{
				k[j]++;
			}

		}

		if (abs(r) > max)
		{
			std::cerr << r << " (" << i << ") ";
			max = abs(r);
		}
	}
	k[0] *= 1+(10/ref);	

	std::cout << "\n the larger element found is " << max << "\n"; 		
	for (int i=0; i<l; i++)
	{
		std::cout << "The number of elements around " << i*pow10/ref << " is \t" << k[i] << "\n";
	}

}