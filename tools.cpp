#include "tools.h"
#include <iostream>


namespace tools
{

	ZmodQ inv_modQ(ZmodQ val)
	{
		int64_t a = val.value, b = Q;
		int64_t b0 = b, t, quot;
		int64_t x0 = 0, x1 = 1;
		if (b == 1) return (ZmodQ) {1};
		while (a > 1) 
		{
			quot = a / b;
			t = b, b = a % b, a = t;
			t = x0, x0 = x1 - quot * x0, x1 = t;
		}
		if (x1 < 0) x1 += b0;
		return ((ZmodQ) {x1});
	}
	

}
