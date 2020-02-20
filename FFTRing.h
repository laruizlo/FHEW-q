#ifndef FFTRING_H
#define FFTRING_H

#include "params.h"

namespace FFTRing
{
	void Setup ();
	void Ring2FFT (Ring_FFT *coVector, const Ring_ModQ *reVector, const int index=0);
	void FFT2Ring (Ring_ModQ *reVector, const Ring_FFT *coVector, const int index=0);
	void Destructor ();
}

#endif