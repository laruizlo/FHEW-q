#include <iostream>
#include "params.h"
#include "LWE.h"
#include "HomACC.h"
#include <cassert>
#include <cstdlib>
  

namespace FHEtools
{

  /*
    \\\\\\\\\\\\\\\\\\\\\\\\ ----- fwrite_ek ----- ////////////////////////

  */

  void fwrite_ek(const HomACC::EvalKey& EK, FILE* f) 
  {
    /* --- Write bootstrapping key --- */
    for (int i = 0; i < n; ++i) 
    {
      for (int j = 1; j < B_r; ++j) 
      {
        for (int k = 0; k < d_r; ++k) 
        {
          assert(fwrite(EK.BSkey[i][j][k], sizeof(HomACC::ct_FFT), 1, f));
        }
      }
    }
    
    /* --- Write switching key --- */
    for (int i = 0; i < phi_N; ++i) 
    {
      for (int j = 0; j < B_ks; ++j) 
      {
        for (int k = 0; k < d_ks; ++k) 
        {
          assert(fwrite(EK.KSkey[i][j][k], sizeof(LWE::CipherTextQ), 1, f));
        }
      }
    }

  }


  /*
    \\\\\\\\\\\\\\\\\\\\\\\\\ ----- fread_ek ----- ////////////////////////

  */

  HomACC::EvalKey* fread_ek(FILE* f) 
  {

    HomACC::EvalKey* EK = new HomACC::EvalKey;
    
    // Read bootstrapping key
    for (int i = 0; i < n; ++i) 
    {
      for (int j = 1; j < B_r; ++j) 
      {
        for (int k = 0; k < d_r; ++k) 
        {
          EK->BSkey[i][j][k] = (HomACC::ct_FFT*) fftw_malloc(sizeof(HomACC::ct_FFT));
          assert(fread(EK->BSkey[i][j][k], sizeof(HomACC::ct_FFT), 1, f));
        }
      }
    }
    std::cout << "BSKey Read. \n";
    
    // Read switching key
    for (int i = 0; i < phi_N; ++i) 
    {
      for (int j = 0; j < B_ks; ++j) 
      {
        for (int k = 0; k < d_ks; ++k) 
        {
          EK->KSkey[i][j][k] = new LWE::CipherTextQ;
          assert(fread(EK->KSkey[i][j][k], sizeof(LWE::CipherTextQ), 1, f));
        }
      }
    }
    std::cout << "KSkey Read. \n";


    return EK;
  }

}