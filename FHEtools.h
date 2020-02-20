#ifndef FHEWTOOLS_H
#define FHEWTOOLS_H

#include "HomACC.h"
#include <cassert>
#include <cstdlib>


namespace FHEtools
{
  void fwrite_ek(const HomACC::EvalKey& EK, FILE* f);
  HomACC::EvalKey* fread_ek(FILE* f);
}

#endif