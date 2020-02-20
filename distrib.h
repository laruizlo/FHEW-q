#ifndef DISTRIB_H
#define DISTRIB_H

#include "params.h"

/* --- Distribution Structure --- */

typedef struct {
  double std_dev;                                         // Standard deviation
  int offset;                                             // Size of table. Set to 0 to ignore table and use rejection sampling
  const float* table;                                     // CDF of Gaussian of standard deviation std_dev 
} Distrib;

int Sample(const Distrib& Chi);                           // sample integer with gaussian distribution

/* ----- Distributional Tables ----- */


const float Chi1_Table[12] = 
{
  6.82100340236e-20, 6.82745420379e-17, 3.70723486968e-14,
  1.09385320287e-11, 1.75776798627e-09, 1.54306108374e-07,
  7.43151128657e-06, 0.000197571144374, 0.00292649511724,
  0.024505041784, 0.11878103083, 0.34690097998
};




const float Binary_Table[1] = {
 .25 
};

const float NoTable[1] = {1};


/* ----- Distributions ----- */

const Distrib Chi1 = {
  1.27,                                   // std_dev, in this case it won't be used!
  12,                                     // offset
  Chi1_Table                              // table
};

const Distrib Chi2 = {
  (double) (1 << 17),                     // std_dev
  0,                                      // offset
  NoTable                                 // table
};

const Distrib Chi3 = {
  (double) 2.9,                           // std_dev
  0,                                      // offset
  NoTable                                 // table
};

const Distrib Chi_Binary = {
  0,                                      // std_dev
  1,                                      // offset
  Binary_Table
};




#endif
