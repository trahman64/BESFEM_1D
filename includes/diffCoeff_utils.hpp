#ifndef DIFFCOEFF_UTILS_HPP
#define DIFFCOEFF_UTILS_HPP

#include "mfem.hpp"

mfem::GridFunction Compute_Damb(mfem::GridFunction &Cn);
mfem::GridFunction Compute_DLi(mfem::GridFunction &Cn);
mfem::GridFunction Compute_Dmp(mfem::GridFunction &De_gf, double tc1);
void ComputeConstantce(double t_minus, double Cst1,
                       double &tc1, double &tc2, double &scaleConst);

#endif