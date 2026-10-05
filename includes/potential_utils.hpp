#ifndef POTENTIAL_UTILS_HPP
#define POTENTIAL_UTILS_HPP

#include "mfem.hpp"

mfem::GridFunction Compute_Kaps(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2);
mfem::GridFunction Compute_i0(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2);
mfem::GridFunction Compute_OCV(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2);
mfem::GridFunction Compute_Kpl(mfem::GridFunction &De_gf, 
							   mfem::GridFunction &Cn_gf,
							   double scaleConst);
#endif