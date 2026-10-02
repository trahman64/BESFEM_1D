#ifndef POTENTIAL_UTILS_HPP
#define POTENTIAL_UTILS_HPP

#include "mfem.hpp"

mfem::GridFunction ComputeKaps(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2);
mfem::GridFunction Compute_i0(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2);
mfem::GridFunction Compute_OCV(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2);

#endif