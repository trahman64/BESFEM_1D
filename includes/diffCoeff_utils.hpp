#ifndef DIFFCOEFF_UTILS_HPP
#define DIFFCOEFF_UTILS_HPP

#include "mfem.hpp"

mfem::GridFunction ComputeDamb(mfem::GridFunction &Cn);
mfem::GridFunction ComputeDLi(mfem::GridFunction &Cn);

#endif