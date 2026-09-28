#include "mfem.hpp"
#include <cmath>


mfem::GridFunction ComputeKaps(mfem::GridFunction &Cn)
{
    mfem::GridFunction Kps(Cn.FESpace());
    for (int i = 0; i < Cn.Size(); i++) {
        Kps(i) = 0.01929 + 0.7045 * std::tanh(2.399 * Cn(i)) - 0.7238 * std::tanh(2.412 * Cn(i));
    }
    return Kps;
}