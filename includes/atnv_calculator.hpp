#ifndef ATNV_CALCULATOR_HPP
#define ATNV_CALCULATOR_HPP

#include "mfem.hpp"

class AtnVCalculator {
public:
    AtnVCalculator(mfem::FiniteElementSpace *fespace, mfem::GridFunction &Dmp);

    mfem::Vector& Compute(mfem::GridFunction &Cn);
    void UpdateDmp();

private:
    mfem::FiniteElementSpace *fespace;
    mfem::GridFunctionCoefficient Dmp_cf;
    mfem::BilinearForm K;
    mfem::Vector C_vec, AtnV;
};

#endif