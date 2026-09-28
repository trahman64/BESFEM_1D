#include "mfem.hpp"
#include "../includes/atnv_calculator.hpp"

AtnVCalculator::AtnVCalculator(mfem::FiniteElementSpace *fespace_, 
	mfem::GridFunction &Dmp)
    : fespace(fespace_), Dmp_cf(&Dmp), K(fespace_)
{
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(Dmp_cf));
    K.Assemble();
    K.Finalize();
    AtnV.SetSize(fespace->GetTrueVSize());
}

mfem::Vector& AtnVCalculator::Compute(mfem::GridFunction &Cn) {
    Cn.GetTrueDofs(C_vec);
    K.SpMat().Mult(C_vec, AtnV);
    return AtnV;
}

void AtnVCalculator::UpdateDmp() {
    K.Update();
    K.Assemble();
    K.Finalize();
}