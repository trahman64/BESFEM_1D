#include "mfem.hpp"
#include "../includes/stat_potential.hpp"
#include <iostream>

StatPotential::StatPotential(mfem::FiniteElementSpace *fespace_,
                              mfem::Coefficient &effective_kappa,
                              mfem::Array<int> &ess_bdr_)
    : fespace(fespace_), ess_bdr(ess_bdr_),
      phi(fespace_), X(fespace_->GetTrueVSize()), B(fespace_->GetTrueVSize())
{
    phi = 0.0;

    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(effective_kappa));
    K->Assemble();
    K->Finalize();

    K_mat = &K->SpMat();

    prec = new mfem::GSSmoother(*K_mat);
    solver.SetOperator(*K_mat);
    solver.SetPreconditioner(*prec);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(0.0);
    solver.SetMaxIter(500);
    solver.SetPrintLevel(0);
}

StatPotential::~StatPotential() {
    delete K;
    delete prec;
}

void StatPotential::Solve(mfem::GridFunction &source, mfem::Vector &Additional, double Bv) {
    mfem::LinearForm b(fespace);

    mfem::GridFunctionCoefficient source_coeff(&source);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(source_coeff));
    b.Assemble();

    mfem::ConstantCoefficient dbc_coeff(Bv);
    phi.ProjectBdrCoefficient(dbc_coeff, ess_bdr);

    K->FormLinearSystem(ess_tdof_list, phi, b, A, X, B);

    B += Additional;

    solver.SetOperator(A);
    solver.Mult(B, X);

    K->RecoverFEMSolution(X, b, phi);
}

mfem::GridFunction& StatPotential::GetPotential() {
    return phi;
}

void StatPotential::SavePote() {
    phi.Save("phi_solid.gf");
    std::cout << "Saved phi_solid.gf" << std::endl;
}