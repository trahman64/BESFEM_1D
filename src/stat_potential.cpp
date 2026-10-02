#include "mfem.hpp"
#include "../includes/stat_potential.hpp"
#include <iostream>

StatPotential::StatPotential(mfem::FiniteElementSpace *fespace_,
                              mfem::Array<int> &ess_bdr_, double BcV_)
    : fespace(fespace_), ess_bdr(ess_bdr_), BcV(BcV_),
      phi(fespace_), X(fespace_->GetTrueVSize()), B(fespace_->GetTrueVSize())
{
    phi = BcV;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
}

StatPotential::~StatPotential() {
    delete K;
    delete prec;
    delete region_weight;
    delete kappa_coeff;
    delete weight_eff_kappa;
}

void StatPotential::SetWeightVector(double eps_sep, double eps_eld,
                                    double tau_sep, double tau_eld) {
    mfem::Mesh *mesh = fespace->GetMesh();

    weight_vector.SetSize(mesh->attributes.Max());
    weight_vector(0) = (eps_sep * eps_sep) / (tau_sep * tau_sep);   // separator
    weight_vector(1) = (eps_eld * eps_eld) / (tau_eld * tau_eld);   // electrode

    region_weight = new mfem::PWConstCoefficient(weight_vector);
}

void StatPotential::SetCoefficient(mfem::GridFunction &Kappa) {
    kappa_coeff = new mfem::GridFunctionCoefficient(&Kappa);
    weight_eff_kappa = new mfem::ProductCoefficient(*region_weight, *kappa_coeff);
}

void StatPotential::BuildOperator() {
    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(*weight_eff_kappa));
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
    solver.iterative_mode = true;
}

void StatPotential::UpdateOperator() {
    K->Update();      // resets K's matrix data, keeps its integrator list intact
    K->Assemble();     // re-evaluates weight_eff_kappa at each quadrature point, rebuilds entries
    K->Finalize();

    K_mat = &K->SpMat();       // re-point, cheap insurance (see below)
    solver.SetOperator(*K_mat); // cheap insurance, always safe
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

void StatPotential::SavePote(const std::string &filename) {
    phi.Save(filename.c_str());
    std::cout << "Saved " << filename << std::endl;
}