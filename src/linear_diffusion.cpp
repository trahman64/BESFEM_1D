#include "mfem.hpp"
#include "../includes/linear_diffusion.hpp"
#include <iostream>

LinearDiffusion::LinearDiffusion(mfem::Mesh *mesh_, mfem::FiniteElementSpace *fespace_,
                                mfem::Array<int> &nbc_bdr_, double C0_, double dt_)
    : mesh(mesh_), fespace(fespace_), nbc_bdr(nbc_bdr_), C0(C0_), dt(dt_),
      C(fespace_), rhs(fespace_->GetTrueVSize()), X(fespace_->GetTrueVSize()),
      C_prev(fespace_->GetTrueVSize()),
      rxn_lf(fespace_), vol_lf(fespace_)
{
    mfem::ConstantCoefficient one(1.0);
    rxn_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    rxn_lf.Assemble();

    C = C0;
    C.GetTrueDofs(C_prev);
}

LinearDiffusion::~LinearDiffusion() {
    delete M;
    delete K;
    delete TmatR;
    delete TmatL;
    delete prec;
    delete epsilon_coeff;
    delete region_weight;
    delete D_coeff;
    delete weight_eff_D;
}

void LinearDiffusion::SetWeightVector(double eps_sep, double eps_eld,
                                      double tau_sep, double tau_eld) {
    epsilon_vector.SetSize(mesh->attributes.Max());
    epsilon_vector(0) = eps_sep;
    epsilon_vector(1) = eps_eld;
    epsilon_coeff = new mfem::PWConstCoefficient(epsilon_vector);

    weight_vector.SetSize(mesh->attributes.Max());
    weight_vector(0) = (eps_sep * eps_sep) / (tau_sep * tau_sep);
    weight_vector(1) = (eps_eld * eps_eld) / (tau_eld * tau_eld);
    region_weight = new mfem::PWConstCoefficient(weight_vector);

    vol_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(*epsilon_coeff));
    vol_lf.Assemble();
}

void LinearDiffusion::SetCoefficient(mfem::GridFunction &D) {
    D_coeff = new mfem::GridFunctionCoefficient(&D);
    weight_eff_D = new mfem::ProductCoefficient(*region_weight, *D_coeff);
}

void LinearDiffusion::BuildOperator() {
    M = new mfem::BilinearForm(fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(*epsilon_coeff));
    M->Assemble();
    M->Finalize();

    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(*weight_eff_D));
    K->Assemble();
    K->Finalize();

    M_mat = &M->SpMat();
    K_mat = &K->SpMat();

    TmatR = Add(1.0, *M_mat, -0.5 * dt, *K_mat);
    TmatL = Add(1.0, *M_mat,  0.5 * dt, *K_mat);

    prec = new mfem::GSSmoother(*TmatL);
    solver.SetOperator(*TmatL);
    solver.SetPreconditioner(*prec);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(0.0);
    solver.SetMaxIter(500);
    solver.SetPrintLevel(0);
    solver.iterative_mode = true;
}

void LinearDiffusion::UpdateOperator() {
    // Only D's values changed (same mesh/sparsity) -- reassemble K, then
    // rebuild the Crank-Nicolson matrices that depend on it. M is untouched
    // since it depends only on epsilon_coeff, which is fixed.
    K->Update();
    K->Assemble();
    K->Finalize();
    K_mat = &K->SpMat();

    delete TmatR;
    delete TmatL;
    TmatR = Add(1.0, *M_mat, -0.5 * dt, *K_mat);
    TmatL = Add(1.0, *M_mat,  0.5 * dt, *K_mat);

    solver.SetOperator(*TmatL);
    // Optionally rebuild prec if convergence degrades over many updates:
    // delete prec;
    // prec = new mfem::GSSmoother(*TmatL);
    // solver.SetPreconditioner(*prec);
}

void LinearDiffusion::Stepping(mfem::GridFunction &source) {
    mfem::GridFunctionCoefficient reaction(&source);
    mfem::LinearForm R_current(fespace);
    R_current.AddDomainIntegrator(new mfem::DomainLFIntegrator(reaction));

    double f_in = rxn_lf(source);
    f_in *= -1.0;
    mfem::ConstantCoefficient nbcCoef(f_in);

    R_current.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(nbcCoef), nbc_bdr);   // uses the member now

    R_current.Assemble();

    rhs = R_current;
    rhs *= dt;

    TmatR->Mult(C_prev, X);
    X += rhs;

    solver.Mult(X, C_prev);
    C.SetFromTrueDofs(C_prev);
}

double LinearDiffusion::GetMeanConcentration() {
    double integral_epsC = vol_lf(C);
    double integral_eps  = vol_lf.Sum();
    return integral_epsC / integral_eps;
}

void LinearDiffusion::SaveMesh() {
    mesh->Save("linear_mesh.mesh");
}

void LinearDiffusion::SaveConc(const std::string &filename) {
    C.Save(filename.c_str());
}

mfem::GridFunction& LinearDiffusion::GetConcentration() {
    return C;
}