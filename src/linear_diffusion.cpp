#include "mfem.hpp"
#include "../includes/linear_diffusion.hpp"
#include <iostream>

LinearDiffusion::LinearDiffusion(mfem::Mesh *mesh_, mfem::FiniteElementSpace *fespace_,
                      double epsilon_sep_, double epsilon_eld_, double De_,
                      double tau_sep_, double tau_eld_, double t_minus_,
                      double C0_, double dt_)
    : mesh(mesh_), fespace(fespace_),
      epsilon_sep(epsilon_sep_), epsilon_eld(epsilon_eld_), De(De_),
      tau_sep(tau_sep_), tau_eld(tau_eld_), t_minus(t_minus_), C0(C0_), dt(dt_),
      C(fespace_), rhs(fespace_->GetTrueVSize()), X(fespace_->GetTrueVSize()),
      C_prev(fespace_->GetTrueVSize()),
      rxn_lf(fespace_), vol_lf(fespace_)
{
    epsilon_vector.SetSize(mesh->attributes.Max());          // MEMBER, not local
    mfem::Vector diffusion_vector(mesh->attributes.Max());   // local is fine -- used once

    epsilon_vector(0) = epsilon_sep;
    diffusion_vector(0) = De * epsilon_vector(0) / (tau_sep * tau_sep);

    epsilon_vector(1) = epsilon_eld;
    diffusion_vector(1) = De * epsilon_vector(1) / (tau_eld * tau_eld);

    epsilon_coeff = new mfem::PWConstCoefficient(epsilon_vector);   // now references the MEMBER
    mfem::PWConstCoefficient diffusion_coeff(diffusion_vector);

    M = new mfem::BilinearForm(fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(*epsilon_coeff));
    M->Assemble();
    M->Finalize();

    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(diffusion_coeff));
    K->Assemble();
    K->Finalize();

    M_mat = &M->SpMat();
    K_mat = &K->SpMat();

    mfem::ConstantCoefficient one(1.0);
    rxn_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    rxn_lf.Assemble();

    vol_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(*epsilon_coeff));
    vol_lf.Assemble();

    C = C0;
    C.GetTrueDofs(C_prev);

    TmatR = Add(1.0, *M_mat, -0.5 * dt, *K_mat);
    TmatL = Add(1.0, *M_mat,  0.5 * dt, *K_mat);

    prec = new mfem::GSSmoother(*TmatL);
    solver.SetOperator(*TmatL);
    solver.SetPreconditioner(*prec);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(0.0);
    solver.SetMaxIter(500);
    solver.SetPrintLevel(0);
}

LinearDiffusion::~LinearDiffusion() {
    delete M;
    delete K;
    delete TmatR;
    delete TmatL;
    delete prec;
    delete epsilon_coeff;   // added
}

void LinearDiffusion::Stepping(mfem::GridFunction source) {
    source *= t_minus;

    mfem::GridFunctionCoefficient reaction(&source);   // FIX: pointer, not object
    mfem::LinearForm R_current(fespace);
    R_current.AddDomainIntegrator(new mfem::DomainLFIntegrator(reaction));

    double f_in = rxn_lf(source);   // now valid -- rxn_lf is a member
    mfem::ConstantCoefficient nbcCoef(f_in);

    mfem::Array<int> nbc_w_bdr(mesh->bdr_attributes.Max());
    nbc_w_bdr = 0;
    nbc_w_bdr[0] = 1;

    R_current.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(nbcCoef), nbc_w_bdr);
    R_current.Assemble();

    rhs = R_current;
    rhs *= dt;

    TmatR->Mult(C_prev, X);   // X is now properly sized
    X += rhs;

    solver.Mult(X, C_prev);
    C.SetFromTrueDofs(C_prev);
}


double LinearDiffusion::GetMeanConcentration() {

    double integral_epsC = vol_lf(C);        // ∫ epsilon * C dx
    double integral_eps  = vol_lf.Sum();      // ∫ epsilon * 1 dx  (since Sum() adds up the assembled entries,
                                                //  each of which is ∫ epsilon * phi_i dx; their sum over ALL i
                                                //  telescopes to ∫ epsilon dx via partition-of-unity of phi_i's)
    return integral_epsC / integral_eps;
}


void LinearDiffusion::SaveMesh() {
    mesh->Save("diffusion_mesh.mesh");
}

void LinearDiffusion::SaveConc() {
    C.Save("diffusion_concentration.gf");
}

mfem::GridFunction& LinearDiffusion::GetConcentration() {
    return C;
}