#include "mfem.hpp"
#include "../includes/solid_potential.hpp"
#include <iostream>

SolidPotential::SolidPotential(
    mfem::Mesh *mesh_,
    mfem::FiniteElementSpace *fespace_)
    : mesh(mesh_),
      fespace(fespace_),
      phi_s(fespace_)
{
    phi_s = 0.0;
}


void SolidPotential::Solve(double rxn)
{
    const double F = 96485.332;

    const double a_rxn = 2.409e3;

    const double kappa_s = 0.075;

    // Electrode
    const double eps_s = 0.699;
    const double tau_s = 1.324;

    // Separator
    const double eps_s_sep = 1e-4;
    const double tau_s_sep = 1e4;


    // ==========================================
    // Effective solid conductivity
    // ==========================================

    double kappa_s_eff =
        kappa_s * eps_s / (tau_s * tau_s);

    double kappa_s_sep =
        kappa_s * eps_s_sep /
        (tau_s_sep * tau_s_sep);


    // ==========================================
    // Boundary condition
    // ==========================================

    mfem::Array<int> ess_bdr(
        mesh->bdr_attributes.Max());

    ess_bdr = 0;

    // Boundary attribute 2
    ess_bdr[1] = 1;

    mfem::Array<int> ess_tdof_list;

    fespace->GetEssentialTrueDofs(
        ess_bdr,
        ess_tdof_list);


    // ==========================================
    // Conductivity coefficient
    // ==========================================

    mfem::Vector kappa_values(
        mesh->attributes.Max());

    kappa_values = 0.0;

    // Region 1: separator
    kappa_values(0) = kappa_s_sep;

    // Region 2: electrode
    kappa_values(1) = kappa_s_eff;

    mfem::PWConstCoefficient kappa(
        kappa_values);


    // ==========================================
    // Source
    // ==========================================

    mfem::Vector source_values(
        mesh->attributes.Max());

    source_values = 0.0;

    // Separator
    source_values(0) = 0.0;

    // Electrode
    source_values(1) =
        a_rxn * rxn * F * eps_s;

    mfem::PWConstCoefficient source(
        source_values);


    mfem::LinearForm b(fespace);

    b.AddDomainIntegrator(
        new mfem::DomainLFIntegrator(source));

    b.Assemble();


    // ==========================================
    // Bilinear form
    // ==========================================

    mfem::BilinearForm a(fespace);

    a.AddDomainIntegrator(
        new mfem::DiffusionIntegrator(kappa));

    a.Assemble();


    // ==========================================
    // Dirichlet BC
    // ==========================================

    mfem::ConstantCoefficient dbc_coeff(0.0);

    phi_s.ProjectBdrCoefficient(
        dbc_coeff,
        ess_bdr);


    // ==========================================
    // Linear system
    // ==========================================

    mfem::SparseMatrix A;
    mfem::Vector X;
    mfem::Vector B;

    a.FormLinearSystem(
        ess_tdof_list,
        phi_s,
        b,
        A,
        X,
        B);


    // ==========================================
    // Solve
    // ==========================================

    mfem::GSSmoother M(A);

    mfem::PCG(
        A,
        M,
        B,
        X,
        1,
        200,
        1e-12,
        0.0);


    a.RecoverFEMSolution(
        X,
        b,
        phi_s);
}


void SolidPotential::Save(
    const std::string &filename)
{
    phi_s.Save(filename.c_str());
}


const mfem::GridFunction&
SolidPotential::GetPotential() const
{
    return phi_s;
}