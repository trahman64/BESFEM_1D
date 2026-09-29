#include "mfem.hpp"
#include "../includes/electrolyte_potential.hpp"
#include <iostream>

LiquidPotential::LiquidPotential(
    mfem::Mesh *mesh_,
    mfem::FiniteElementSpace *fespace_,
    LinearDiffusion *diffusion_,
    double a_,
    double F_,
    double R_,
    double T_,
    double epsilon_s_,
    double epsilon_e_,
    double tau_e_,
    double D_plus_,
    double D_minus_
)
    : mesh(mesh_),
      fespace(fespace_),
      diffusion(diffusion_),
      a(a_),
      F(F_),
      R(R_),
      T(T_),
      epsilon_s(epsilon_s_),
      epsilon_e(epsilon_e_),
      tau_e(tau_e_),
      D_plus(D_plus_),
      D_minus(D_minus_),
      phi_e(fespace_)
{
    phi_e = 0.0;
}


void LiquidPotential::Solve(double rxn)
{
    // Separator properties
    const double epsilon_e_sep = 1.0;
    const double tau_e_sep = 1.0;

    const double z_plus = 1.0;
    const double z_minus = -1.0;


    mfem::GridFunction &C =
        diffusion->GetConcentration();



    mfem::Array<int> ess_bdr(
        mesh->bdr_attributes.Max());

    ess_bdr = 0;
    ess_bdr[0] = 1;

    mfem::Array<int> ess_tdof_list;

    fespace->GetEssentialTrueDofs(
        ess_bdr,
        ess_tdof_list);



    double m_plus =
        D_plus / (R * T);

    double m_minus =
        D_minus / (R * T);


    double mobility_factor =
        ((z_plus * m_plus) -
         (z_minus * m_minus)) * F;


    mfem::GridFunctionCoefficient C_coeff(&C);


    mfem::Vector conductivity_values(
        mesh->attributes.Max());

    conductivity_values = 0.0;

    // Region 1: separator
    conductivity_values(0) =
        mobility_factor *
        epsilon_e_sep /
        (tau_e_sep * tau_e_sep);

    // Region 2: electrode
    conductivity_values(1) =
        mobility_factor *
        epsilon_e /
        (tau_e * tau_e);


    mfem::PWConstCoefficient conductivity_factor(
        conductivity_values);


    mfem::ProductCoefficient kappa_e_bar(
        conductivity_factor,
        C_coeff);


    double D_pm =
        D_plus - D_minus;


    mfem::Vector Dpm_values(
        mesh->attributes.Max());

    Dpm_values = 0.0;

    // Region 1: separator
    Dpm_values(0) =
        D_pm *
        epsilon_e_sep /
        (tau_e_sep * tau_e_sep);

    // Region 2: electrode
    Dpm_values(1) =
        D_pm *
        epsilon_e /
        (tau_e * tau_e);


    mfem::PWConstCoefficient Dpm_bar(
        Dpm_values);



    mfem::BilinearForm form(fespace);

    form.AddDomainIntegrator(
        new mfem::DiffusionIntegrator(kappa_e_bar));

    form.AddDomainIntegrator(
        new mfem::DiffusionIntegrator(Dpm_bar));

    form.Assemble();



    mfem::Vector source_values(
        mesh->attributes.Max());

    source_values = 0.0;

    // Separator
    source_values(0) = 0.0;

    // Electrode
    source_values(1) =
        a * rxn * epsilon_s;


    mfem::PWConstCoefficient source(
        source_values);


    mfem::LinearForm b(fespace);

    b.AddDomainIntegrator(
        new mfem::DomainLFIntegrator(source));

    b.Assemble();


    // Weak form sign
    b *= -1.0;


 

    mfem::ConstantCoefficient dbc_coeff(0.0);

    phi_e.ProjectBdrCoefficient(
        dbc_coeff,
        ess_bdr);


    // Linear system

    mfem::SparseMatrix A;
    mfem::Vector X;
    mfem::Vector B;

    form.FormLinearSystem(
        ess_tdof_list,
        phi_e,
        b,
        A,
        X,
        B);

    // Solve

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


    form.RecoverFEMSolution(
        X,
        b,
        phi_e);
}


void LiquidPotential::Save(
    const std::string &filename)
{
    phi_e.Save(filename.c_str());
}


const mfem::GridFunction&
LiquidPotential::GetPotential() const
{
    return phi_e;
}