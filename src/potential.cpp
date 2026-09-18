#include "mfem.hpp"
#include <iostream>
#include <fstream>

#include "../includes/potential.hpp"
#include "../includes/diffusion.hpp"

using namespace std;
using namespace mfem;

Potential::Potential(mfem::Mesh *mesh, mfem::FiniteElementSpace *fec, Diffusion *diffusion) : mesh(mesh), fespace(fec), diffusion(diffusion), phi_s(fec), phi_e(fec){
    phi_s = 0.0;
    phi_e = 0.0;
}

void Potential::SolveSolid(double rxn){
    const double F = 96485.332;
    const double a_rxn = 2.409e3;
    const double kappa_s = 0.075;
    const double eps_s = 0.699;
    const double tau_s = 1.324;
    const double eps_s_sep = 1e-4;
    const double tau_s_sep = 1e4;
    double kappa_s_eff = kappa_s * eps_s / (tau_s * tau_s) * eps_s;
    double kappa_s_sep = kappa_s * eps_s_sep / (tau_s_sep * tau_s_sep);
    
    double BvP = 0.0;
    
    Array<int> ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[1] = 1;
    Array<int> ess_tdof_list;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    Vector kappa_values(mesh->attributes.Max());
    kappa_values = 0.0;
    kappa_values(0) = kappa_s_sep;    
    kappa_values(1) = kappa_s_eff;   

    PWConstCoefficient kappa(kappa_values);

    Vector source_values(mesh->attributes.Max());

    source_values = 0.0;
    source_values(0) = 0.0;  
//     source_values(1) = a_rxn * rxn * F * eps_s; 
    source_values(1) = a_rxn * rxn * F;  

    PWConstCoefficient source(source_values);

    LinearForm b(fespace);
    b.AddDomainIntegrator(new DomainLFIntegrator(source));
    b.Assemble();


    BilinearForm a(fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(kappa));
    a.Assemble();

    SparseMatrix A;
    Vector X;
    Vector B;
    
    mfem::ConstantCoefficient dbc_coeff(BvP);
	phi_s.ProjectBdrCoefficient(dbc_coeff, ess_bdr);

    a.FormLinearSystem(ess_tdof_list,phi_s,b,A,X,B);
    GSSmoother M(A);
    PCG(A,M, B,X,1,200,1e-12,0.0);
    a.RecoverFEMSolution(X,b,phi_s);

}

void Potential::SolveElectrolyte(double rxn){
    const double R = 8.314;
    const double T = 300.0;
    const double F = 96485.332;

    const double a_rxn = 2.409e3;
    const double eps_s = 0.699;

    const double eps_e_sep = 1.0;
    const double eps_e_electrode = 0.301;

    const double tau_e_sep = 1.0;
    const double tau_e_electrode = 1.521;

    const double D_plus = 0.125e-5;
    const double D_minus = 0.4e-5;

    const double z_plus = 1.0;
    const double z_minus = -1.0;


    mfem::GridFunction &C = diffusion->GetConcentration();
    Array<int> ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[0] = 1;
    Array<int> ess_tdof_list;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    double m_plus = D_plus / (R * T);
    double m_minus = D_minus / (R * T);
    double mobility_factor = ((z_plus * m_plus) - (z_minus * m_minus)) * F;

    mfem::GridFunctionCoefficient C_coeff(&C);
    mfem::Vector conductivity_values(mesh->attributes.Max());
    conductivity_values = 0.0;

    conductivity_values(0) = mobility_factor * eps_e_sep / (tau_e_sep * tau_e_sep) ;
    conductivity_values(1) = mobility_factor * eps_e_electrode / (tau_e_electrode * tau_e_electrode);

    mfem::PWConstCoefficient conductivity_factor(conductivity_values);
    mfem::ProductCoefficient kappa_e_bar(conductivity_factor, C_coeff);

    mfem::BilinearForm a(fespace);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(kappa_e_bar));

    double D_pm = D_plus - D_minus;
    mfem::Vector Dpm_values(mesh->attributes.Max());
    Dpm_values = 0.0;
    Dpm_values(0) = D_pm * (eps_e_sep) / (tau_e_sep * tau_e_sep);
    Dpm_values(1) = D_pm * (eps_e_electrode) / (tau_e_electrode * tau_e_electrode);
    mfem::PWConstCoefficient Dpm_bar(Dpm_values);
    
    a.AddDomainIntegrator(
    new mfem::DiffusionIntegrator(Dpm_bar));

    a.Assemble();
    
    mfem::Vector source_values(mesh->attributes.Max());
    source_values = 0.0;
    source_values(0) = 0.0;
    source_values(1) = a_rxn * rxn * eps_s;

    mfem::PWConstCoefficient source(source_values);

    mfem::LinearForm b(fespace);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(source));
    b.Assemble();

    b*= -1.0;

    mfem::SparseMatrix A;
    mfem::Vector X;
    mfem::Vector B;
    
    a.FormLinearSystem(ess_tdof_list,phi_e,b,A,X,B);
    GSSmoother M(A);
    PCG(A,M, B,X,1,200,1e-12,0.0);
    a.RecoverFEMSolution(X,b,phi_e);
    
}
void Potential::Solve(double rxn){
    std::cout << "  Entering SolveSolid" << std::endl;
    SolveSolid(rxn);
    std::cout << "  Finishing SolveSolid" << std::endl;
    
    std::cout << "  Entering SolveLiquid" << std::endl;
    SolveElectrolyte(rxn);
    std::cout << "  Finishing SolveLiquid" << std::endl;
    
}
   
void Potential::Save()
{
    mesh->Save("mesh.mesh");
    phi_s.Save("solid_potential.gf");
    phi_e.Save("electrolyte_potential.gf");
}

const mfem::GridFunction& Potential::GetSolidPotential() const
{
    return phi_s;
}

const mfem::GridFunction& Potential::GetLiquidPotential() const
{
    return phi_e;
}

