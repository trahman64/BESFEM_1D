#include "mfem.hpp"
#include "../includes/diffusion.hpp"
#include "../includes/diffusion.hpp"
#include <iostream>


Diffusion::Diffusion(mfem::Mesh *mesh,mfem::FiniteElementSpace *fec):mesh(mesh),fespace(fec),C(fec),MC(fec->GetTrueVSize()),KC(fec->GetTrueVSize()),rhs(fec->GetTrueVSize()),C_prev(fec->GetTrueVSize()),C_new(fec->GetTrueVSize()){

    mfem::Vector epsilon_vector(mesh->attributes.Max());
    mfem::Vector diffusion_vector(mesh->attributes.Max());
    mfem::Vector reaction_vector(mesh->attributes.Max());

    //Region 1
    double De = 0.25e-5;
    const double tau_sep = 1.0;
    const double tau_electrode = 1.521;
    epsilon_vector(0) = 1.0;
    diffusion_vector(0) =  De * epsilon_vector(0) / (tau_sep * tau_sep);;

    //Region 2
    epsilon_vector(1) = 0.301;
    diffusion_vector(1) = epsilon_vector(1) * De /(tau_electrode * tau_electrode);

    mfem::PWConstCoefficient epsilon(epsilon_vector);
    mfem::PWConstCoefficient diffusion(diffusion_vector);

    //M
    M = new mfem::BilinearForm(fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(epsilon));
    M->Assemble();
    M->Finalize();

    //K
    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(diffusion));
    K->Assemble();
    K->Finalize();

    

    
    // mfem::ConstantCoefficient rxn_coeff(a * t_minus * rxn * 0.5);
    // mfem::Array<int> bdr(mesh->bdr_attributes.Max());
    // bdr[0]=1;
    // R->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(rxn_coeff),bdr);

    
    //In matrix
    C = 0.001;
    M_mat = &M->SpMat();
    K_mat = &K->SpMat();

    C.GetTrueDofs(C_prev);
}

void Diffusion::Stepping(double dt,double rxn){
    const double a_rxn = 2.409e3;
    const double t_minus = 0.7619;
    const double eps_s = 0.699;

    mfem::Vector reaction_vector(mesh->attributes.Max());
    reaction_vector = 0.0;

    //Separator
    reaction_vector(0) = 0.0;
    
    //Electrode
    reaction_vector(1) = a_rxn * rxn * t_minus * eps_s;

    mfem::PWConstCoefficient reaction(reaction_vector);
    mfem::LinearForm R_current(fespace);
    R_current.AddDomainIntegrator(new mfem::DomainLFIntegrator(reaction));
    R_current.Assemble();

    
    double electrode_length = 0.5;
    double f_in = a_rxn * rxn * t_minus * eps_s * electrode_length;

    mfem::Array<int> ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[0] = 1;
    ess_bdr[1] = 1;

    mfem::Array<int> ess_tdof_list;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    mfem::Vector boundary_values(mesh->bdr_attributes.Max());
    boundary_values(0) = f_in;
    boundary_values(1) = 0.0;

    mfem::PWConstCoefficient boundary_coeff(boundary_values);
    C.ProjectBdrCoefficient(boundary_coeff,ess_bdr);
    C.GetTrueDofs(C_prev);

    M_mat->Mult(C_prev, MC);
    K_mat->Mult(C_prev, KC);
    rhs = R_current;
    rhs -= KC;
    rhs *=dt;
    rhs += MC;

    mfem::SparseMatrix A;
    mfem::Vector X;
    mfem::Vector B;

    M->FormLinearSystem(ess_tdof_list,C,rhs,A,X,B);
    
    mfem::GSSmoother M_prec(A);
    mfem::CGSolver solver;

    solver.SetOperator(A);
    solver.SetPreconditioner(M_prec);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(0.0);
    solver.SetMaxIter(500);
    solver.SetPrintLevel(0);


    
    
    solver.Mult(B,X);
        
    C.GetTrueDofs(C_prev);
}
void Diffusion::Save(){
    mesh->Save("diffusion_mesh.mesh");
    C.Save("diffusion_concentration.gf");
}

mfem::GridFunction& Diffusion::GetConcentration(){
    return C;
}
    
    
