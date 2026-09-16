#include "mfem.hpp"
#include "../includes/radial_diffusion.hpp"
#include <iostream>

double r_square(const mfem::Vector &x){
    double r = x(0);
    return r * r;
}

double diff_coeff(const mfem::Vector &x){
    double d = 1.5e-10;
    double r = x(0);
    return d * r * r;
}

Radial_Diffusion::Radial_Diffusion(mfem::Mesh *mesh, mfem::FiniteElementSpace *fec) : mesh(mesh),fespace(fec),C(fec),MC(fec->GetTrueVSize()),KC(fec->GetTrueVSize()),rhs(fec->GetTrueVSize()),C_new(fec->GetTrueVSize()){

    C = 0.3;

    mfem::FunctionCoefficient r_coeff(r_square);
    mfem::FunctionCoefficient diffusion_coeff(diff_coeff);

    // M 
    M = new mfem::BilinearForm(fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(r_coeff));
    M->Assemble();
    M->Finalize();

    // K
    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(diffusion_coeff));
    K->Assemble();
    K->Finalize();


    //In matrix
    M_mat = &M->SpMat();
    K_mat = &K->SpMat();
}

    
void Radial_Diffusion::Stepping(double dt, double rxn){

     // R
    std::cout << "Radial 1" << std::endl;
    const double d = 1.5e-10;
    const double rho = 0.0312;
    const double particle_radius = 1.0;
    double surface_flux = d * particle_radius  * particle_radius* rxn / rho;
    std::cout << "Radial surface_flux" << std::endl;
    
    mfem::ConstantCoefficient flux_coeff(surface_flux);
    mfem::Array<int> out_bdr(mesh->bdr_attributes.Max());
    out_bdr = 0;
    out_bdr[1] = 1;
    std::cout << "Radial flux coeff" << std::endl;

    mfem::LinearForm R_current(fespace);
    R_current.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(flux_coeff),out_bdr);
    R_current.Assemble();
    std::cout << "Radial Assemble" << std::endl;

    M_mat->Mult(C,MC);
    K_mat->Mult(C,KC);

    rhs = R_current;
    rhs -= KC;
    rhs *= dt;
    rhs += MC;
    std::cout << "Radial rhs" <<std::endl;
    
    mfem::GSSmoother M_prec(*M_mat);
    mfem::CGSolver solver;

    solver.SetOperator(*M_mat);
    solver.SetPreconditioner(M_prec);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(0.0);
    solver.SetMaxIter(500);
    solver.SetPrintLevel(0);

    

    C_new = C;
    solver.Mult(rhs, C_new);
    C = C_new;
}

void Radial_Diffusion::Save(){
    mesh->Save("radial_mesh.mesh");
    C.Save("radial_concentration.gf");
}

double Radial_Diffusion::GetSurfaceConc(){
    return C[C.Size()-1];
}

    
    