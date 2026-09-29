#include "mfem.hpp"
#include "../includes/spherical_diffusion.hpp"
#include <iostream>



SphericalDiffusion::SphericalDiffusion(double radius, int n_elements, double diffusion_coef, int x_idx_, double dt_, int fe_order, double C0) : R(radius), D(diffusion_coef), order(fe_order), x_idx(x_idx_), dt(dt_){

    BuildMesh(n_elements);
    fec = new mfem::H1_FECollection(order, 1);
    fespace = new mfem::FiniteElementSpace(mesh, fec);

    C.SetSpace(fespace);
    C= C0;
    C_prev.SetSize(fespace->GetTrueVSize());
    C.GetTrueDofs(C_prev);

    BuildOperators();

    outer_bdr_marker.SetSize(mesh->bdr_attributes.Max());
    outer_bdr_marker = 0;
    outer_bdr_marker[1] = 1;

    Initialize();

    // mfem::FunctionCoefficient r_coeff(r_square);
    // mfem::FunctionCoefficient diffusion_coeff(diff_coeff);

    
}

SphericalDiffusion::~SphericalDiffusion(){
    delete M;
    delete K;
    delete Tmat;
    delete prec;
    delete fespace;
    delete fec;
    delete mesh;
}


void SphericalDiffusion::BuildMesh(int n_elements){
    mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian1D(n_elements,R));
}

void SphericalDiffusion::BuildOperators(){
    mfem::FunctionCoefficient r2_coeff([](const mfem::Vector &x){
        return x(0) * x(0);});
    
    mfem::ConstantCoefficient D_coeff(D);
    mfem::ProductCoefficient D_r2(D_coeff, r2_coeff);

    // M 
    M = new mfem::BilinearForm(fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(r2_coeff));
    M->Assemble();
    M->Finalize();

    // K
    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(D_r2));
    K->Assemble();
    K->Finalize();

     //In matrix
    M_mat = &M->SpMat();
    K_mat = &K->SpMat();
}

void SphericalDiffusion::Initialize(){
    prec = new mfem::GSSmoother(*M_mat);
    solver.SetOperator(*M_mat);
    solver.SetPreconditioner(*prec);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(0.0);
    solver.SetMaxIter(500);
    solver.SetPrintLevel(0);

    Tmat = Add(1.0, *M_mat, -dt, *K_mat);
}


void SphericalDiffusion::Stepping(double surface_flux){

     // R
    
    std::cout << "Radial surface_flux" << std::endl;
    
    mfem::ConstantCoefficient flux_coeff(surface_flux * R * R);
    
    mfem::LinearForm R_bc(fespace);
    R_bc.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(flux_coeff),outer_bdr_marker);
    R_bc.Assemble();
    std::cout << "Radial Assemble" << std::endl;

    mfem::Vector rhs(R_bc);
    rhs *= dt;

    mfem::Vector X(fespace->GetTrueVSize());
    Tmat->Mult(C_prev,X);
    X += rhs;
    
    solver.Mult(X, C_prev);
    C.SetFromTrueDofs(C_prev);
}

double SphericalDiffusion::GetMeanConcentration(){
    mfem::FunctionCoefficient r2_coeff([](const mfem::Vector &x){
        return x(0) * x(0);
    });

    mfem::LinearForm vol_lf(fespace);
    vol_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(r2_coeff));
    vol_lf.Assemble();

    double integral_Cr2 = vol_lf(C);
    double volume = R * R * R / 2.0;

    return integral_Cr2 / volume;

    
}

void SphericalDiffusion::SaveConc(const std::string &prefix) {
    C.Save((prefix + "_concentration.gf").c_str());
}

void SphericalDiffusion::SaveMesh(const std::string &prefix) {
    mesh->Save((prefix + "_mesh.mesh").c_str());
}

double SphericalDiffusion::GetConcentrationAt(double x){
    mfem::Array<int> elem_ids;
    mfem::Array<mfem::IntegrationPoint> ips;
    mfem::DenseMatrix point_mat(1,1);

    mesh->FindPoints(point_mat, elem_ids, ips);

    if (elem_ids[0] < 0) {
        std::cerr << "GetConcentrationAt: x =" << x << " is outside the mesh [0, " << R << "]" << std::endl;
        return 0.0;
    }
    return C.GetValue(elem_ids[0], ips[0]);
}

    
    