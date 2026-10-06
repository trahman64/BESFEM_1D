#include "../includes/spherical_diffusion.hpp"
#include "../includes/diffCoeff_utils.hpp"   
#include <iostream>

SphericalDiffusion::SphericalDiffusion(int x_idx_, 
									   double radius,
                                       int n_elements,
                                       int fe_order, 
                                       double C0,
                                       double dt_,
                                       TabulatedData &DLi_table_)
    : R(radius), order(fe_order), x_idx(x_idx_), dt(dt_), DLi_table(DLi_table_)
{
    BuildMesh(n_elements);

    fec = new mfem::H1_FECollection(order, fe_order);   // 1D H1 elements
    fespace = new mfem::FiniteElementSpace(mesh, fec);

    C.SetSpace(fespace);
    C = C0;
    C_prev.SetSize(fespace->GetTrueVSize());
    C.GetTrueDofs(C_prev);
    D_li.SetSpace(fespace); 
    vol_lf.Update(fespace);

    BuildOperators();

    // MakeCartesian1D convention: left boundary (r=0) = attribute 1,
    // right boundary (r=R) = attribute 2.
    outer_bdr_marker.SetSize(mesh->bdr_attributes.Max());
    outer_bdr_marker = 0;
    outer_bdr_marker[1] = 1;   // attribute 2 -> index 1

	partiVolume = R * R * R / 3.0;   // ∫_0^R r^2 dr
	partiSurfArea = R * R; // "raw" r^2 at the surface, consistent with GetParticleVolume()'s un-4π convention
}

SphericalDiffusion::~SphericalDiffusion() {
    delete M;
    delete K;
    delete Tmat;
    delete prec;
    delete r2_coeff;
    delete D_li_coeff;
    delete D_r2;
    delete fespace;
    delete fec;
    delete mesh;
}

void SphericalDiffusion::BuildMesh(int n_elements) {
    mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian1D(n_elements, R));
}

void SphericalDiffusion::BuildOperators() {
    r.SetSpace(fespace);
    mesh->GetNodes(r);
    r *= r;
    r2_coeff = new mfem::GridFunctionCoefficient(&r);

    D_li = Compute_DLi(C, DLi_table);
    D_li_coeff = new mfem::GridFunctionCoefficient(&D_li);
    D_r2 = new mfem::ProductCoefficient(*D_li_coeff, *r2_coeff);

    M = new mfem::BilinearForm(fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(*r2_coeff));
    M->Assemble();
    M->Finalize();

    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(*D_r2));
    K->Assemble();
    K->Finalize();

    M_mat = &M->SpMat();
    K_mat = &K->SpMat();
    
    // --- merged in from Initialize() ---
    prec = new mfem::GSSmoother(*M_mat);
    solver.SetOperator(*M_mat);
    solver.SetPreconditioner(*prec);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(0.0);
    solver.SetMaxIter(500);
    solver.SetPrintLevel(0);
    
    Tmat = Add(1.0, *M_mat, -dt, *K_mat);
    
	// in BuildOperators(), after r2_coeff is built:
	vol_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(*r2_coeff));
	vol_lf.Assemble();        
}


void SphericalDiffusion::UpdateOperator() {
    D_li = Compute_DLi(C, DLi_table);   // D_li_coeff already points at this member -- sees new values automatically

    K->Update();
    K->Assemble();
    K->Finalize();
    K_mat = &K->SpMat();

    delete Tmat;
    Tmat = Add(1.0, *M_mat, -dt, *K_mat);
}


void SphericalDiffusion::Stepping(double surface_flux) {
    // Weak-form boundary term at r=R: [D r^2 dC/dr v]_{r=R} = R^2 * surface_flux * v(R)
    mfem::ConstantCoefficient flux_coeff(surface_flux * R * R);

    mfem::LinearForm R_bc(fespace);
    R_bc.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(flux_coeff), outer_bdr_marker);
    R_bc.Assemble();

    mfem::Vector rhs(R_bc);
    rhs *= dt;
        
    mfem::Vector X(fespace->GetTrueVSize());
    Tmat->Mult(C_prev, X);   // X = (M - dt*K) * C^n   -- correct now, Tmat has the minus sign
    X += rhs;                // X = (M - dt*K)*C^n + dt*f

    solver.Mult(X, C_prev);  // solves M * C^{n+1} = X
    
    C.SetFromTrueDofs(C_prev);
}

double SphericalDiffusion::GetMeanConcentration() {
    double integral_Cr2 = vol_lf(C);
    return integral_Cr2 / partiVolume;
}

double SphericalDiffusion::GetParticleVolume() {
    return partiVolume;
}

double SphericalDiffusion::GetParticleSurfaceArea() {
    return partiSurfArea;   
}


double SphericalDiffusion::GetParticleTotalLi() {
    double integral_Cr2 = vol_lf(C);
    return integral_Cr2;
}   

double SphericalDiffusion::GetConcentrationAt(double x) {
    mfem::Array<int> elem_ids;
    mfem::Array<mfem::IntegrationPoint> ips;
    mfem::DenseMatrix point_mat(1, 1);
    point_mat(0, 0) = x;

    mesh->FindPoints(point_mat, elem_ids, ips);

    if (elem_ids[0] < 0) {
        std::cerr << "GetConcentrationAt: x=" << x << " is outside the mesh [0, "
                  << R << "]" << std::endl;
        return 0.0;
    }

    return C.GetValue(elem_ids[0], ips[0]);
}

int SphericalDiffusion::GetParticleID() const {
    return x_idx;
}

mfem::GridFunction& SphericalDiffusion::GetConcentration() {
    return C;
}

void SphericalDiffusion::SaveConc(const std::string &prefix) {
    C.Save((prefix + std::to_string(x_idx) + "_concentration.gf").c_str());
}

void SphericalDiffusion::SaveMesh(const std::string &prefix) {
    mesh->Save((prefix + std::to_string(x_idx) + "_mesh.mesh").c_str());
}