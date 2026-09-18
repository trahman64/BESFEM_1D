#include "../includes/spherical_diffusion.hpp"
#include <iostream>

SphericalDiffusion::SphericalDiffusion(double radius,
                                        int n_elements,
                                        double diffusion_coef,
                                        int x_idx_,
                                        double dt_,
                                        int fe_order,
                                        double C0)
    : R(radius), D(diffusion_coef), order(fe_order), x_idx(x_idx_), dt(dt_)
{
    BuildMesh(n_elements);

    fec = new mfem::H1_FECollection(order, 1);   // 1D H1 elements
    fespace = new mfem::FiniteElementSpace(mesh, fec);

    C.SetSpace(fespace);
    C = C0;
    C_prev.SetSize(fespace->GetTrueVSize());
    C.GetTrueDofs(C_prev);

    BuildOperators();

    // MakeCartesian1D convention: left boundary (r=0) = attribute 1,
    // right boundary (r=R) = attribute 2.
    outer_bdr_marker.SetSize(mesh->bdr_attributes.Max());
    outer_bdr_marker = 0;
    outer_bdr_marker[1] = 1;   // attribute 2 -> index 1

    Initialize();   // build Tmat + solver now that M_mat, K_mat, dt are all set
}

SphericalDiffusion::~SphericalDiffusion() {
    delete M;
    delete K;
    delete Tmat;
    delete prec;
    // M_mat / K_mat point into M / K's own storage -- not deleted separately
    delete fespace;
    delete fec;
    delete mesh;
}

void SphericalDiffusion::BuildMesh(int n_elements) {
    mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian1D(n_elements, R));
}

void SphericalDiffusion::BuildOperators() {
    mfem::FunctionCoefficient r2_coeff([](const mfem::Vector &x) {
        return x(0) * x(0);
    });

    mfem::ConstantCoefficient D_coeff(D);
    mfem::ProductCoefficient D_r2(D_coeff, r2_coeff);   // D * r^2  (stiffness weight)

    M = new mfem::BilinearForm(fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(r2_coeff));
    M->Assemble();
    M->Finalize();

    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(D_r2));
    K->Assemble();
    K->Finalize();

    M_mat = &M->SpMat();
    K_mat = &K->SpMat();
    
    std::cout << D << std::endl;    
    
}

void SphericalDiffusion::Initialize() {
    // Forward Euler: solve M * C^{n+1} = (M - dt*K) C^n + dt*f
    // Operator is just M -- doesn't depend on dt at all.
    prec = new mfem::GSSmoother(*M_mat);
    solver.SetOperator(*M_mat);
    solver.SetPreconditioner(*prec);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(0.0);
    solver.SetMaxIter(500);
    solver.SetPrintLevel(0);

    // Explicit RHS matrix, built once since dt is fixed:
    Tmat = Add(1.0, *M_mat, -dt, *K_mat);   // Tmat = M - dt*K  (note: MINUS now)
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
    mfem::FunctionCoefficient r2_coeff([](const mfem::Vector &x) {
        return x(0) * x(0);
    });

    mfem::LinearForm vol_lf(fespace);
    vol_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(r2_coeff));
    vol_lf.Assemble();

    double integral_Cr2 = vol_lf(C);
    double volume = R * R * R / 3.0;   // ∫_0^R r^2 dr

    return integral_Cr2 / volume;
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
    C.Save((prefix + "_concentration.gf").c_str());
}

void SphericalDiffusion::SaveMesh(const std::string &prefix) {
    mesh->Save((prefix + "_mesh.mesh").c_str());
}