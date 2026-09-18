#include "mfem.hpp"
#include "../includes/linear_diffusion.hpp"
#include <iostream>

LinearDiffusion::LinearDiffusion(mfem::Mesh *mesh_, mfem::FiniteElementSpace *fespace_,
                      double epsilon_, double D_, double a_, double t_minus_,
                      double tau_electrode_, double C0_, double dt_)
    : mesh(mesh_), fespace(fespace_),
      epsilon(epsilon_), D(D_), a(a_), t_minus(t_minus_), tau_electrode(tau_electrode_), C0(C0_), dt(dt_),
      C(fespace_), MC(fespace_->GetTrueVSize()), KC(fespace_->GetTrueVSize()),
      rhs(fespace_->GetTrueVSize()),
      C_prev(fespace_->GetTrueVSize())
{
    mfem::Vector epsilon_vector(mesh->attributes.Max());
    mfem::Vector diffusion_vector(mesh->attributes.Max());

    const double tau_sep = 1.0;

    // Region 1 (separator) -- fully permeable
    epsilon_vector(0) = 1.0;
    diffusion_vector(0) = D * epsilon_vector(0) / (tau_sep * tau_sep);

    // Region 2 (electrode)
    epsilon_vector(1) = epsilon;
    diffusion_vector(1) = epsilon * D / (tau_electrode * tau_electrode);

//     mfem::PWConstCoefficient epsilon_coeff(epsilon_vector);
    epsilon_coeff = new mfem::PWConstCoefficient(epsilon_vector); 
    mfem::PWConstCoefficient diffusion_coeff(diffusion_vector);

    // M
    M = new mfem::BilinearForm(fespace);
    M->AddDomainIntegrator(new mfem::MassIntegrator(*epsilon_coeff));
    M->Assemble();
    M->Finalize();

    // K
    K = new mfem::BilinearForm(fespace);
    K->AddDomainIntegrator(new mfem::DiffusionIntegrator(diffusion_coeff));
    K->Assemble();
    K->Finalize();

    M_mat = &M->SpMat();
    K_mat = &K->SpMat();

    // Compute electrode (region 2) length from the mesh
    mfem::Array<int> electrode_marker(mesh->attributes.Max());
    electrode_marker = 0;
    electrode_marker[1] = 1;   // attribute 2 -> electrode

    mfem::ConstantCoefficient one(1.0);
    mfem::LinearForm length_lf(fespace);
    length_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(one), electrode_marker);
    length_lf.Assemble();
    electrode_length = length_lf.Sum();

    // Initial condition
    C = C0;
    C.GetTrueDofs(C_prev);
    
    // Build Crank-Nicolson matrices ONCE, now that dt, M_mat, K_mat are all fixed
	// Crank-Nicolson: (M + 0.5*dt*K) C^{n+1} = (M - 0.5*dt*K) C^n + dt*f
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
}

void LinearDiffusion::Stepping(double rxn) {
    mfem::Vector reaction_vector(mesh->attributes.Max());
    reaction_vector = 0.0;

    // Separator
    reaction_vector(0) = 0.0;

    // Electrode
    reaction_vector(1) = -1.0 * (a * rxn * t_minus);

    mfem::PWConstCoefficient reaction(reaction_vector);
    mfem::LinearForm R_current(fespace);
    R_current.AddDomainIntegrator(new mfem::DomainLFIntegrator(reaction));

    double f_in = a * rxn * t_minus * electrode_length;

    mfem::ConstantCoefficient nbcCoef(f_in);

    mfem::Array<int> nbc_w_bdr(mesh->bdr_attributes.Max());
    nbc_w_bdr = 0;
    nbc_w_bdr[0] = 1;   // Neumann (flux) BC on the west boundary

    R_current.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(nbcCoef), nbc_w_bdr);
    R_current.Assemble();

    rhs = R_current;
    rhs *= dt;

    mfem::Vector X(fespace->GetTrueVSize());
    TmatR->Mult(C_prev, X);
    X += rhs;

    solver.Mult(X, C_prev);

    C.SetFromTrueDofs(C_prev);
}


double LinearDiffusion::GetMeanConcentration() {
    mfem::LinearForm vol_lf(fespace);
    vol_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(*epsilon_coeff));
    vol_lf.Assemble();

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