#include "mfem.hpp"
#include "../includes/diffusion.hpp"
#include "../includes/diffusion.hpp"
#include <iostream>


Diffusion::Diffusion(mfem::Mesh *mesh,mfem::FiniteElementSpace *fec):mesh(mesh),fespace(fec),
C(fec),MC(fec->GetTrueVSize()),KC(fec->GetTrueVSize()),rhs(fec->GetTrueVSize()),
C_prev(fec->GetTrueVSize()),C_new(fec->GetTrueVSize()){

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
    reaction_vector(1) = -1.0 * (a_rxn * rxn * t_minus);

    mfem::PWConstCoefficient reaction(reaction_vector);
    mfem::LinearForm R_current(fespace);
    R_current.AddDomainIntegrator(new mfem::DomainLFIntegrator(reaction));

    double electrode_length = 60e-4 ;    //0.5;
    double f_in = a_rxn * rxn * t_minus * electrode_length; 

 
	mfem::ConstantCoefficient nbcCoef(f_in);

	mfem::Array<int> boundary_dofs;					// nature boundary	
	// Neumann BC on the west boundary. CnE
	mfem::Array<int> nbc_w_bdr(mesh->bdr_attributes.Max());
	nbc_w_bdr = 0; nbc_w_bdr[0] = 1;	
	R_current.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(nbcCoef), nbc_w_bdr);	
	R_current.Assemble();

	mfem::Vector X, rhs;
	mfem::SparseMatrix A;
	K->FormLinearSystem(boundary_dofs, C, R_current, A, X, rhs);
	rhs *= dt;	
	
	
	mfem::SparseMatrix *TmatR, *TmatL;
	// Crank-Nicolson matrices
	TmatR = Add(1.0, *M_mat, -0.5*dt, *K_mat);		
	TmatL = Add(1.0, *M_mat,  0.5*dt, *K_mat);		
			
	
	
	TmatR->Mult(C_prev, X);
	X += rhs;
	
	
	// solver
    mfem::GSSmoother M_prec(*TmatL);
    mfem::CGSolver solver;

    solver.SetOperator(*TmatL);
    solver.SetPreconditioner(M_prec);
    solver.SetRelTol(1e-12);
    solver.SetAbsTol(0.0);
    solver.SetMaxIter(500);
    solver.SetPrintLevel(0);		
	
	// time stepping
	solver.Mult(X, C_prev) ;
	
	// recover
	C.SetFromTrueDofs(C_prev);  
	
}
void Diffusion::Save(){
    mesh->Save("diffusion_mesh.mesh");
    C.Save("diffusion_concentration.gf");
}

mfem::GridFunction& Diffusion::GetConcentration(){
    return C;
}


    
    
