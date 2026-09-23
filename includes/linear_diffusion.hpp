#ifndef LINEAR_DIFFUSION_HPP
#define LINEAR_DIFFUSION_HPP

#include "mfem.hpp"

class LinearDiffusion {
public:
    LinearDiffusion(mfem::Mesh *mesh,
              mfem::FiniteElementSpace *fespace,
              double epsilon_sep,
              double epsilon_eld,              
              double De,
              double tau_sep,
              double tau_eld,
              double t_minus,
              double C0,
              double dt);

    ~LinearDiffusion();

    void Stepping(mfem::GridFunction source);
    void SaveMesh();
    void SaveConc();
    mfem::GridFunction& GetConcentration();
    double GetMeanConcentration();

private:
    mfem::Mesh *mesh;
    mfem::FiniteElementSpace *fespace;

	double epsilon_sep;
	double epsilon_eld;             
	double De;
	double tau_sep;
	double tau_eld;
	double t_minus;
	double C0;
	double dt;

    mfem::GridFunction C;
    mfem::Vector rhs, X;
    mfem::Vector C_prev;

    mfem::BilinearForm *M, *K;
    mfem::SparseMatrix *M_mat, *K_mat;  

    mfem::SparseMatrix *TmatR = nullptr;
    mfem::SparseMatrix *TmatL = nullptr;
    mfem::GSSmoother   *prec = nullptr;
    mfem::CGSolver     solver;    
    
    mfem::Vector epsilon_vector;              // stored so PWConstCoefficient stays valid
    mfem::PWConstCoefficient *epsilon_coeff;  // built once, reused    

    mfem::LinearForm rxn_lf, vol_lf;          // NEW: assembled once in constructor, reused later
};

#endif