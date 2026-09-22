#ifndef LINEAR_DIFFUSION_HPP
#define LINEAR_DIFFUSION_HPP

#include "mfem.hpp"

class LinearDiffusion {
public:
    LinearDiffusion(mfem::Mesh *mesh,
              mfem::FiniteElementSpace *fespace,
              double epsilon,
              double D,
              double a,
              double t_minus,
              double tau_electrode,
              double C0,
              double dt);

    ~LinearDiffusion();

    void Stepping(double rxn);
    void SaveMesh();
    void SaveConc();
    mfem::GridFunction& GetConcentration();
    double GetMeanConcentration();

private:
    mfem::Mesh *mesh;
    mfem::FiniteElementSpace *fespace;

    double epsilon;
    double D;
    double a;
    double t_minus;
    double tau_electrode;
    double C0;
    double electrode_length;   // computed from mesh (region 2)
    double dt;  // fixed timestep, set at construction

    mfem::GridFunction C;
    mfem::Vector MC, KC, rhs;
    mfem::Vector C_prev;

    mfem::BilinearForm *M, *K;
    mfem::SparseMatrix *M_mat, *K_mat;  

    mfem::SparseMatrix *TmatR = nullptr;
    mfem::SparseMatrix *TmatL = nullptr;
    mfem::GSSmoother   *prec = nullptr;
    mfem::CGSolver      solver;    
    
    mfem::Vector epsilon_vector;              // stored so PWConstCoefficient stays valid
    mfem::PWConstCoefficient *epsilon_coeff;  // built once, reused    
    
};

#endif