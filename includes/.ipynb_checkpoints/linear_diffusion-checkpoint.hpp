#ifndef LINEAR_DIFFUSION_HPP
#define LINEAR_DIFFUSION_HPP

#include "mfem.hpp"

class LinearDiffusion{
    private:
        mfem::Mesh *mesh;
        mfem::FiniteElementSpace *fespace;

        double epsilon;
        double D;
        double a;
        double t_minus;
        double tau_electrode;
        double C0;
        double electrode_length;
        double dt;
        mfem::GridFunction C;
        mfem::Vector C_prev;
        mfem::Vector KC;
        mfem::Vector MC;
        mfem::Vector rhs;
        mfem::PWConstCoefficient *epsilon_coeff;

        mfem::BilinearForm *M;
        mfem::BilinearForm *K;
    
        mfem::SparseMatrix *M_mat;
        mfem::SparseMatrix *K_mat;
    
        mfem::SparseMatrix *TmatR = nullptr;
        mfem::SparseMatrix *TmatL = nullptr;
    
        // Solver
        mfem::GSSmoother *prec = nullptr;
        mfem::CGSolver solver;
        mfem::Vector epsilon_vector;

    public:
        LinearDiffusion(mfem::Mesh *mesh_,mfem::FiniteElementSpace *fespace_,
              double epsilon,
              double D,
              double a,
              double t_minus,
              double tau_electrode,
              double C0,
              double dt);

        ~LinearDiffusion();
        void Stepping(double rxn);
        double GetMeanConcentration();
        void SaveMesh();
        void SaveConc();
        mfem::GridFunction& GetConcentration();
};

#endif
    