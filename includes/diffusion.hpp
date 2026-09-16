#ifndef DIFFUSION_HPP
#define DIFFUSION_HPP

#include "mfem.hpp"

class Diffusion{
    private:
        mfem::Mesh *mesh;
        mfem::FiniteElementSpace *fespace;
        mfem::GridFunction C;
        mfem::BilinearForm *M;
        mfem::BilinearForm *K;
        mfem::SparseMatrix *M_mat;
        mfem::SparseMatrix *K_mat;
        mfem::Vector C_prev;
        mfem::Vector C_new;
        mfem::Vector KC;
        mfem::Vector MC;
        mfem::Vector rhs;

    public:
        Diffusion(mfem::Mesh *mesh, mfem::FiniteElementSpace *fec);
        void Stepping(double dt,double rxn);
        void Save();
        mfem::GridFunction& GetConcentration();
};

#endif
    