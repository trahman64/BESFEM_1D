#ifndef RADIAL_DIFFUSION_HPP
#define RADIAL_DIFFUSION_HPP

#include "mfem.hpp"
class Radial_Diffusion{
    private:
        mfem::Mesh *mesh;
        mfem::FiniteElementSpace *fespace;
        mfem::GridFunction C;
        mfem::BilinearForm *M;
        mfem::BilinearForm *K;
        mfem::LinearForm *R;
        mfem::SparseMatrix *M_mat;
        mfem::SparseMatrix *K_mat;
        mfem::Vector C_new;
        mfem::Vector KC;
        mfem::Vector MC;
        mfem::Vector rhs;

    public:
        Radial_Diffusion(mfem::Mesh *mesh, mfem::FiniteElementSpace *fec);
        void Stepping(double dt,double rxn);
        void Save();
        double GetSurfaceConc();
};

#endif