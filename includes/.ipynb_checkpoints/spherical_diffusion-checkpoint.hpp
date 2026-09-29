#ifndef SPHERICAL_DIFFUSION_HPP
#define SPHERICAL_DIFFUSION_HPP

#include "mfem.hpp"
class SphericalDiffusion{
    private:
        double R;
        double D;
        int order;
        int x_idx;
        double dt;

        mfem::Mesh *mesh;
        mfem::FiniteElementCollection *fec;
        mfem::FiniteElementSpace *fespace;

        mfem::GridFunction C;
        mfem::Vector C_prev;

        mfem::BilinearForm *M, *K;
        mfem::SparseMatrix *M_mat, *K_mat;

        mfem::Array<int> outer_bdr_marker;
        mfem::SparseMatrix *Tmat = nullptr;
        mfem::GSSmoother *prec = nullptr;
        mfem::CGSolver solver;

        void BuildMesh(int n_elements);
        void BuildOperators();
        void Initialize();


    public:
        SphericalDiffusion(double radius, int n_elements, double diffusion_coef, int x_idx, double dt, int fe_order = 1, double C0 = 0.0);
        ~SphericalDiffusion();
        
        void Stepping(double surface_flux);
        double GetMeanConcentration();
        double GetConcentrationAt(double x);
        int GetParticleID() const;

        void SaveConc(const std::string &prefix = "sphere");
        void SaveMesh(const std::string &prefix = "sphere");
};

#endif