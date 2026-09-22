#ifndef SPHERICAL_DIFFUSION_HPP
#define SPHERICAL_DIFFUSION_HPP

#include "mfem.hpp"
#include <string>

class SphericalDiffusion {
public:
    // radius         : particle radius R
    // n_elements     : number of elements along the radial mesh [0, R]
    // diffusion_coef : D (assumed spatially constant here)
    // x_idx          : index/ID of this particle within the host domain
    // dt             : fixed timestep used for all Stepping() calls
    // fe_order       : polynomial order of the FE space (default linear)
    // C0             : initial (uniform) concentration
    SphericalDiffusion(double radius,
                        int n_elements,
                        double diffusion_coef,
                        int x_idx,
                        double dt,
                        int fe_order = 1,
                        double C0 = 0.0);

    ~SphericalDiffusion();

    // Advance one timestep with a prescribed surface flux
    // (D * dC/dr at r = R; positive = flux INTO the sphere)
    void Stepping(double surface_flux);

    mfem::GridFunction& GetConcentration();
    double GetMeanConcentration();        // volume-averaged C (r^2-weighted)
    double GetConcentrationAt(double x);  // value at a specific radial location x in [0, R]
    int    GetParticleID() const;         // FIX: was declared as returning GridFunction&

    void SaveConc(const std::string &prefix = "sphere");
    void SaveMesh(const std::string &prefix = "sphere");

private:
    double R;      // particle radius
    double D;       // diffusion coefficient
    int order;      // FE order
    int x_idx;       // index/ID of this particle within the host domain
    double dt;       // fixed timestep, set at construction

    mfem::Mesh *mesh;
    mfem::FiniteElementCollection *fec;
    mfem::FiniteElementSpace *fespace;

    mfem::GridFunction C;
    mfem::Vector C_prev;

    mfem::BilinearForm *M, *K;
    mfem::SparseMatrix *M_mat, *K_mat;

    mfem::Array<int> outer_bdr_marker;  // marks r = R boundary attribute

    // Cached forward-Euler operator/solver (built once, during construction)
    mfem::SparseMatrix *Tmat = nullptr;
    mfem::GSSmoother   *prec = nullptr;
    mfem::CGSolver      solver;

    void BuildMesh(int n_elements);
    void BuildOperators();
    void Initialize();   // builds Tmat + solver, using this->dt -- called once, in the constructor
};

#endif