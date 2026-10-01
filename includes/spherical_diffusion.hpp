#ifndef SPHERICAL_DIFFUSION_HPP
#define SPHERICAL_DIFFUSION_HPP

#include "mfem.hpp"
#include <string>

class SphericalDiffusion {
public:
    SphericalDiffusion(int x_idx,
                        double radius,
                        int n_elements,
                        int fe_order,
                        double C0,
                        double dt);

    ~SphericalDiffusion();

    void Stepping(double surface_flux);
    void UpdateOperator();   // call after Stepping() updates C, to refresh D_li-dependent K/Tmat

    mfem::GridFunction& GetConcentration();
    double GetMeanConcentration();
    double GetConcentrationAt(double x);
    int GetParticleID() const;

    void SaveConc(const std::string &prefix = "sphere");
    void SaveMesh(const std::string &prefix = "sphere");

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

    mfem::GridFunction D_li;
    mfem::GridFunction r;                                 // now a member -- referenced by r2_coeff
    mfem::GridFunctionCoefficient *r2_coeff = nullptr;     // now a member
    mfem::GridFunctionCoefficient *D_li_coeff = nullptr;   // now a member (was "DC")
    mfem::ProductCoefficient *D_r2 = nullptr;              // now a member

    mfem::Array<int> outer_bdr_marker;

    mfem::BilinearForm *M = nullptr;
    mfem::BilinearForm *K = nullptr;
    mfem::SparseMatrix *M_mat, *K_mat;
    mfem::SparseMatrix *Tmat = nullptr;

    mfem::GSSmoother *prec = nullptr;
    mfem::CGSolver solver;

    void BuildMesh(int n_elements);
    void BuildOperators();
//     void Initialize();
};

#endif