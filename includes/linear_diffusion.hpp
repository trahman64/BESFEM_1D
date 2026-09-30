#ifndef LINEAR_DIFFUSION_HPP
#define LINEAR_DIFFUSION_HPP

#include "mfem.hpp"

class LinearDiffusion {
public:
	LinearDiffusion(mfem::Mesh *mesh, mfem::FiniteElementSpace *fespace,
                 mfem::Array<int> &nbc_bdr,   // NEW: Neumann boundary marker, caller-supplied
                 double C0, double dt);

    ~LinearDiffusion();

    void SetWeightVector(double eps_sep, double eps_eld, double tau_sep, double tau_eld);
    void SetCoefficient(mfem::GridFunction &D);   // D: diffusivity field (e.g. De_gf)
    void BuildOperator();
    void UpdateOperator();   // call after D's VALUES change (same mesh/space)

    void Stepping(mfem::GridFunction &source);

    double GetMeanConcentration();
    void SaveMesh();
    void SaveConc(const std::string &filename = "elyConc.gf");
    mfem::GridFunction& GetConcentration();

private:
    mfem::Mesh *mesh;
    mfem::FiniteElementSpace *fespace;

    double C0, dt;

    mfem::GridFunction C;
    mfem::Vector rhs, X, C_prev;

    mfem::Vector epsilon_vector;    // porosity per region -- for M
    mfem::Vector weight_vector;     // eps/tau^2 per region -- scales D for K
    mfem::PWConstCoefficient *epsilon_coeff = nullptr;
    mfem::PWConstCoefficient *region_weight = nullptr;
    mfem::GridFunctionCoefficient *D_coeff = nullptr;
    mfem::ProductCoefficient *weight_eff_D = nullptr;

    mfem::BilinearForm *M = nullptr;
    mfem::BilinearForm *K = nullptr;
    mfem::SparseMatrix *M_mat, *K_mat;

    mfem::SparseMatrix *TmatR = nullptr;
    mfem::SparseMatrix *TmatL = nullptr;
    mfem::GSSmoother   *prec = nullptr;
    mfem::CGSolver solver;

    mfem::Array<int> nbc_bdr;
    mfem::LinearForm rxn_lf, vol_lf;
};

#endif