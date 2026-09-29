#ifndef STAT_POTENTIAL_HPP
#define STAT_POTENTIAL_HPP

#include "mfem.hpp"

class StatPotential {
public:
    StatPotential(
            mfem::FiniteElementSpace *fespace,
            mfem::Array<int> &ess_bdr
        );
    ~StatPotential();

    void Solve(mfem::GridFunction &source, mfem::Vector &Additional, double Bv);

    mfem::GridFunction& GetPotential();
    void SavePote();
    void SetWeightVector(double eps_sep, double eps_eld, double tau_sep,double tau_eld);

    void SetCoefficient(mfem::GridFunction &Kappa);
    void BuildOperator();
    void UpdateOperator();

private:
    mfem::FiniteElementSpace *fespace;

    mfem::GridFunction phi;
    mfem::Vector X, B;

    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;

    mfem::BilinearForm *K;
    mfem::SparseMatrix *K_mat;
    mfem::SparseMatrix A;   // FIX: object, not pointer

    mfem::GSSmoother *prec;
    mfem::CGSolver solver;
    mfem::PWConstCoefficient *region_weight = nullptr;
    mfem::GridFunctionCoefficient *kappa_coeff = nullptr;
    mfem::ProductCoefficient *weight_eff_kappa = nullptr;
    mfem::Vector weight_vector;
};

#endif