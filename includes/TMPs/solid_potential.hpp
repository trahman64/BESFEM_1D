#ifndef SOLID_POTENTIAL_HPP
#define SOLID_POTENTIAL_HPP

#include "stat_potential.hpp"

class SolidPotential : public StatPotential {
public:
    SolidPotential(mfem::FiniteElementSpace *fespace,
                    mfem::Coefficient &effective_kappa,
                    mfem::Array<int> &ess_bdr);

    ~SolidPotential();

    mfem::Coefficient& Conductivity(double eps_s_sep, double eps_s_eld,
                                     double tau_s_sep, double tau_s_eld,
                                     double kappa_s);

    mfem::Vector& GetAdditional();

private:
    mfem::Vector Additional;

    // Persistent storage for the coefficient chain built by Conductivity() --
    // ALL of these must outlive weight_eff_kappa, since ProductCoefficient
    // only stores references to its inputs.
    mfem::Vector weight_vector;
    mfem::PWConstCoefficient *region_weight = nullptr;
    mfem::GridFunction Kappa;
    mfem::GridFunctionCoefficient *kappa_coeff = nullptr;
    mfem::ProductCoefficient *effective_kappa_inner = nullptr;
    mfem::ProductCoefficient *weight_eff_kappa = nullptr;
};

#endif