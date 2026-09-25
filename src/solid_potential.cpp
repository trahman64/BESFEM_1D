#include "../includes/solid_potential.hpp"

SolidPotential::SolidPotential(mfem::FiniteElementSpace *fespace,
                                mfem::Coefficient &effective_kappa,
                                mfem::Array<int> &ess_bdr)
    : StatPotential(fespace, effective_kappa, ess_bdr),
      Additional(fespace->GetTrueVSize()),
      Kappa(fespace)
{
    Additional = 0.0;
}

SolidPotential::~SolidPotential() {
    delete region_weight;
    delete kappa_coeff;
    delete effective_kappa_inner;
    delete weight_eff_kappa;
    // fespace/mesh are NOT owned here -- do not delete them
}

mfem::Coefficient& SolidPotential::Conductivity(double eps_s_sep, double eps_s_eld,
                                                  double tau_s_sep, double tau_s_eld,
                                                  double kappa_s)
{
    mfem::Mesh *mesh = fespace->GetMesh();

    weight_vector.SetSize(mesh->attributes.Max());
    weight_vector(0) = eps_s_sep / (tau_s_sep * tau_s_sep);   // separator
    weight_vector(1) = eps_s_eld / (tau_s_eld * tau_s_eld);   // electrode

    region_weight = new mfem::PWConstCoefficient(weight_vector);

    Kappa = kappa_s;   // uniform value; use ProjectCoefficient if it should vary spatially
    kappa_coeff = new mfem::GridFunctionCoefficient(&Kappa);

    effective_kappa_inner = new mfem::ProductCoefficient(*region_weight, *kappa_coeff);
    weight_eff_kappa = new mfem::ProductCoefficient(*region_weight, *effective_kappa_inner);

    return *weight_eff_kappa;
}

mfem::Vector& SolidPotential::GetAdditional() {
    return Additional;
}