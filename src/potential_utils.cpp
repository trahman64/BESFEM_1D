#include "mfem.hpp"

// in a small header, e.g. potential_utils.hpp, or just above SolidPotential's definition
mfem::Vector ComputeWeightVector(mfem::FiniteElementSpace *fespace,
                                   double eps_s_sep, double eps_s_eld,
                                   double tau_s_sep, double tau_s_eld)
{
    mfem::Vector weight_vector(fespace->GetMesh()->attributes.Max());
    weight_vector(0) = eps_s_sep / (tau_s_sep * tau_s_sep);   // separator
    weight_vector(1) = eps_s_eld / (tau_s_eld * tau_s_eld);   // electrode
    return weight_vector;
}
