#include "mfem.hpp"
#include "../includes/cell_kinetics.hpp"
#include "../includes/potential_utils.hpp"   // assuming Compute_i0/Compute_pOCV live here, adjust if not

CellKinetics::CellKinetics(mfem::FiniteElementSpace *fespace_,
                            mfem::Array<int> &region2_dofs,
                            mfem::Array<int> &is_in_region2_,
                            StatPotential &solid_potential_,
                            StatPotential &liquid_potential_,
                            double aPv, double alpha_a, double alpha_c,
                            double Cst1, double F)
    : fespace(fespace_), is_in_region2(is_in_region2_),
      solid_potential(solid_potential_), liquid_potential(liquid_potential_),
      OCV_gf(fespace_), i0_gf(fespace_), phs_gf(fespace_), phl_gf(fespace_),
      cellRxn(fespace_, region2_dofs, aPv, alpha_a, alpha_c, Cst1, F)
{
    OCV_gf = 0.0;
    i0_gf = 0.0;
    phs_gf = 0.0;
    phl_gf = 0.0;
}

mfem::GridFunction& CellKinetics::Compute(mfem::GridFunction &Ce, mfem::GridFunction &Cp_surf) {
    i0_gf = Compute_i0(Cp_surf, is_in_region2);
    OCV_gf = Compute_OCV(Cp_surf, is_in_region2);
    phs_gf = solid_potential.GetPotential();
    phl_gf = liquid_potential.GetPotential();

    return cellRxn.Compute(Ce, Cp_surf, i0_gf, OCV_gf, phs_gf, phl_gf);
}