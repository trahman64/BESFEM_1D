#ifndef CELL_KINETICS_HPP
#define CELL_KINETICS_HPP

#include "mfem.hpp"
#include "stat_potential.hpp"
#include "butler_volmer.hpp"

class CellKinetics {
public:
    CellKinetics(mfem::FiniteElementSpace *fespace,
                 mfem::Array<int> &region2_dofs,
                 mfem::Array<int> &is_in_region2,
                 StatPotential &solid_potential,
                 StatPotential &liquid_potential,
                 double aPv, double alpha_a, double alpha_c,
                 double Cst1, double F);

    // Computes i0, OCV from Cp_surf, pulls current phs/phl from the potential
    // solvers, and returns the resulting reaction rate field.
    mfem::GridFunction& Compute(mfem::GridFunction &Ce, mfem::GridFunction &Cp_surf);

    mfem::GridFunction& GetI0() { return i0_gf; }
    mfem::GridFunction& GetOCV() { return OCV_gf; }
	mfem::GridFunction& GetPhs() { return phs_gf; }
	mfem::GridFunction& GetPhl() { return phl_gf; } 
	double GetTotalRxnCurrent() { return cellRxn.GetTotalRxnCurrent(); }	   

private:
    mfem::FiniteElementSpace *fespace;
    mfem::Array<int> is_in_region2;
    StatPotential &solid_potential;
    StatPotential &liquid_potential;

    mfem::GridFunction OCV_gf, i0_gf, phs_gf, phl_gf;
    ButlerVolmer cellRxn;
};

#endif