#include "mfem.hpp"
#include "../includes/butler_volmer.hpp"
#include <cmath>

ButlerVolmer::ButlerVolmer(mfem::FiniteElementSpace *fespace_,
                            mfem::Array<int> &region2_dofs_,
                            double alpha_a_, double alpha_c_,
                            double Cst1_)
    : fespace(fespace_), region2_dofs(region2_dofs_),
      alpha_a(alpha_a_), alpha_c(alpha_c_), Cst1(Cst1_),
      rxn(fespace_)
{
    rxn = 0.0;
}

mfem::GridFunction& ButlerVolmer::Compute(mfem::GridFunction &i0,
                                           mfem::GridFunction &OCV,
                                           mfem::GridFunction &phs,
                                           mfem::GridFunction &phl)
{
    rxn = 0.0;

    for (int i = 0; i < region2_dofs.Size(); i++) {
        int p_id = region2_dofs[i];
        double eta = phs(p_id) - phl(p_id) - OCV(p_id);
        rxn(p_id) = i0(p_id) * (std::exp(-alpha_a * Cst1 * eta)
                              - std::exp( alpha_c * Cst1 * eta));
    	std::cout << rxn(p_id) << " " << p_id << std::endl; 
    }
    return rxn;
}

double ButlerVolmer::GetTotalRxnCurrent()
{
	return rxn.Sum();
}

void ButlerVolmer::Save(const std::string &filename)
{
	rxn.Save(filename.c_str());
}