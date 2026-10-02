#include "mfem.hpp"
#include "../includes/butler_volmer.hpp"
#include <cmath>

ButlerVolmer::ButlerVolmer(mfem::FiniteElementSpace *fespace_,
                            mfem::Array<int> &region2_dofs_,
                            double aPv_, double alpha_a_, double alpha_c_,
                            double Cst1_, double F_)
    : fespace(fespace_), region2_dofs(region2_dofs_), aPv(aPv_),
      alpha_a(alpha_a_), alpha_c(alpha_c_), Cst1(Cst1_), F(F_),
      rxn(fespace_), rxn_lf(fespace_)
{
    rxn = 0.0;

    mfem::Array<int> region2_marker(fespace->GetMesh()->attributes.Max());
    region2_marker = 0;
    region2_marker[1] = 1;

    mfem::ConstantCoefficient one(1.0);
    rxn_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(one), region2_marker);
    rxn_lf.Assemble();
}

mfem::GridFunction& ButlerVolmer::Compute(mfem::GridFunction &Ce_gf,
										  mfem::GridFunction &Cp_gf,
										  mfem::GridFunction &i0,
                                          mfem::GridFunction &OCV,
                                          mfem::GridFunction &phs,
                                          mfem::GridFunction &phl)
{
    rxn = 0.0;

    for (int i = 0; i < region2_dofs.Size(); i++) {
        int p_id = region2_dofs[i];

        double eta = phs(p_id) - phl(p_id);
        double Kfw = i0(p_id)/(F*0.001      ) * std::exp( alpha_c*Cst1*OCV(p_id));
        double Kfb = i0(p_id)/(F*Cp_gf(p_id)) * std::exp(-alpha_a*Cst1*OCV(p_id));        
        
        rxn(p_id) = Kfw*Ce_gf(p_id) * std::exp(-alpha_a*Cst1*eta) -
                    Kfb*Cp_gf(p_id) * std::exp( alpha_c*Cst1*eta) ;   	
        rxn(p_id) *= aPv;            	
    }
    return rxn;
}

double ButlerVolmer::GetTotalRxnCurrent()
{
	return rxn_lf(rxn);   // ∫_region2 rxn dx
// 	return rxn.Sum();
}

void ButlerVolmer::Save(const std::string &filename)
{
	rxn.Save(filename.c_str());
}