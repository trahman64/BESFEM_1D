#ifndef BUTLER_VOLMER_HPP
#define BUTLER_VOLMER_HPP

#include "mfem.hpp"

class ButlerVolmer {
public:
    ButlerVolmer(mfem::FiniteElementSpace *fespace,
                 mfem::Array<int> &region2_dofs,
                 double alpha_a, double alpha_c,
                 double Cst1, double F);

    mfem::GridFunction& Compute(mfem::GridFunction &Ce_gf,
    						    mfem::GridFunction &Cp_gf,
    							mfem::GridFunction &i0,
                                mfem::GridFunction &OCV,
                                mfem::GridFunction &phs,
                                mfem::GridFunction &phl);
    void Save(const std::string &filename = "Rxn.gf");   
    
    double GetTotalRxnCurrent();                           

private:
    mfem::FiniteElementSpace *fespace;
    mfem::Array<int> region2_dofs;
    double alpha_a, alpha_c, Cst1, F;
    mfem::GridFunction rxn;
};

#endif