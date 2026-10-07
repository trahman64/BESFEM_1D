#ifndef REACTION_HPP
#define REACTION_HPP

#include "mfem.hpp"

class ButlerVolmer
{
public:

    ButlerVolmer(
        mfem::FiniteElementSpace *fespace,
        double i0,
        double F,
        double R,
        double T,
        double alpha,
        double a
    );


    void Compute(
        const mfem::GridFunction &Ce,
        const mfem::GridFunction &Cp_surf,
        const mfem::GridFunction &phi_e,
        const mfem::GridFunction &phi_p
    );


    mfem::GridFunction& GetReaction();


private:

    mfem::FiniteElementSpace *fespace;

    // Constants
    double i0;
    double F;
    double R;
    double T;
    double alpha;
    double a;


    mfem::GridFunction OCV;
    mfem::GridFunction Kfw;
    mfem::GridFunction Kbw;
    mfem::GridFunction dPHE;
    mfem::GridFunction rxn;


    double CalculateOCV(double x);
};

#endif