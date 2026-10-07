#include "../includes/reaction.hpp"

#include <cmath>


ButlerVolmer::ButlerVolmer(
    mfem::FiniteElementSpace *fespace_,
    double i0_,
    double F_,
    double R_,
    double T_,
    double alpha_,
    double a_
)
    :
    fespace(fespace_),

    i0(i0_),
    F(F_),
    R(R_),
    T(T_),
    alpha(alpha_),
    a(a_),

    OCV(fespace_),
    Kfw(fespace_),
    Kbw(fespace_),
    dPHE(fespace_),
    rxn(fespace_)
{
    OCV = 0.0;
    Kfw = 0.0;
    Kbw = 0.0;
    dPHE = 0.0;
    rxn = 0.0;
}



double ButlerVolmer::CalculateOCV(double x)
{
    return 1.095 * x * x
         - 8.234e-7 * std::exp(14.32 * x)
         + 4.692 * std::exp(-0.5389 * x);
}


void ButlerVolmer::Compute(
    const mfem::GridFunction &Ce,
    const mfem::GridFunction &Cp_surf,
    const mfem::GridFunction &phi_e,
    const mfem::GridFunction &phi_p
)
{
    OCV = 0.0;
    Kfw = 0.0;
    Kbw = 0.0;
    dPHE = 0.0;
    rxn = 0.0;


    // F / RT
    double Cst1 = F / (R * T);




    for (int i = 20; i < 80; i++)
    {
        double cp = Cp_surf(i);


        

        OCV(i) = CalculateOCV(cp);



        Kfw(i) =
            i0
            /
            (F * 0.001)
            *
            std::exp(
                alpha
                * Cst1
                * OCV(i)
            );


        Kbw(i) =
            i0
            /
            (F * cp)
            *
            std::exp(
                -alpha
                * Cst1
                * OCV(i)
            );



        dPHE(i) =
            phi_p(i) - phi_e(i);



        rxn(i) =
            a
            *
            (
                Kfw(i)
                * Ce(i)
                * std::exp(
                    -alpha
                    * Cst1
                    * dPHE(i)
                )

                -

                Kbw(i)
                * cp
                * std::exp(
                    alpha
                    * Cst1
                    * dPHE(i)
                )
            );
    }
}




mfem::GridFunction& ButlerVolmer::GetReaction()
{
    return rxn;
}