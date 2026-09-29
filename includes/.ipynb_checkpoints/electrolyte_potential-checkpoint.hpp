#ifndef LIQUID_POTENTIAL_HPP
#define LIQUID_POTENTIAL_HPP

#include "mfem.hpp"
#include "linear_diffusion.hpp"
#include <string>

class LiquidPotential
{
private:
    mfem::Mesh *mesh;
    mfem::FiniteElementSpace *fespace;
    LinearDiffusion *diffusion;

    double a;
    double F;
    double R;
    double T;

    double epsilon_s;
    double epsilon_e;
    double tau_e;

    double D_plus;
    double D_minus;

    mfem::GridFunction phi_e;

public:
    LiquidPotential(
        mfem::Mesh *mesh_,
        mfem::FiniteElementSpace *fespace_,
        LinearDiffusion *diffusion_,
        double a_,
        double F_,
        double R_,
        double T_,
        double epsilon_s_,
        double epsilon_e_,
        double tau_e_,
        double D_plus_,
        double D_minus_
    );

    void Solve(double rxn);

    void Save(
        const std::string &filename = "liquid_potential.gf"
    );

    const mfem::GridFunction& GetPotential() const;
};

#endif