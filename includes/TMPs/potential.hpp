#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

#include "mfem.hpp"

class Diffusion;
class Potential{
    private:
        mfem::Mesh *mesh;
        mfem::FiniteElementSpace *fespace;
        mfem::GridFunction phi_e;
        mfem::GridFunction phi_s;
        void SolveSolid(double rxn);
        void SolveElectrolyte(double rxn);
        Diffusion *diffusion;  

    public:
        Potential(mfem::Mesh *mesh, mfem::FiniteElementSpace *fec, Diffusion *diffusion);
        void Solve(double rxn);
        void Save();
        const mfem::GridFunction& GetSolidPotential() const;
        const mfem::GridFunction& GetLiquidPotential() const;
};

#endif