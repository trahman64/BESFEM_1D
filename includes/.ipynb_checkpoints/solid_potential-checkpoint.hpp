#ifndef SOLID_POTENTIAL_HPP
#define SOLID_POTENTIAL_HPP

#include "mfem.hpp"
#include <string>

class SolidPotential
{
private:
    mfem::Mesh *mesh;
    mfem::FiniteElementSpace *fespace;

    mfem::GridFunction phi_s;

public:
    SolidPotential(mfem::Mesh *mesh_,
                   mfem::FiniteElementSpace *fespace_);

    void Solve(double rxn);

    void Save(const std::string &filename = "solid_potential.gf");

    const mfem::GridFunction& GetPotential() const;
};

#endif