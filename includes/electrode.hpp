#ifndef ELECTRODE_HPP
#define ELECTRODE_HPP

#include "mfem.hpp"
#include "spherical_diffusion.hpp"
#include <vector>

class Electrode {
public:
    Electrode(mfem::FiniteElementSpace *fespace,
              mfem::Array<int> &region2_dofs,
              mfem::GridFunction &parti_radii,
              int n_elements, int fe_order,
              double Cp0, double dt);

    void Stepping(mfem::GridFunction &flux_gf);
    void UpdateOperators();

    mfem::GridFunction& GetSurfaceConcentration();
    mfem::GridFunction& GetPartiMeanConcentration();
    mfem::GridFunction& GetParticleVolume();
    mfem::GridFunction& GetParticleTotalLi();

    double GetTotalVolume();
    double GetTotalLi();
    double GetDoD();
    
    void SaveAllConc(const std::string &prefix = "particle");
    void SaveConcByID(int particle_id, const std::string &prefix = "particle");

private:
    mfem::FiniteElementSpace *fespace;
    mfem::Array<int> region2_dofs;
    mfem::GridFunction parti_radii;
    std::vector<SphericalDiffusion> particles;

    mfem::GridFunction Cp_surf, Cp_mConc, part_totLi, part_volume;
};

#endif