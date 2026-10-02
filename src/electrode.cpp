#include "../includes/electrode.hpp"

Electrode::Electrode(mfem::FiniteElementSpace *fespace_,
                      mfem::Array<int> &region2_dofs_,
                      mfem::GridFunction &parti_radii_, 
                      int n_elements, int fe_order,
                      double Cp0, double dt)
    : fespace(fespace_), region2_dofs(region2_dofs_), parti_radii(parti_radii_),
      Cp_surf(fespace_), Cp_mConc(fespace_), part_totLi(fespace_), part_volume(fespace_)
{
    int n_elde_nodes = region2_dofs.Size();
    particles.reserve(n_elde_nodes);

    for (int i = 0; i < n_elde_nodes; i++) {
        int p_id = region2_dofs[i];
        particles.emplace_back(p_id, parti_radii(p_id), n_elements, fe_order, Cp0, dt);
    }

    Cp_surf = 0.0;
    Cp_mConc = 0.0;
    part_totLi = 0.0;
    part_volume = 0.0;

    // Volume doesn't change over time -- compute once, here
    for (int p = 0; p < n_elde_nodes; p++) {
        int p_id = particles[p].GetParticleID();
        part_volume(p_id) = particles[p].GetParticleVolume();
    }
}

void Electrode::Stepping(mfem::GridFunction &flux_gf) {
    for (auto &particle : particles) {
        int p_id = particle.GetParticleID();
        particle.Stepping(flux_gf(p_id));
    }
//     flux_gf.Print();
}

void Electrode::UpdateOperators() {
    for (auto &particle : particles) {
        particle.UpdateOperator();
    }
}

mfem::GridFunction& Electrode::GetSurfaceConcentration() {
    for (auto &particle : particles) {
        int p_id = particle.GetParticleID();
        Cp_surf(p_id) = particle.GetConcentrationAt(parti_radii(p_id));
    }
    return Cp_surf;
}

mfem::GridFunction& Electrode::GetPartiMeanConcentration() {
    for (auto &particle : particles) {
        int p_id = particle.GetParticleID();
        Cp_mConc(p_id) = particle.GetMeanConcentration();
    }
    return Cp_mConc;
}

mfem::GridFunction& Electrode::GetParticleVolume() {
    return part_volume;   // already computed once in the constructor
}

mfem::GridFunction& Electrode::GetParticleTotalLi() {
    for (auto &particle : particles) {
        int p_id = particle.GetParticleID();
        part_totLi(p_id) = particle.GetParticleTotalLi();
    }
    return part_totLi;
}

double Electrode::GetTotalVolume() {
    return part_volume.Sum();
}

double Electrode::GetTotalLi() {
    GetParticleTotalLi();   // refresh part_totLi first
    return part_totLi.Sum();
}

double Electrode::GetDoD() {
    GetParticleTotalLi();   // refresh part_totLi first
    return part_totLi.Sum()/part_volume.Sum();
}

void Electrode::SaveAllConc(const std::string &prefix) {
    for (auto &particle : particles) {
        particle.SaveConc(prefix);
    }
}

void Electrode::SaveConcByID(int particle_id, const std::string &prefix) {
    for (auto &particle : particles) {
        if (particle.GetParticleID() == particle_id) {
            particle.SaveConc(prefix);
            return;
        }
    }
    std::cerr << "SaveConcByID: no particle with ID " << particle_id << " found" << std::endl;
}