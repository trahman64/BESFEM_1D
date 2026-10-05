#include "../includes/electrode.hpp"

Electrode::Electrode(mfem::FiniteElementSpace *fespace_,
                      mfem::Array<int> &region2_dofs_,                    
                      mfem::GridFunction &parti_radii_, 
                      int n_elements, int fe_order,
                      mfem::GridFunction &parti_initC_, double dt)
    : fespace(fespace_), region2_dofs(region2_dofs_), 
      parti_initC(parti_initC_), parti_radii(parti_radii_),
      Cp_surf(fespace_), Cp_mConc(fespace_), part_totLi(fespace_), 
      part_volume(fespace_), part_surfArea(fespace_)
{
    int n_elde_nodes = region2_dofs.Size();

    particles.reserve(n_elde_nodes);

    for (int i = 0; i < n_elde_nodes; i++) {
        int p_id = region2_dofs[i];
        particles.push_back(std::make_unique<SphericalDiffusion>(p_id,
            parti_radii(p_id), n_elements, fe_order, 
            parti_initC(p_id), dt));
    }

    Cp_surf = 0.0;
    Cp_mConc = 0.0;
    part_totLi = 0.0;
    part_volume = 0.0;
    part_surfArea = 0.0;

    for (int p = 0; p < n_elde_nodes; p++) {
        int p_id = particles[p]->GetParticleID();
        part_volume(p_id) = particles[p]->GetParticleVolume();
        part_surfArea(p_id) = particles[p]->GetParticleSurfaceArea();
    }
}

void Electrode::Stepping(mfem::GridFunction &flux_gf) {
    for (auto &particle : particles) {
        int p_id = particle->GetParticleID();
        particle->Stepping(flux_gf(p_id));
    }
}

void Electrode::UpdateOperators() {
    for (auto &particle : particles) {
        particle->UpdateOperator();
    }
}

mfem::GridFunction& Electrode::GetSurfaceConcentration() {
    Cp_surf = 0.0;	
    for (auto &particle : particles) {
        int p_id = particle->GetParticleID();
        Cp_surf(p_id) = particle->GetConcentrationAt(parti_radii(p_id));
    }
    return Cp_surf;
}

mfem::GridFunction& Electrode::GetPartiMeanConcentration() {
    Cp_mConc = 0.0;
    for (auto &particle : particles) {
        int p_id = particle->GetParticleID();
        Cp_mConc(p_id) = particle->GetMeanConcentration();
    }
    return Cp_mConc;
}

mfem::GridFunction& Electrode::GetParticleVolume() {
    return part_volume;   // already computed once in the constructor
}

mfem::GridFunction& Electrode::GetParticleSurfArea() {
    return part_surfArea;   // already computed once in the constructor
}

mfem::GridFunction& Electrode::GetParticleTotalLi() {
    part_totLi = 0.0;
    for (auto &particle : particles) {
        int p_id = particle->GetParticleID();
        part_totLi(p_id) = particle->GetParticleTotalLi();
    }
    return part_totLi;
}

double Electrode::GetTotalVolume() {
    return part_volume.Sum();
}

double Electrode::GetTotalSurfArea() {
    return part_surfArea.Sum();
}

double Electrode::GetTotalLi() {
    GetParticleTotalLi();   // refresh part_totLi first
    return part_totLi.Sum();
}

double Electrode::GetDoD() {
    GetParticleTotalLi();   // refresh part_totLi first
    return part_totLi.Sum()/part_volume.Sum();
}


double Electrode::GetElectrodeLength() {
    mfem::Mesh *mesh = fespace->GetMesh();

    mfem::Array<int> region2_marker(mesh->attributes.Max());
    region2_marker = 0;
    region2_marker[1] = 1;

    mfem::ConstantCoefficient one(1.0);
    mfem::LinearForm length_lf(fespace);
    length_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(one), region2_marker);
    length_lf.Assemble();

    return length_lf.Sum();
}


void Electrode::SaveAllConc(const std::string &prefix) {
    for (auto &particle : particles) {
        particle->SaveConc(prefix);
    }
}

void Electrode::SaveConcByID(int particle_id, const std::string &prefix) {
    for (auto &particle : particles) {
        if (particle->GetParticleID() == particle_id) {
            particle->SaveConc(prefix);
            return;
        }
    }
    std::cerr << "SaveConcByID: no particle with ID " << particle_id << " found" << std::endl;
}
