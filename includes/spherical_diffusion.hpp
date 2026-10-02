#ifndef SPHERICAL_DIFFUSION_HPP
#define SPHERICAL_DIFFUSION_HPP

#include "mfem.hpp"
#include <string>

class SphericalDiffusion {
public:
    SphericalDiffusion(int x_idx,
                        double radius,
                        int n_elements,
                        int fe_order,
                        double C0,
                        double dt);

    ~SphericalDiffusion();

    void UpdateOperator();   // call after Stepping() updates C, to refresh D_li-dependent K/Tmat
    void Stepping(double surface_flux);

    double GetMeanConcentration();
	double GetParticleVolume();
	double GetParticleSurfaceArea();	
	double GetParticleTotalLi();    
    double GetConcentrationAt(double x);
    int GetParticleID() const;



    mfem::GridFunction& GetConcentration();
    void SaveConc(const std::string &prefix = "sphere");
    void SaveMesh(const std::string &prefix = "sphere");

private:
    double R;
    int order;
    int x_idx;
    double dt;

    mfem::Mesh *mesh;
    mfem::FiniteElementCollection *fec;
    mfem::FiniteElementSpace *fespace;

    mfem::GridFunction C;
    mfem::Vector C_prev;

    mfem::GridFunction D_li;
    mfem::GridFunction r;                                 // now a member -- referenced by r2_coeff
    mfem::GridFunctionCoefficient *r2_coeff = nullptr;     // now a member
    mfem::GridFunctionCoefficient *D_li_coeff = nullptr;   // now a member (was "DC")
    mfem::ProductCoefficient *D_r2 = nullptr;              // now a member

    mfem::Array<int> outer_bdr_marker;

    mfem::BilinearForm *M = nullptr;
    mfem::BilinearForm *K = nullptr;
    mfem::SparseMatrix *M_mat, *K_mat;
    mfem::SparseMatrix *Tmat = nullptr;

    mfem::GSSmoother *prec = nullptr;
    mfem::CGSolver solver;

	mfem::LinearForm vol_lf;     
	
    void BuildMesh(int n_elements);
    void BuildOperators();
   
	double partiVolume;
    double partiSurfArea;
    


};

#endif