#include "mfem.hpp"
#include "../includes/diffusion.hpp"
#include "../includes/potential.hpp"
#include "../includes/radial_diffusion.hpp"

const double F = 96485.332;
const double R = 8.314;
const double T = 300.0;

const double alpha = 0.5;

const double i0 = 0.5e-3;   

const double a = 2.409e3;    
const double t_minus = 0.7619;

const double rho = 0.0312;  

double ocv(double x){
    return 1.095 * x * x - 8.234e-7 * std::exp(14.32 * x) + 4.692 * std::exp(-0.5389 * x);
}

double butlerVolmer(double phi_s, double phi_e, double x){
    double phi_ocv = ocv(x);
    double n = (phi_s - phi_e) - phi_ocv;
    double rxn = (i0/F) * (std::exp((-alpha*F/(R*T))*n) - std::exp(((1-alpha)*F/(R*T))*n));
    return rxn;
}



int main(){
    mfem::Mesh mesh("../inputs/Mesh_80_1D01.mesh");
    int dim = mesh.Dimension();
    int order = 1;

    mfem::H1_FECollection fec(order, dim);
    mfem::FiniteElementSpace fespace(&mesh, &fec);

    std::cout << "Creating diffusion" << std::endl;
    Diffusion salt_electrolyte(&mesh, &fespace);

    std::cout << "Creating radial" << std::endl;
    Radial_Diffusion radial_diffusion(&mesh, &fespace);

    std::cout << "Creating potential" << std::endl;
    Potential poission(&mesh, &fespace,&salt_electrolyte);
    
    
    
    std::cout << "Starting loop" << std::endl;

    double dt = 1e-4;
    int num_steps = 10000;
    double rxn = 0.02e-6;

    for (int i = 0; i <num_steps; i++){
        std::cout << "Step " << i << ": diffusion" << std::endl;
        salt_electrolyte.Stepping(dt,rxn);
        
//         std::cout << "Step " << i << ": radial_diffusion" << std::endl;
//         radial_diffusion.Stepping(dt,rxn);
// 
//         std::cout << "11. Radial finished" << std::endl;
//         std::cout << "Step " << i << ": potential" << std::endl;
//         poission.Solve(rxn);
// 
//         double x_surface = radial_diffusion.GetSurfaceConc();
//         const mfem::GridFunction &phi_s = poission.GetSolidPotential();
//         const mfem::GridFunction &phi_e = poission.GetLiquidPotential();
// 
//         double phi_s_value = phi_s[phi_s.Size() - 1];
//         double phi_e_value = phi_e[phi_e.Size() - 1];
// 
//         rxn = butlerVolmer(phi_s_value, phi_e_value, x_surface);
// 
//         std::cout
//             << "Step: " << i
//             << "  X_surface: " << x_surface
//             << "  phi_s: " << phi_s_value
//             << "  phi_e: " << phi_e_value
//             << "  rxn: " << rxn
//             << std::endl;
    }
    salt_electrolyte.Save();
    
    
    mfem::GridFunction C(&fespace);
    C = salt_electrolyte.GetConcentration();
    std::cout << (C[1]-C[0])/1e-5 << " " << C[0] << " " << C[1] << std::endl;
    
    mfem::Vector epsilon_vector(mesh->attributes.Max());
    epsilon_vector(0) = 1.0;
    
	mfem::ConstantCoefficient one(1.0);
	mfem::LinearForm mass_lf(&fespace);
	mass_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
	mass_lf.Assemble();
	
	double integral_u = mass_lf(C);      // ∫ u dx  (LinearForm::operator() computes the dot product)
	double volume      = mass_lf.Sum();   // ∫ 1 dx = domain volume
	
	double mean = integral_u / volume; 
   std::cout << mean << std::endl;
//     radial_diffusion.Save();
//     poission.Save();
    

    return 0;
}