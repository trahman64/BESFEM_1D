#include "mfem.hpp"
// #include "../includes/diffusion.hpp"
// #include "../includes/potential.hpp"
// #include "../includes/radial_diffusion.hpp"
#include "../includes/spherical_diffusion.hpp"
#include "../includes/linear_diffusion.hpp"
#include "../includes/stat_potential.hpp"

const double F = 96485.332;
const double R = 8.314;
const double T = 300.0;

const double alpha = 0.5;

const double i0 = 0.5e-3;   

const double a = 2.409e3;    
const double t_minus = 0.7619;

const double rho = 0.0312; 


// volume fraction of liquid 
const double eps_l_eld = 0.301;
const double eps_l_sep = 1.0;
// tortuousity of liquid
const double tau_l_eld = 1.521; 
const double tau_l_sep = 1.0;  
const double C0 = 0.001;
const double De = 0.25e-5;

const double dt = 1.0e-4; 

double kappa_s = 0.075;
// volume fraction of solid 
const double eps_s_eld = 0.699;
const double eps_s_sep = 1.0e0;
// tortuousity of liquid
const double tau_s_eld = 1.324; 
const double tau_s_sep = 1.0e-3;
double BvP = 0.0;



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
    mfem::Mesh mesh("../inputs/Mesh_80_1D02.mesh");
    int dim = mesh.Dimension();
    int order = 1;

    mfem::H1_FECollection fec(order, dim);
    mfem::FiniteElementSpace fespace(&mesh, &fec);

	mfem::GridFunction rxn(&fespace);;
	rxn = 0.0;
	for (int i = 20; i <= 80; i++) {
		rxn(i) = a*0.02e-6;
	}
	
	mfem::GridFunction source_ely(&fespace);
	source_ely = rxn;
	source_ely *= -1.0;
	source_ely *= t_minus;
	
    std::cout << "Creating linear diffusion" << std::endl;
	LinearDiffusion salt_electrolyte(&mesh, &fespace, eps_l_sep, eps_l_eld,
		De, tau_l_sep, tau_l_eld, C0, dt);

//     std::cout << "Creating radial" << std::endl;
//     Radial_Diffusion radial_diffusion(&mesh, &fespace);
//     SphericalDiffusion particle_1(4.0e-4, 40, 1.5e-10, 60, 1e-2, 1, 0.3);    

//     std::cout << "Creating potential" << std::endl;
//     Potential poission(&mesh, &fespace,&salt_electrolyte);

    // Region-wise porosity/tortuosity weighting, same pattern as LinearDiffusion
    mfem::Vector weight_vector(mesh.attributes.Max());
//     weight_vector(0) = eps_s_sep / (tau_s_sep * tau_s_sep);   // separator
//     weight_vector(1) = eps_s_eld / (tau_s_eld * tau_s_eld);   // electrode
//     mfem::PWConstCoefficient region_weight(weight_vector);
//   

	mfem::Vector wv_s = ComputeWeightVector(&fespace, eps_s_sep, eps_s_eld, tau_s_sep, tau_s_eld);
	mfem::PWConstCoefficient region_weight_s(wv_s);
	
	mfem::GridFunction Kappa_s(&fespace);
	Kappa_s = kappa_s;
	mfem::GridFunctionCoefficient kappa_coeff_s(&Kappa_s);
	mfem::ProductCoefficient effective_kappa_s(region_weight_s, kappa_coeff_s);
  
//     mfem::GridFunction Kappa(&fespace);
//     Kappa = kappa_s;
//     mfem::GridFunctionCoefficient kappa_coeff(&Kappa);          // Kappa itself, spatially varying
//     mfem::ProductCoefficient effective_kappa(region_weight, kappa_coeff);
    mfem::ProductCoefficient wt_eff_kappa_s(region_weight_s, effective_kappa_s);    

    mfem::Array<int> ess_bdr_s(mesh.bdr_attributes.Max());
    ess_bdr_s = 0;
    ess_bdr_s[1] = 1;

    StatPotential sold_potent(&fespace, wt_eff_kappa_s, ess_bdr_s);

	mfem::Vector Additional_0(fespace.GetTrueVSize());
	Additional_0 = 0.0e0;
	

    
	mfem::GridFunction source_phs(&fespace);
	source_phs = rxn;
	source_phs *= F;
	
   
	mfem::GridFunction C(&fespace);
	C = salt_electrolyte.GetConcentration();  

// 	Vector
	    
    std::cout << "Starting loop" << std::endl;

    double dt = 1e-2;
    int num_steps = 10000;
//     double rxn = 0.02e-6;
    double frx_p = 0.0;

    for (int i = 0; i <num_steps; i++){
        std::cout << "Step " << i << ": diffusion" << std::endl;
		salt_electrolyte.Stepping(source_ely);
		
		sold_potent.Solve(source_phs, Additional, BvP);
				     
//         salt_electrolyte.Stepping(dt,rxn);
               
//         std::cout << "Step " << i << ": radial_diffusion" << std::endl;
//         particle_1.Stepping(rxn);
        
//         frx_p = particle_1.GetMeanConcentration();
        
//         std::cout << i << "--" << frx_p << std::endl;
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
    C = salt_electrolyte.GetConcentration(); 
    C.Print();
    salt_electrolyte.SaveConc();
    double MnConc = salt_electrolyte.GetMeanConcentration();
    std::cout << "MC = " << MnConc << std::endl;
 
 
	mfem::GridFunction Php(&fespace);    
    Php = sold_potent.GetPotential();
    sold_potent.SavePote();
    
// 	particle_1.SaveConc();
    

// 	C = salt_electrolyte.GetConcentration();  
// 	C.Print();  

// 	integral_u = mass_lf(C);      // ∫ u dx  (LinearForm::operator() computes the dot product)
// 	volume     = mass_lf.Sum();   // ∫ 1 dx = domain volume
// 	
// 	mean = integral_u / volume; 
// 	std::cout << mean << std::endl;
// //     radial_diffusion.Save();
//     poission.Save();
//     
    

    return 0;
}