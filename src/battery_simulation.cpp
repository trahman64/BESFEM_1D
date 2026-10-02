#include "mfem.hpp"
#include "../includes/spherical_diffusion.hpp"
#include "../includes/linear_diffusion.hpp"
#include "../includes/stat_potential.hpp"
#include "../includes/atnv_calculator.hpp"
#include "../includes/diffCoeff_utils.hpp"
#include "../includes/potential_utils.hpp"
#include "../includes/electrode.hpp"
#include "../includes/butler_volmer.hpp"

const double F = 96485.332;
const double R = 8.314;
const double T = 300.0;

const double alpha_a = 0.5;
const double alpha_c = 0.5;

// const double i0 = 0.5e-3;   

const double aPv = 2.409e3;    
const double t_minus = 0.7619;

const double rho = 0.0312; 


// volume fraction of liquid 
const double eps_l_eld = 0.301;
const double eps_l_sep = 1.0;
// tortuousity of liquid
const double tau_l_eld = 1.521; 
const double tau_l_sep = 1.0;  
const double Ce0 = 0.001;
const double De = 0.25e-5;

const double dt = 1.0e-2; 

double kappa_s = 0.075;
// volume fraction of solid 
const double eps_s_eld = 0.699;
const double eps_s_sep = 1.0e0;
// tortuousity of liquid
const double tau_s_eld = 1.324; 
const double tau_s_sep = 1.0e-3;
double BvP = 2.99;
double BvE = -1.0;


double rad = 4.0e-4;
double Cp0 = 0.3;


int main(){
    mfem::Mesh mesh("../inputs/Mesh_80_1D02.mesh");
    int dim = mesh.Dimension();
    int order = 1;

    mfem::H1_FECollection fec(order, dim);
    mfem::FiniteElementSpace fespace(&mesh, &fec);
    
    // information of the electrode region
	mfem::Array<int> region2_dofs;
	mfem::Array<int> is_in_region2(fespace.GetVSize());
	is_in_region2 = 0;
	
	for (int e = 0; e < mesh.GetNE(); e++) {
		if (mesh.GetAttribute(e) == 2) {   // element belongs to region 2
			mfem::Array<int> vdofs;
			fespace.GetElementVDofs(e, vdofs);
			for (int i = 0; i < vdofs.Size(); i++) {
				is_in_region2[vdofs[i]] = 1;
			}
		}
	}
	
	// Convert the boolean mask into an actual list of DOF indices
	for (int i = 0; i < is_in_region2.Size(); i++) {
		if (is_in_region2[i]) {
			region2_dofs.Append(i);
		}
	}    
	int n_elde_nodes = region2_dofs.Size();
    

	mfem::GridFunction rxn(&fespace);;
	rxn = 0.0;

	for (int i = 0; i < n_elde_nodes; i++) {
		rxn(region2_dofs[i]) = 0.2e-6;
	}

// 	for (int i = 20; i <= 80; i++) {
// 		rxn(i) = 1e-6;
// 	}
// 	rxn.Print();



	// ======================================================	
	//   ______ _           _             _       _       
	//  |  ____| |         | |           | |     | |      
	//  | |__  | | ___  ___| |_ _ __ ___ | |_   _| |_ ___ 
	//  |  __| | |/ _ \/ __| __| '__/ _ \| | | | | __/ _ \
	//  | |____| |  __/ (__| |_| | | (_) | | |_| | ||  __/
	//  |______|_|\___|\___|\__|_|  \___/|_|\__, |\__\___|
	//                                       __/ |        
	//                                      |___/         
	// ======================================================	
	
	
	mfem::Array<int> nbc_bdr(mesh.bdr_attributes.Max());
	nbc_bdr = 0;
	nbc_bdr[0] = 1;   // mark whichever boundary attribute is the Neumann flux boundary	
	
    std::cout << "Creating linear diffusion" << std::endl;

	LinearDiffusion salt_electrolyte(&mesh, &fespace, nbc_bdr, Ce0, dt);
	mfem::GridFunction Ce(&fespace);
	Ce = salt_electrolyte.GetConcentration();	
	mfem::GridFunction De_gf(&fespace);
	De_gf = ComputeDamb(Ce);
	
	salt_electrolyte.SetWeightVector(eps_l_sep, eps_l_eld, tau_l_sep, tau_l_eld);
	salt_electrolyte.SetCoefficient(De_gf);
	salt_electrolyte.BuildOperator();

	
	
	// ======================================================		
	//   ______ _           _                 _      
	//  |  ____| |         | |               | |     
	//  | |__  | | ___  ___| |_ _ __ ___   __| | ___ 
	//  |  __| | |/ _ \/ __| __| '__/ _ \ / _` |/ _ \
	//  | |____| |  __/ (__| |_| | | (_) | (_| |  __/
	//  |______|_|\___|\___|\__|_|  \___/ \__,_|\___|
	// ======================================================	                                              
                                                                                       
                                         
    std::cout << "Creating particle diffusion" << std::endl;    // Reserving the size of the vector

	mfem::GridFunction parti_radii(&fespace);
	parti_radii = rad;
	
	Electrode NMC_electrode(&fespace, region2_dofs, parti_radii, 40, 1, Cp0, dt);
 	mfem::GridFunction Cp_surf(&fespace);
 	mfem::GridFunction Cp_mConc(&fespace);
//  	Cp_surf = 0.0;
//  	Cp_mConc = 0.0; 
 	
	Cp_surf = NMC_electrode.GetSurfaceConcentration();
	Cp_mConc = NMC_electrode.GetPartiMeanConcentration();
	
 	double Vol = NMC_electrode.GetTotalVolume();
 	double totX = NMC_electrode.GetTotalLi();
	double Xfr = NMC_electrode.GetDoD();
 	std::cout << Vol << " " << totX << " " << Xfr << std::endl;
 	

// 	SphericalDiffusion p_test(60, rad, 40, 1, Cp0, dt);


 
//     std::cout << "Creating potential" << std::endl;
// 
       
	// ============================================
	//   _____  _     _    _____ 
	//  |  __ \| |   (_)  / ____|
	//  | |__) | |__  _  | (___  
	//  |  ___/| '_ \| |  \___ \ 
	//  | |    | | | | |  ____) |
	//  |_|    |_| |_|_| |_____/ 
	// ============================================
                          

    mfem::Array<int> ess_bdr_s(mesh.bdr_attributes.Max());
    ess_bdr_s = 0;
    ess_bdr_s[1] = 1;
    
    mfem::GridFunction Kappa(&fespace);
    Kappa = ComputeKaps(Cp_mConc, is_in_region2);
//     Kappa = kappa_s;

	StatPotential solid_potential(&fespace, ess_bdr_s, BvP);
	solid_potential.SetWeightVector(eps_s_sep, eps_s_eld, tau_s_sep, tau_s_eld);   // sets region_weight
	solid_potential.SetCoefficient(Kappa);                                   // needs region_weight -- must come after SetWeightVector
	solid_potential.BuildOperator();                                         // needs weight_eff_kappa -- must come after SetCoefficient
	
	mfem::Vector AtnV_0(fespace.GetTrueVSize());
	AtnV_0 = 0.0e0;	


	// ============================================	
	//   _____  _     _   _      
	//  |  __ \| |   (_) | |     
	//  | |__) | |__  _  | |     
	//  |  ___/| '_ \| | | |     
	//  | |    | | | | | | |____ 
	//  |_|    |_| |_|_| |______|
	// ============================================
							  
	mfem::Array<int> ess_bdr_l(mesh.bdr_attributes.Max());   // -> if mesh is a pointer, else keep .
	ess_bdr_l = 0;
	ess_bdr_l[0] = 1;
	
	double Cst1 = F / R / T;
	double tc1 = (2 * t_minus - 1.0) / (2 * t_minus * (1.0 - t_minus));
	double tc2 = 1.0 / (2 * t_minus * (1.0 - t_minus)) * Cst1;
	double scaleConst = tc2 / tc1 * Cst1;
	
	mfem::GridFunction Dmp(&fespace);
	mfem::GridFunction Kpl(&fespace);
	
	Dmp = De_gf;
	Dmp *= tc1;
	
	Kpl = De_gf;
	Kpl *= scaleConst;
	Kpl *= Ce;
	
	StatPotential liquid_potential(&fespace, ess_bdr_l, BvE);
	liquid_potential.SetWeightVector(eps_l_sep, eps_l_eld, tau_l_sep, tau_l_eld);
	liquid_potential.SetCoefficient(Kpl);
	liquid_potential.BuildOperator();
	
	AtnVCalculator AtnVCalt(&fespace, Dmp);
	mfem::Vector &AtnV = AtnVCalt.Compute(Ce);   // reference, no unnecessary copy
		

	//  =============================================================  
	//   ____        _ _         __      __   _                     
	//  |  _ \      | | |        \ \    / /  | |                    
	//  | |_) |_   _| | |_ ___ _ _\ \  / /__ | |_ __ ___   ___ _ __ 
	//  |  _ <| | | | | __/ _ \ '__\ \/ / _ \| | '_ ` _ \ / _ \ '__|
	//  | |_) | |_| | | ||  __/ |   \  / (_) | | | | | | |  __/ |   
	//  |____/ \__,_|_|\__\___|_|    \/ \___/|_|_| |_| |_|\___|_|   
	//  =============================================================                                                            

	mfem::GridFunction OCV_gf(&fespace);
	mfem::GridFunction i0_gf(&fespace);	
	mfem::GridFunction phs_gf(&fespace);
	mfem::GridFunction phl_gf(&fespace);
	OCV_gf = 0.0;
	i0_gf = 0.0;
	phs_gf = 0.0;
	phl_gf = 0.0;
 
	i0_gf = Compute_i0(Cp_surf, is_in_region2);
	OCV_gf = Compute_pOCV(Cp_surf, is_in_region2);
	phs_gf = solid_potential.GetPotential();
	phl_gf = liquid_potential.GetPotential();	 
	
// 	i0_gf.Print();
	
	ButlerVolmer cellRxn(&fespace, region2_dofs, alpha_a, alpha_c, Cst1, F);                                              
	rxn = cellRxn.Compute(Ce, Cp_surf, i0_gf, OCV_gf, phs_gf, phl_gf);
	
	mfem::GridFunction source_ely(&fespace);
	source_ely = rxn;
	source_ely *= (aPv*t_minus);
	source_ely.Neg();

	
 	mfem::GridFunction source_eld(&fespace);
 	source_eld = rxn;
	source_eld /= rho;		
	source_eld.Print();

	mfem::GridFunction source_phs(&fespace);
	source_phs = rxn;
	source_phs *= (aPv*F);
	
	mfem::GridFunction source_phl(&fespace);
	source_phl = rxn;
	source_phl *= aPv;
	source_phl.Neg();
	
	double totCrnt = cellRxn.GetTotalRxnCurrent();
	std::cout << totCrnt << "---" << std::endl;
	

    std::cout << "Starting loop" << std::endl;

    int num_steps = 10000;
//     double rxn = 0.02e-6;
    double frx_p = 0.0;

    for (int iter = 0; iter <num_steps; iter++){
    	if (iter % 20 == 0 ){
			std::cout << "Step " << iter << ": diffusion" << std::endl;
    	}


// 		p_test.UpdateOperator();
// 		p_test.Stepping(source_eld(60));
// 		std::cout << source_eld(60) << "xxx" << std::endl;


// 		Ce = salt_electrolyte.GetConcentration();	
// 		De_gf = ComputeDamb(Ce);
// 		salt_electrolyte.UpdateOperator();
// 		salt_electrolyte.Stepping(source_ely);

		NMC_electrode.UpdateOperators();		
		NMC_electrode.Stepping(source_eld);
		

// 		
// 		Cp_surf = NMC_electrode.GetSurfaceConcentration();

// 		std::cout << p_test.GetConcentrationAt(rad) << " --->> " << 
// 			Cp_surf(60) << std::endl;		
		
// 		Cp_mConc = NMC_electrode.GetPartiMeanConcentration();
				
		// Recompute Kpl's values from the NEW De_gf/Ce, same object:
// 		Kpl = De_gf;
// 		Kpl *= scaleConst;
// 		Kpl *= Ce;
// 		liquid_potential.UpdateOperator();   // now correctly reflects the new Kpl
	
		// Similarly, Dmp needs recomputing if it should track the new De_gf too:
// 		Dmp = De_gf;
// 		Dmp *= tc1;
// 		AtnVCalt.UpdateDmp();
// 	
// 		mfem::Vector &AtnV = AtnVCalt.Compute(Ce);
				
// 		Kappa = ComputeKaps(Cp_mConc, is_in_region2);			
// 		solid_potential.UpdateOperator();
// 		liquid_potential.UpdateOperator();		


// 		i0_gf = Compute_i0(Cp_surf, is_in_region2);
// 		OCV_gf = Compute_pOCV(Cp_surf, is_in_region2);
// 		phs_gf = solid_potential.GetPotential();
// 		phl_gf = liquid_potential.GetPotential();
// 		
// 		rxn = cellRxn.Compute(i0_gf, OCV_gf, phs_gf, phl_gf);	
					

// 		solid_potential.Solve(source_phs, AtnV_0, BvP);
// 		liquid_potential.Solve(source_phl, AtnV, BvE);	


			
				     
//             << std::endl;
    }

// 	NMC_electrode.SaveAllConc();
	NMC_electrode.SaveConcByID(60);
// 	std::cout << NMC_electrode.GetDoD() << std::endl;
	
	

//     p1.SaveConc();
// 	mfem::GridFunction C(&fespace);
//     C = salt_electrolyte.GetConcentration(); 
//     C.Print();
//     salt_electrolyte.SaveConc("elyConc.gf");
//     double MnConc = salt_electrolyte.GetMeanConcentration();
//     std::cout << "MC = " << MnConc << std::endl;
 
//  	std::cout << particles[3].GetMeanConcentration() 
//  		<< "  " << particles[57].GetMeanConcentration() << std::endl;
//  	particles[3].SaveConc();
//  	particles[57].SaveConc();	 
//  	p_test.SaveConc();
 	
// 	mfem::GridFunction Phi(&fespace);    
//     Phi = solid_potential.GetPotential();
    solid_potential.SavePote("solid_phi.gf");
//     Phi = solid_potential.GetPotential();
    liquid_potential.SavePote("liquid_phi.gf");    
    cellRxn.Save();
//     rxn.Print();
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