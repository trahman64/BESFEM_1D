#include "mfem.hpp"
#include "../includes/spherical_diffusion.hpp"
#include "../includes/linear_diffusion.hpp"
#include "../includes/stat_potential.hpp"
#include "../includes/atnv_calculator.hpp"
#include "../includes/diffCoeff_utils.hpp"
#include "../includes/potential_utils.hpp"
#include "../includes/electrode.hpp"
#include "../includes/butler_volmer.hpp"
#include "../includes/cell_kinetics.hpp"
#include <fstream>
#include <iomanip>



const double F = 96485.332;
const double R = 8.314;
const double T = 300.0;

const double alpha_a = 0.5;
const double alpha_c = 0.5;

// volume fraction of liquid 
const double eps_l_eld = 0.301;
const double eps_l_sep = 1.0;
// tortuousity of liquid
const double tau_l_eld = 1.521; 
const double tau_l_sep = 1.0;  
const double Ce0 = 0.001;
// const double De0 = 0.25e-5;
const double t_minus = 0.7619;


// volume fraction of solid 
const double eps_s_eld = 0.699;
const double eps_s_sep = 1.0e-3;
// tortuousity of solid
const double tau_s_eld = 1.324; 
const double tau_s_sep = 1.0e-3;

const double aPv = 2.409e3;    
const double rho = 0.0312; 

double rad = 4.0e-4;
double Cp0 = 0.3;
double X_e = 0.3;
double X_f = 0.95;

// double BvP =  3.0; 
// double BvE = -1.081745;
double BvP =  3.081745; 
double BvE = -1.0;

double dCV = 0.0;
double cut_off = 3.09;
double C_rate = 0.5;
double CV_sr = 1e-3;
double tols = 1e-12;
double toll = 1e-12;
int internal_maxiter = 200;

const double dt = 1.0e-2; 
double tm = 0.0;
int num_steps = (3600.0*2/dt); // 10; //  





int main(){

	// ==================================
	//                       _     
	//                      | |    
	//   _ __ ___   ___  ___| |__  
	//  | '_ ` _ \ / _ \/ __| '_ \ 
	//  | | | | | |  __/\__ \ | | |
	//  |_| |_| |_|\___||___/_| |_|
	// ==================================
                            
                                   
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

// 	for (int i = 0; i < n_elde_nodes; i++) {
// 		rxn(region2_dofs[i]) = 0.2e-6;
// 	}


	std::ofstream outfile("DoD_output.csv");
	outfile << "time,DoD,totCrnt,CellVoltage\n";

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
	De_gf = Compute_Damb(Ce);

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

	mfem::GridFunction parti_initC(&fespace);
	parti_initC = Cp0;
	
	mfem::GridFunction parti_radii(&fespace);
	parti_radii = rad;
	
	Electrode NMC_electrode(&fespace, region2_dofs, 
		parti_radii, 40, 1, parti_initC, dt);
 	mfem::GridFunction Cp_surf(&fespace);
 	mfem::GridFunction Cp_mConc(&fespace);
 	Cp_surf = 0.0;
 	Cp_mConc = 0.0; 
 	
	Cp_surf = NMC_electrode.GetSurfaceConcentration();
	Cp_mConc = NMC_electrode.GetPartiMeanConcentration();
	
 	double Vol = NMC_electrode.GetTotalVolume();
 	double totX = NMC_electrode.GetTotalLi();
	double Xfr = NMC_electrode.GetDoD();
 	std::cout << "initial:" << Vol << " " << totX << " " << Xfr << std::endl;

	double next_output_DoD = Cp0;   // the next threshold to trigger on
	const double DoD_step = 0.005; 	

// 	SphericalDiffusion p_test(60, rad, 40, 1, Cp0, dt);

 
    std::cout << "Creating potentials" << std::endl; 
       
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
	Kappa = Compute_Kaps(Cp_mConc, is_in_region2);    

	StatPotential solid_potential(&fespace, ess_bdr_s, BvP);
	solid_potential.SetWeightVector(eps_s_sep, eps_s_eld, tau_s_sep, tau_s_eld);   // sets region_weight
	solid_potential.SetCoefficient(Kappa);                                   // needs region_weight -- must come after SetWeightVector
	solid_potential.BuildOperator();                                         // needs weight_eff_kappa -- must come after SetCoefficient
	
	mfem::Vector AtnV_0(fespace.GetTrueVSize());
	AtnV_0 = 0.0e0;	

	mfem::GridFunction phs_old_gf(&fespace);
	mfem::GridFunction phs_tmp_gf(&fespace);	

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
	double tc1, tc2, scaleConst;
	ComputeConstantce(t_minus, Cst1, tc1, tc2, scaleConst);	
	
	mfem::GridFunction Dmp(&fespace);
	mfem::GridFunction Kpl(&fespace);
	
	Dmp = Compute_Dmp(De_gf, tc1);
	Kpl = Compute_Kpl(De_gf, Ce, scaleConst);
	
	StatPotential liquid_potential(&fespace, ess_bdr_l, BvE);
	liquid_potential.SetWeightVector(eps_l_sep, eps_l_eld, tau_l_sep, tau_l_eld);
	liquid_potential.SetCoefficient(Kpl);
	liquid_potential.BuildOperator();
	
	AtnVCalculator AtnVCalt(&fespace, Dmp);
	mfem::Vector &AtnV = AtnVCalt.Compute(Ce);   // reference, no unnecessary copy
		
	mfem::GridFunction phl_old_gf(&fespace);
	mfem::GridFunction phl_tmp_gf(&fespace);	
	
    std::cout << "Creating reaction" << std::endl;
    
	//  =============================================================  
	//   ____        _ _         __      __   _                     
	//  |  _ \      | | |        \ \    / /  | |                    
	//  | |_) |_   _| | |_ ___ _ _\ \  / /__ | |_ __ ___   ___ _ __ 
	//  |  _ <| | | | | __/ _ \ '__\ \/ / _ \| | '_ ` _ \ / _ \ '__|
	//  | |_) | |_| | | ||  __/ |   \  / (_) | | | | | | |  __/ |   
	//  |____/ \__,_|_|\__\___|_|    \/ \___/|_|_| |_| |_|\___|_|   
	//  =============================================================                                                            

		
	CellKinetics EC_Rxn(&fespace, region2_dofs, is_in_region2, 
		solid_potential, liquid_potential, aPv, alpha_a, alpha_c,
		Cst1, F);
	rxn = EC_Rxn.Compute(Ce, Cp_surf);		

	mfem::GridFunction source_phs(&fespace);
	mfem::GridFunction source_phl(&fespace);

	
	double dphs_2n = 0.0;
	double dphl_2n = 0.0;


	mfem::GridFunction source_ely(&fespace);
 	mfem::GridFunction source_eld(&fespace);
	
		
	double eld_length = NMC_electrode.GetElectrodeLength();
	double cap_gl = rho * (X_f-X_e) * eld_length * eps_s_eld;
	double pVA_scale =  (NMC_electrode.GetTotalVolume() / 
		NMC_electrode.GetTotalSurfArea())* aPv; // / eps_s_eld;
		
	double I_trgt =  cap_gl/(3600.0/C_rate) * pVA_scale;
	
// 	double totCrnt = cellRxn.GetTotalRxnCurrent();
	double totCrnt = EC_Rxn.GetTotalRxnCurrent();
	std::cout << BvP << " : " << totCrnt << " ---> " << I_trgt << std::endl;
// 	std::cout << NMC_electrode.GetTotalVolume() / 
// 		NMC_electrode.GetTotalSurfArea()  << std::endl;
	

    std::cout << "Starting loop" << std::endl;

  
// 	num_steps = 28000;
    for (int t_step = 0; t_step <= num_steps; t_step++){
    	if (t_step % 1000 == 0 ){
			std::cout << "Step " << t_step << ": diffusion" << "  " 
			<< NMC_electrode.GetDoD() << "  " << tm << " || " <<
			(NMC_electrode.GetDoD()-0.3)/(X_f-X_e) << "  " 
			<< tm/(3600.0/C_rate) << std::endl;
    	}


// 		p_test.UpdateOperator();
// 		p_test.Stepping(source_eld(60));
// 		std::cout << source_eld(60) << "xxx" << std::endl;

		source_ely = rxn;
		source_ely *= t_minus;
		source_ely.Neg();
		Ce = salt_electrolyte.GetConcentration();	
		De_gf = Compute_Damb(Ce);
		salt_electrolyte.UpdateOperator();
		salt_electrolyte.Stepping(source_ely);

		source_eld = rxn;
		source_eld /= (aPv * rho * eps_s_eld);	
		NMC_electrode.UpdateOperators();		
		NMC_electrode.Stepping(source_eld);
		Cp_surf = NMC_electrode.GetSurfaceConcentration();
		Cp_mConc = NMC_electrode.GetPartiMeanConcentration();

// 		std::cout << p_test.GetConcentrationAt(rad) << " --->> " << 
// 			Cp_surf(60) << std::endl;		
		Kappa = Compute_Kaps(Cp_mConc, is_in_region2);			
		solid_potential.UpdateOperator();		

				
		// Recompute Kpl's values from the NEW De_gf/Ce, same object:
		Kpl = Compute_Kpl(De_gf, Ce, scaleConst);
		liquid_potential.UpdateOperator();   // now correctly reflects the new Kpl
	
		// Similarly, Dmp needs recomputing if it should track the new De_gf too:
		Dmp = Compute_Dmp(De_gf, tc1);		
		AtnVCalt.UpdateDmp();
	
		mfem::Vector &AtnV = AtnVCalt.Compute(Ce);
		liquid_potential.UpdateOperator();	
				

		// internal loop for rxn, phs, and phl
		dphs_2n = 1.0;
		dphl_2n = 1.0;
		for (int internal = 0; internal < internal_maxiter; internal++) {
			phs_old_gf = EC_Rxn.GetPhs();
			phl_old_gf = EC_Rxn.GetPhl();
		
			rxn = EC_Rxn.Compute(Ce, Cp_surf);
		
			source_phs = rxn; source_phs *= F;
			source_phl = rxn; source_phl.Neg();
		
			solid_potential.Solve(source_phs, AtnV_0, BvP);
			liquid_potential.Solve(source_phl, AtnV, BvE);   // only ONE call now
		
			phs_tmp_gf = solid_potential.GetPotential();
			phs_tmp_gf -= phs_old_gf;
			dphs_2n = phs_tmp_gf.Norml2();
		
			phl_tmp_gf = liquid_potential.GetPotential();   // no second Solve() before this
			phl_tmp_gf -= phl_old_gf;
			dphl_2n = phl_tmp_gf.Norml2();
		
			if (dphs_2n < tols && dphl_2n < toll) {
				break;
			}
		}
		

		totCrnt = EC_Rxn.GetTotalRxnCurrent();
		dCV = std::copysign(CV_sr, I_trgt - totCrnt);
		dCV *= dt;
		BvE += dCV;

		tm += dt; 
		
		if (NMC_electrode.GetDoD() >= next_output_DoD) {
			// Save/output data here
			std::cout << "DoD reached " << next_output_DoD << 
				" at step " << t_step << std::endl;
	
			NMC_electrode.SaveAllConc("output/NMC_DoD_" + 
				std::to_string(next_output_DoD) + "gf");
			// ... any other output you want (voltage, current, etc.) ...

			salt_electrolyte.SaveConc("output/elyConc_" + 
				std::to_string(next_output_DoD) + ".gf");
			
			liquid_potential.SavePote("output/liquid_phi_" + 
				std::to_string(next_output_DoD) + ".gf");  	
			
			solid_potential.SavePote("output/solid_phi_" + 
				std::to_string(next_output_DoD) + ".gf");
			
		    outfile << std::setprecision(10)
					<< tm << ","
					<< NMC_electrode.GetDoD() << ","
					<< totCrnt << ","
					<< BvP - BvE << std::endl;;
	
			next_output_DoD += DoD_step;   // advance to the NEXT threshold
		}

		if (BvP-BvE < cut_off) {break;}
		
		
		if (t_step % 1000 == 0 ){
			std::cout << NMC_electrode.GetDoD() <<  " ^^^ " << 
				BvP - BvE << " __ " << totCrnt << " --> " << 
				I_trgt << std::endl; 
		}


    }
    outfile.close();

// 	NMC_electrode.SaveAllConc();
// 	NMC_electrode.SaveConcByID(60);
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
//     solid_potential.SavePote("solid_phi.gf");
//     Phi = solid_potential.GetPotential();
//     liquid_potential.SavePote("liquid_phi.gf");    
//     cellRxn.Save();
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