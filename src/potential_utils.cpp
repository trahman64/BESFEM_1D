#include "mfem.hpp"
#include <cmath>


mfem::GridFunction ComputeKaps(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2)
{
    mfem::GridFunction Kps(Cn.FESpace());

	double C_ref = 0.5;
	double Kps_ref = 0.01929 + 0.7045 * std::tanh(2.399 * C_ref) - 
		0.7238 * std::tanh(2.412 * C_ref);    
		
    for (int i = 0; i < is_in_region2.Size(); i++) {
    	if (is_in_region2[i] == 1){
			Kps(i) = 0.01929 + 0.7045 * std::tanh(2.399 * Cn(i)) - 
				0.7238 * std::tanh(2.412 * Cn(i));
		} else {
			Kps(i) = Kps_ref;
    	}

    }
    return Kps;
}


mfem::GridFunction Compute_i0(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2)
{	
	double C_ref = 0.5;
	double pwV = -0.2 * (C_ref-0.37) - 0.9376 * std::tanh(8.961*C_ref- 3.195) - 1.559;
	double i0_ref = std::exp(pwV);
	
    mfem::GridFunction i0_gf(Cn.FESpace());
    
    for (int i = 0; i < is_in_region2.Size(); i++) {
    	if (is_in_region2[i] == 1){
			pwV = -0.2 * (Cn(i)-0.37) - 0.9376 * std::tanh(8.961*Cn(i)- 3.195) - 1.559;
			i0_gf(i) = std::exp(pwV);    	
    	} else {
    		i0_gf(i) = i0_ref;
    	}
    } 
    return i0_gf;
}


mfem::GridFunction Compute_pOCV(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2)
{	
	double C_ref = 0.5;
	double OCV_ref = 1.095 * C_ref*C_ref - 8.234e-7*std::exp(14.32*C_ref) + 
		4.692*std::exp(-0.5389*C_ref);
	
    mfem::GridFunction OCV_gf(Cn.FESpace());
    
    for (int i = 0; i < is_in_region2.Size(); i++) {
    	if (is_in_region2[i] == 1){
			OCV_gf(i) = 1.095 * Cn(i)*Cn(i) - 8.234e-7*std::exp(14.32*Cn(i)) + 
				4.692*std::exp(-0.5389*C_ref);    	
    	} else {
    		OCV_gf(i) = OCV_ref;
    	}
    } 
    return OCV_gf;
}