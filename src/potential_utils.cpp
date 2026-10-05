#include "mfem.hpp"
#include <cmath>


mfem::GridFunction Compute_Kaps(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2)
{
    mfem::GridFunction Kps(Cn.FESpace());

	double C_ref = 0.5;
	double Kps_ref = 0.01929 + 0.7045 * std::tanh(2.399 * C_ref) - 
		0.7238 * std::tanh(2.412 * C_ref);    
		
//     for (int i = 0; i < is_in_region2.Size(); i++) {
//     	if (is_in_region2[i] == 1){
// 			Kps(i) = 0.01929 + 0.7045 * std::tanh(2.399 * Cn(i)) - 
// 				0.7238 * std::tanh(2.412 * Cn(i));
// 		} else {
// 			Kps(i) = Kps_ref;
//     	}
//     }
	Kps = Kps_ref;
    return Kps;
}


mfem::GridFunction Compute_i0(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2)
{	
	double pwV;
    mfem::GridFunction i0_gf(Cn.FESpace());
    
    for (int i = 0; i < is_in_region2.Size(); i++) {
    	if (is_in_region2[i] == 1){
// 			pwV = -0.2 * (Cn(i)-0.37) - 0.9376 * std::tanh(8.961*Cn(i)- 3.195) - 1.559;
// 			i0_gf(i) = pow(10.0, pwV) * 1.0e-3;    	
			i0_gf(i) = 0.2*1e-3;
    	} else {
    		i0_gf(i) = 0.0;
    	}
    } 
    return i0_gf;
}


mfem::GridFunction Compute_OCV(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2)
{		
    mfem::GridFunction OCV_gf(Cn.FESpace());
    
    for (int i = 0; i < is_in_region2.Size(); i++) {
    	if (is_in_region2[i] == 1){
			OCV_gf(i) = 1.095 * Cn(i)*Cn(i) - 8.234e-7*std::exp(14.32*Cn(i)) + 
				4.692*std::exp(-0.5389*Cn(i));    	
    	} else {
    		OCV_gf(i) = 0.0;
    	}
    } 
    return OCV_gf;
}

mfem::GridFunction Compute_Kpl(mfem::GridFunction &De_gf, 
							   mfem::GridFunction &Cn_gf, 
							   double scaleConst)
{
	mfem::GridFunction Kpl_gf(De_gf.FESpace());
	Kpl_gf = De_gf;
	Kpl_gf *= scaleConst;	
	Kpl_gf *= Cn_gf;
	
	return Kpl_gf;
}