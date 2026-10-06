#include "mfem.hpp"
#include <cmath>


mfem::GridFunction Compute_Kaps(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2)
{
    mfem::GridFunction Kps(Cn.FESpace());
	Kps = 0.0075;
    return Kps;
}


mfem::GridFunction Compute_i0(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2)
{	
	double pwV;
    mfem::GridFunction i0_gf(Cn.FESpace());
    
    for (int i = 0; i < is_in_region2.Size(); i++) {
    	if (is_in_region2[i] == 1){   	
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