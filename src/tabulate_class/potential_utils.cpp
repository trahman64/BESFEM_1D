#include "mfem.hpp"
#include <cmath>


mfem::GridFunction Compute_Kaps(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2,
							    TabulatedData &KLi_table)
{
    mfem::GridFunction Kps(Cn.FESpace());
	Kps = KLi_table.Interp(&Cn);
	double Kps_ref = 0.0075;  

    for (int i = 0; i < is_in_region2.Size(); i++) {
    	if (is_in_region2[i] == 0){
			Kps(i) = Kps_ref;
		}
    }
    return Kps;
}


mfem::GridFunction Compute_i0(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2,
 							  TabulatedData &i0_table)
{	
    mfem::GridFunction i0_gf(Cn.FESpace());
	i0_gf = i0_table.Interp(&Cn);
    
    for (int i = 0; i < is_in_region2.Size(); i++) {
    	if (is_in_region2[i] == 0){   	
    		i0_gf(i) = 0.0;
    	}
    } 
    return i0_gf;
}


mfem::GridFunction Compute_OCV(mfem::GridFunction &Cn, mfem::Array<int> is_in_region2,
							   TabulatedData &OCV_table)
{		
    mfem::GridFunction OCV_gf(Cn.FESpace());
	OCV_gf = iOCV_table.Interp(&Cn);    
    
    for (int i = 0; i < is_in_region2.Size(); i++) {
    	if (is_in_region2[i] == 0){
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