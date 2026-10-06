#include "mfem.hpp"
#include <cmath>

// in a small header, e.g. potential_utils.hpp, or just above SolidPotential's definition
mfem::GridFunction Compute_Damb(mfem::GridFunction &Cn)
{
	mfem::GridFunction De(Cn.FESpace());
	De = 2.5e-6;
 	return De;
}

void ComputeConstantce(double t_minus, double Cst1,
                        double &tc1, double &tc2, double &scaleConst)
{
    tc1 = (2 * t_minus - 1.0) / (2 * t_minus * (1.0 - t_minus));
    tc2 = 1.0 / (2 * t_minus * (1.0 - t_minus)) * Cst1;
    scaleConst = tc2 / tc1 * Cst1;
}


mfem::GridFunction Compute_Dmp(mfem::GridFunction &De_gf, double tc1)
{
	mfem::GridFunction Dmp_gf(De_gf.FESpace());
	Dmp_gf = De_gf;
	Dmp_gf *= tc1;	
	
	return Dmp_gf;
}

mfem::GridFunction Compute_DLi(mfem::GridFunction &Cn)
{
	mfem::GridFunction DLi(Cn.FESpace());
	DLi = 1.5e-10;
 	return DLi;
}

