#include "mfem.hpp"
#include <math>

// in a small header, e.g. potential_utils.hpp, or just above SolidPotential's definition
mfem::GridFuction ComputeDamb(mfem::GridFunction &Cn)
{
	mfem::GridFunction De(Cn.FESpace());
	const double D0 = 4.886868634906806e-03;
	
	for (int i = 0, i < Cn.Size(), i++){
		De(i) = D0 * exp(-7.02 - 0.83e3*Cn(i) + 0.05*(1.0e3*Cn(i)*Cn(i)));
	}
 	return De;
}

mfem::GridFuction ComputeDLi(mfem::GridFunction &Cn)
{
	mfem::GridFunction DLi(Cn.FESpace());
	for (int i = 0, i < Cn.Size(), i++){
		DLi(i) = (0.0277-0.0840*Cn(i)+0.1003*Cn(i)*Cn(i))*1.0e-8;
	}
 	return DLi;
}

